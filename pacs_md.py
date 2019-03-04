import numpy as np
import os, glob
from multiprocessing import Pool
from viz_util import *
import argparse

MIN_RMSD = 0.1

parser = argparse.ArgumentParser()
parser.add_argument('--reactant',  '-r',                   default = '0')
parser.add_argument('--target',    '-t',                   default = 'target_processed')
parser.add_argument('--topol',     '-top',                 default = 'topol')
parser.add_argument('--steps',     '-s',     type = int,   default = 1000)
parser.add_argument('--k',         '-k',     type = int,   default = 5)
parser.add_argument('--continue_', '-cn',    type = int,   default = 0)
parser.add_argument('--ntmpi',     '-ntmpi', type = int,   default = 1)
parser.add_argument('--ntomp',     '-ntomp', type = int,   default = 10) 
parser.add_argument('--delete',    '-del'  , type = int,   default = 0)
args = parser.parse_args()
reactant = args.reactant
target   = args.target
topol    = args.topol
ntmpi    = args.ntmpi
ntomp    = args.ntomp
delete   = args.delete


def first_cycle(n):
    # print(n)
    os.system('gmx grompp -f md.mdp -c %s.gro -t %s.cpt -p %s.top -o tmp_0_%d.tpr -maxwarn 5' % (reactant, reactant, topol,n))
    os.system('gmx mdrun -deffnm tmp_0_%d -ntmpi %s -ntomp %s -dlb auto' % (n, ntmpi, ntomp))
    os.system("echo 4 4 | gmx rms -s %s.gro -f tmp_0_%d.trr  -o rmsd_0_%d.xvg -tu ns" % (target, n, n))

def md_cycle(args):
    step = args[0]
    n = args[1]
    md_ind = args[2]
    last = args[3]
    input  = 'tmp_%d_%d' % (step-1, md_ind[0])
    md   = 'md_%d_%d' % (step-1, n)
    tmp     = 'tmp_%d_%d' % (step, n)
    rmsd   = 'rmsd_%d_%d' % (step, n)
    # 全ステップで計算したトラジェクトリを指定の部分まで切り取る
    os.system('echo 0 | gmx trjconv -s %s.tpr -f %s.trr -o %s.trr -e %d' % (input, input, md, md_ind[1]))
    os.system('echo 0 | gmx trjconv -s %s.tpr -f %s.trr -o %s.gro -e %d' % (input, input, md, md_ind[1]))
    # 切り取った箇所から新たにMD
    if not last:
        os.system('gmx grompp -f md.mdp -t %s.trr -o %s.tpr -c %s.gro -maxwarn 5' %(md, tmp, md))
        os.system('gmx mdrun -deffnm %s -ntmpi %d -ntomp %d -dlb auto' % (tmp, ntmpi, ntomp))
        os.system("echo 4 4 | gmx rms -s %s.gro -f %s.trr  -o %s.xvg -tu ns" % (target, tmp, rmsd))

def pacs_md(MAX_CYCLE, n_para, continue_flag):
    dt = 1 # ps for each step
    nsteps = 101 # ショートMDのステップ数+1(時刻0を含む)
    if continue_flag:
        edge_log = np.loadtxt("edge_log.csv", delimiter = ",")
        edge_log = [list(l) for l in edge_log]
        edge_log = edge_log[0:len(edge_log)-1]
        cycle_step = len(edge_log)
        cycle_num = -1
    else:
        cycle_step = 0
        o = open('log_pacs.txt','w')
        o.close()
        edge_log = []  # Treeの遷移関係のlog
        cycle_num = 0 # 現実行時におけるステップ数
    min_rmsd = 10000 # 初期値

    # while cycle_num < MAX_CYCLE and min_rmsd >= MIN_RMSD:
    while cycle_num < MAX_CYCLE:
        if not continue_flag:
            if cycle_step == 0:
                for i in range(n_para):
                    first_cycle(i)
            else:
                for i in range(n_para):
                    md_cycle([cycle_step, i, min_rmsd_idx[i],False])
        else:
            continue_flag = 0
        # 最小RMSDを記録
        rmsd_list = np.zeros(n_para * nsteps)
        for i in range(n_para):
            rmsd_list[i * nsteps : (i+1) * nsteps] = read_rmsd('rmsd_%d_%d.xvg'%(cycle_step, i))
        min_rmsd = min(rmsd_list)

        # 最初の一回だけ初期値を記録
        if cycle_step == 0:
            first_rmsd = rmsd_list[0]
            write_log(first_rmsd)
        # RMSDを記録
        if cycle_num != -1:
            write_log(min_rmsd)
        # Treeの遷移関係を記録
        min_rmsd_idx = rmsd_list.argsort()[:n_para]
        min_rmsd_idx = [(r // nsteps, r % nsteps) for r in min_rmsd_idx]
        edge_log.append([r[0] for r in min_rmsd_idx])


        # 不要なファイルを削除
        # for file in glob.glob('md_[0-9]*') + glob.glob("*#"):
        #     os.remove(file)
       
        for file in glob.glob("*#"):
            os.remove(file)
        for ext in ['trr', 'tpr', 'edr', 'log','gro', 'cpt']:
            for file in glob.glob('tmp_%d_*.%s' % (cycle_step-1, ext)):
                os.remove(file)
        if delete == 1:
            for file in glob.glob('tmp_%d_*.trr' % cycle_step-1):
                os.remove(file)

        cycle_step += 1
        cycle_num  += 1


    # 最後のMDのトラジェクリについて切り取り処理をする
    for i in range(n_para):
        md_cycle([cycle_step, i, min_rmsd_idx[i],True])
    # logをファイルに保存(make_reactive)に渡す
    np.savetxt('edge_log.csv', np.array(edge_log), delimiter=',')

    make_reactive('edge_log.csv')
    draw_pacs_tree_colored('edge_log.csv', 'pacs_tree')

# log.csvを元にshort MDのトラジェクリを繋げる
def make_reactive(edge_log):
    log = np.loadtxt(edge_log, delimiter = ",")
    back_list = [0]
    step = log.shape[0] - 1
    log_index = 0
    # 最もRMSDが小さくなったものについて、backtrackする
    for i in range(step,0, -1):
        log_index = int(log[i, log_index])
        back_list.insert(0, log_index)
    # print(back_list)
    # トラジェクトリを結合
    trajs = []
    for i, j in enumerate(back_list):
        traj = "md_%d_%d"%(i,j)
        print(traj)
        trajs.append(traj)
    traj_str = ""
    for traj in trajs:
        traj_str += (traj + " ")
    # print(traj_str)
    os.system( "gmx trjcat -f " + traj_str + " -o merged_pacs.trr -cat")
    # os.system("echo 4 4| gmx rms -s target_processed.gro -f merged_pacs.trr -tu ns -o rmsd_pacs_tmp.xvg")
    os.system("echo 4 4| gmx rms -s %s.gro -f merged_pacs.trr -tu ns -o rmsd_pacs_tmp.xvg"%target)
    modify_rmsd('rmsd_pacs_tmp.xvg', 'rmsd_pacs.xvg') # ’ダブり’があるので除去
    os.remove('rmsd_pacs_tmp.xvg')

# 各ステップでのrmsdを記録
def write_log(rmsd):
    o = open('log_pacs.txt','a')
    o.write(str(rmsd) + '\n')
    o.close()


if __name__ == '__main__':
    M = args.steps # サイクル数
    k = args.k  # 並列数
    cn = args.continue_
    pacs_md(M, k, cn)
    # for file in (glob.glob("*#") + glob.glob("md_[0-9]*") + glob.glob("res_[0-9]*") + glob.glob("rmsd_[0-9]*")):
    #     os.remove(file)
