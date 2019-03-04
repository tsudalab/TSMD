import numpy as np
import os, glob
from math import *
from viz_util import *
import random
import copy
import argparse
import pickle
MAX_child = 3
MAX_try = 5
MIN_RMSD = 0.1
FIRST_FLAG = 1
INF = 10000

parser = argparse.ArgumentParser()
parser.add_argument('--reactant',  '-r',                   default = '0')
parser.add_argument('--target',    '-t',                   default = 'target_processed')
parser.add_argument('--topol',     '-top',                 default = 'topol')
parser.add_argument('--steps',     '-s',     type = int,   default = 1000)
parser.add_argument('--c',         '-c',     type = float, default = 0.05)
parser.add_argument('--continue_', '-cn',    type = int,   default = 0)
parser.add_argument('--ntmpi',     '-ntmpi', type = int,   default = 1)
parser.add_argument('--ntomp',     '-ntomp', type = int,   default = 10)
parser.add_argument('--delete',    '-del'  , type = int,   default = 0)
parser.add_argument('--thresh',    '-th'  , type = float,   default = 0.1)
parser.add_argument('--alpha',    '-alp'  , type = float,   default = 1.1)
parser.add_argument('--ctype',    '-ctype' ,               default = 'normal')
args = parser.parse_args()
reactant = args.reactant
target   = args.target
c_       = args.c
topol    = args.topol
ntmpi    = args.ntmpi
ntomp    = args.ntomp
delete   = args.delete
th       = args.thresh
alpha    = args.alpha
ctype    = args.ctype

class Node:
    def __init__(self, move = None, parent = None, state = None, c = c_, depth = 0):
        self.parentNode = parent # "None" for the root node
        self.childNodes = []
        self.depth = depth
        self.rmsd_sum = 0
        self.rmsd_max = -INF
        self.visits = 0
        self.state = state
        self.untriedMoves = MAX_child
        self.c = c
        self.rmsd = INF
        self.try_num = 0
        self.n_sim = 0
        self.alpha = alpha
        self.J = 1

    def uct_select_child(self):
        if ctype == 'normal':
            s = sorted(self.childNodes, key = lambda ch: ch.rmsd_max + self.c * sqrt(2*log(self.visits)/ch.visits))[-1] 
        elif ctype == 'adaptive':
            child_rmsds = [ch.rmsd_max for ch in self.childNodes]
            rmsd_diff = max(child_rmsds) - min(child_rmsds)
            c_adap = rmsd_diff * self.c + 0.0001
            s = sorted(self.childNodes, key = lambda ch: ch.rmsd_max + c_adap * sqrt(2*log(self.visits)/ch.visits))[-1] 
        elif ctype == "adaptive2":
            child_rmsds = [ch.rmsd_max for ch in self.childNodes]
            rmsd_diff = max(child_rmsds) - min(child_rmsds)
            c_adap = np.sqrt(2)*self.J/4 * rmsd_diff
            s = sorted(self.childNodes, key = lambda ch: ch.rmsd_max + c_adap * sqrt(2*log(self.visits)/ch.visits))[-1] 
        return s

    def calc_uct(self):
        pnd = self.parentNode
        if pnd == None:
            return -1
        if ctype == "normal":
            uct = self.rmsd_max + self.c * sqrt(2*log(pnd.visits) / self.visits)
        elif ctype == "adaptive":
            child_rmsds = [ch.rmsd_max for ch in pnd.childNodes]
            rmsd_diff = max(child_rmsds) - min(child_rmsds)
            c_adap = rmsd_diff * self.c + 0.0001
            uct = self.rmsd_max + c_adap * sqrt(2*log(pnd.visits) / self.visits)
        elif ctype == "adaptive2":
            child_rmsds = [ch.rmsd_max for ch in pnd.childNodes]
            rmsd_diff = max(child_rmsds) - min(child_rmsds)
            c_adap = np.sqrt(2)*pnd.J/4 * rmsd_diff
            uct = self.rmsd_max + c_adap * sqrt(2*log(pnd.visits) / self.visits)
        return uct

    def make_child(self, s, d):
        n = Node(parent = self, state = s, c = self.c, depth = d)
        return n

    def add_child(self, n):
        self.untriedMoves -= 1
        self.childNodes.append(n)
        return n

    def update(self, result, similarity_list, dec_flag):
        self.visits += 1
        if not dec_flag:
            self.J += 0.1
        # if result > self.rmsd_max:
        #     self.rmsd_max = result
        # if self.state != 0 and self.state-1 < len(similarity_list):
        #     self.n_sim = similarity_list[self.state - 1]

    def MDrun(self):
        global FIRST_FLAG
        state = self.state
        pstate = self.parentNode.state
        self.parentNode.try_num += 1
        tmp = str(state) + '_tmp'
        if pstate == 0:
            os.system('gmx grompp -f md.mdp -c %s.gro -t %s.cpt -p %s.top -o %s.tpr -maxwarn 5' % (reactant, reactant, topol, tmp))
        else:
            os.system('gmx grompp -f md.mdp -t md_%d.trr -o %s.tpr -c md_%d.gro -maxwarn 5' %(pstate, tmp, pstate))
        os.system('gmx mdrun -deffnm %s  -ntmpi %d  -ntomp %d -dlb auto' % (tmp,  ntmpi, ntomp))

        os.system("echo 4 4 | gmx rms -s %s.gro -f %s.trr  -o rmsd_%d.xvg -tu ns" % (target, tmp, state)) # rmsdを測定
        rmsds = np.array(read_rmsd('rmsd_%d.xvg'%state))
        # 初期RMSDを書き込み
        if FIRST_FLAG:
            first_rmsd = rmsds[0]
            o = open('log_pats.txt','w')
            o.write(str(first_rmsd) + '\n')
            o.close()
            FIRST_FLAG = 0

        min_rmsd = np.min(rmsds)
        min_i = rmsds.argsort()[0]
        os.system('echo 0 | gmx trjconv -s %s.tpr -f %s.trr -o md_%s.trr -e %d' % (tmp, tmp, state, min_i)) # 最小値までのトラジェクトリーを切り出し
        os.system('echo 0 | gmx trjconv -s %s.tpr -f %s.trr -o md_%s.gro -e %d -b %d' % (tmp, tmp, state, min_i, min_i))
        os.system('echo 4 | gmx trjconv -s %s.tpr -f %s.trr -o md_bb_%s.gro -e %d -b %d' % (tmp, tmp, state, min_i, min_i))

        for file in glob.glob("*#"):
            os.remove(file)
        for ext in ['trr', 'tpr', 'edr', 'log','gro', 'cpt']:
            for file in glob.glob('%s.%s' % (tmp,ext)):
                os.remove(file)
        if delete == 1:
            for file in glob.glob('%s*' % tmp):
                os.remove(file)
        self.rmsd = min_rmsd
        return min_rmsd

    def prog_widenning(self):
        return sqrt(self.visits) < (3/2) * len(self.childNodes)

def update_rmsd_max(node, similarity_list):
    state = node.state
    if state != 0:
        node.n_sim = similarity_list[state-1]
    if len(node.childNodes) == 0:
        rmsd_cor = node.rmsd * (node.alpha ** node.n_sim)
        node.rmsd_max = -rmsd_cor
        return node.rmsd_max
    rmsd_cors = [update_rmsd_max(ch, similarity_list) for ch in node.childNodes]
    rmsd_max_cor = max(rmsd_cors)
    node.rmsd_max = rmsd_max_cor
    return node.rmsd_max


# 全ての構造との距離を計算し、θ nm未満のものをカウント
def check_similarity(nd, similarity_list, dec_flag):
    if nd.state <= 1:
        return [0]
    os.system('echo 4 4 |gmx rms -s md_bb_%s.gro -f all_structure.gro -o rmsd_tmp.xvg'%(nd.state))
    rmsd_tmp = np.array(read_rmsd('rmsd_tmp.xvg'))
    bool_tmp = rmsd_tmp < th
    n_similar = np.sum(bool_tmp)
    for idx in np.where(bool_tmp)[0]:
        similarity_list[idx] += 1
    if dec_flag:    
      similarity_list.append(n_similar)
    for file in glob.glob("*#"):
        os.remove(file)
    return similarity_list

def make_reactive(pkl):
  with open(pkl, 'rb') as f:
      l = pickle.load(f)

  rootnode = l[0]
  rmsd_dic = {}
  node_dic = {}
  dfs_rmsd(rootnode, rmsd_dic, node_dic)
  max_state = min(rmsd_dic, key=rmsd_dic.get)
  max_node = node_dic[max_state]

  trjs = ""
  node = max_node
  while True:
      trjs = "md_" + str(node.state) + ".trr " + trjs
      node = node.parentNode
      if node.parentNode == None:
          break
  print(trjs)
  os.system('gmx trjcat -f {0} -o reactive.trr -cat'.format(trjs))


def UCT(rootstate):
    steps     = args.steps
    c_        = args.c
    cn        = args.continue_

    succeed = 0
    if cn:
        with open('vars.pickle','rb') as f:
            var_list = pickle.load(f)
        rootnode = var_list[0]
        n_state = var_list[1]
        best_rmsd = var_list[2]
        max_node = var_list[3]
        similarity_list = var_list[4]
        global FIRST_FLAG
        FIRST_FLAG = 0
    else:
        os.system('rm all_structure.gro')
        rootnode = Node(state = rootstate, c = c_)
        n_state = rootstate
        best_rmsd = INF
        max_node    = rootnode
        similarity_list = []
        o = open('log_pats.txt','w')
        o.close()
    for i in range(steps):
        o = open('log_pats.txt','a')
        node = rootnode
        state = rootstate
        # Select
        # while (node.untriedMoves == 0 or node.prog_widenning()) and node.childNodes != []: # progressive widennning
        while (node.untriedMoves == 0 or node.try_num >= MAX_try) and node.childNodes != []: # node is fully expanded and non-terminal
            # node = node.uct_select_childByDepth()
            node = node.uct_select_child()
            state = node.state

        # Expand
        parent_node =node
        parent_rmsd = node.rmsd
        parent_depth = node.depth
        state = n_state + 1
        depth = parent_depth + 1
        # node = node.add_child(state) # add child and descend tree
        node = node.make_child(s = state, d = depth)
        min_rmsd = node.MDrun()
        dec_flag = parent_rmsd - min_rmsd > 0.0001
        similarity_list = check_similarity(node, similarity_list, dec_flag)
        if dec_flag: # RMSDが減少した場合のみexpandする
            parent_node.add_child(node)
            os.system('cat md_bb_%s.gro >> all_structure.gro'%state) # 構造を保存
            n_state += 1
        # Backpropagate
        result = -1 * min_rmsd
        if min_rmsd < best_rmsd:
            best_rmsd = min_rmsd
            max_node = node
        while node != None:
            node.update(result, similarity_list, dec_flag)
            node = node.parentNode

        o.write(str(best_rmsd) + '\n')
        o.close()
        # 全ノードのrmsd_maxをupdate
        update_rmsd_max(rootnode, similarity_list)
        if best_rmsd < MIN_RMSD:
            succeed = 1
            # break

    # 途中経過をpickleに保存
    var_list = []
    var_list.append(rootnode)
    var_list.append(n_state)
    var_list.append(best_rmsd)
    var_list.append(max_node)
    var_list.append(similarity_list)
    with open('vars.pickle', mode = 'wb') as f:
        pickle.dump(var_list, f)

    # make reactive trajectory and tree graph
    make_reactive('vars.pickle')
    draw_pats_tree_colored('vars.pickle', 'tree_graph', col_style='order')


if __name__ == "__main__":
    UCT(0)
