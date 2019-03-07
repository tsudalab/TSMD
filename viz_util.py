import sys
import os
import numpy as np
import pandas as pd
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
import pickle
from pylab import *
from graphviz import Graph
import mdtraj as md
from itertools import combinations
from sklearn.metrics import confusion_matrix, f1_score
# from tsmd import Node

class Node:
    def __init__(self, move = None, parent = None, state = None):
        self.parentNode = parent # "None" for the root node
        self.childNodes = []
        self.rmsd_sum = 0
        self.visits = 0
        self.state = state
        self.try_num = 0
        self.n_sim = 0
        self.J = 1




#---------------------------------------------------------------------------------------
# Utility functions
#---------------------------------------------------------------------------------------
def read_rmsd(name):
    f = open(name)
    rmsd = []
    while True:
        l = f.readline()
        if l[0] != '#' and l[0] != '@':
            break
    while l != '':
        rmsd.append(float(l.split()[1]))
        l = f.readline()
    f.close()
    return rmsd

def modify_rmsd(ip, op):
    f = open(ip)
    o = open(op, 'w')
    while True:
        l = f.readline()
        if l[0] != '#' and l[0] != '@':
            break
        else:
            o.write(l)

    cum_time = 0.0
    rmsd = l.split()[1]
    o.write(str(cum_time) + "    " + str(rmsd) + "\n")
    l = f.readline()
    cum_time += 1.0

    while l != '':
        time = float(l.split()[0])
        rmsd = l.split()[1]
        if time != 0.0:
            o.write(str(cum_time) + "    " + str(rmsd) + "\n")
            cum_time += 1.0
        l = f.readline()
    f.close()
    o.close()

#---------------------------------------------------------------------------------------
# Draw graphs
#---------------------------------------------------------------------------------------
def draw_pacs_tree_colored(log_file, out):
    log = pd.read_csv(log_file, header = None)
    log = np.array(log)
    n_cycles = log.shape[0]
    # G = Graph(format='pdf')
    G = Graph(format='svg')
    G.attr('node', shape='circle', style='filled')
    G.graph_attr.update(size="30")
    color_hex = value2hex(0) 
    G.node('0-0', '0', fillcolor=color_hex)
    # G.node('0-0', '', fillcolor=color_hex, color='white', width='12')
    color_hex = value2hex(255/n_cycles * 1) 
    for i in range(len(log[0])):
        state = '1-' + str(i)
        G.node(state, str(state), fillcolor=color_hex)
        # G.node(state, '', fillcolor=color_hex, color='white', width='12')
        G.edge('0-0', state, color='black')
    for i in range(1, len(log) + 1):
        log_i = log[i-1]
        color_hex = value2hex(255/n_cycles * (i+1))
        for j, l in enumerate(log_i):
            state = str(i + 1) + '-' + str(j)
            pstate = str(i) + '-' + str(int(l))
            G.node(state, str(state), fillcolor=color_hex)
            # G.node(state, '', fillcolor=color_hex, color='white', width='12')
            G.edge(pstate, state, color='black')
    G.render(out)
    make_colorbar(0, n_cycles * 5, out + '_colorbar')


def draw_pats_tree_colored(pkl,out, col_style='contact', **kwarg):
    with open(pkl, 'rb') as f:
        var_list = pickle.load(f)
        rootnode = var_list[0]
    if col_style == 'order':
        num_state = dfs(rootnode)
        print('num_state  ' + str(num_state))
        values = np.arange(num_state)
    elif col_style == 'rmsd':
        values = None 
    elif col_style == 'contact':
        native = kwarg['native']
        traj   = kwarg['traj']
        values = frac_native_contacts(native, traj)
    elif col_style == 'f1':
        native = kwarg['native']
        traj   = kwarg['traj']
        values = calc_f1(native, traj)
    else:
        print('Error!: you can not use such color style ' + col_style)
        exit(1)
    if values is not None:
        make_colorbar(min(values), max(values), out + '_colorbar')
    else:
        make_colorbar(0, rootnode.childNodes[0].rmsd, out + '_colorbar')
    # G = Graph(format='pdf')
    G = Graph(format='svg')
    G.attr("node", shape="circle", style="filled")
    G.graph_attr.update(size="30")
    make_graph(G, rootnode, values)
    G.render(out)


def dfs(nd):
    state = nd.state
    if len(nd.childNodes) == 0:
        return 1
    return sum([dfs(ch) for ch in nd.childNodes]) + 1    



def make_graph(G,nd, values):
    state = nd.state
    if values is not None:
        v = 255 / (max(values) - min(values)) * (values[state] - min(values))
    else: # RMSD
        if state == 0: # rootnodeのRMSDは記録されていない
            v = 255 / (0.7 - 0) * nd.childNodes[0].rmsd 
            print(nd.childNodes[0].rmsd)
        else:
            v = 255 / (0.7 - 0) * nd.rmsd 
    color_hex = value2hex(v)
    # G.node(str(state), fillcolor=color_hex, color='white', width='12')
    G.node(str(state),str(state), fillcolor=color_hex)
    parent_node = nd.parentNode
    if parent_node != None:
        parent_state = parent_node.state
        G.edge(str(parent_state), str(state), color='black')
    for child_node in nd.childNodes:
        make_graph(G,child_node, values)

def value2hex(v):
    r,g,b, _  = [int(c * 255) for c in cm.plasma_r(int(v))]
    # r,g,b, _  = [int(c * 255) for c in cm.brg(int(v))]
    r_hex,g_hex,b_hex = hex(r), hex(g), hex(b)
    color_hex = "#" + r_hex[2:].zfill(2) + g_hex[2:].zfill(2) + b_hex[2:].zfill(2)
    return color_hex

def make_colorbar(vmin, vmax, name):
   fig = plt.figure(figsize=(5,1)) 
   ax = fig.add_axes([0.05, 0.5, 0.9, 0.15])
   cmap = cm.plasma_r
   norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax)
   cb = mpl.colorbar.ColorbarBase(ax, cmap=cmap, norm=norm, orientation='horizontal')
   # cb.set_label(name) 
   plt.savefig(name + '.pdf')
   plt.close()

#---------------------------------------------------------------------------------------
# Generate a reactive trajectory
#---------------------------------------------------------------------------------------

def dfs_rmsd(nd, r_dic, n_dic):
    state = nd.state
    rmsd = nd.rmsd
    r_dic[state] = rmsd
    n_dic[state] = nd
    for ch in nd.childNodes:
      dfs_rmsd(ch,r_dic, n_dic)

#---------------------------------------------------------------------------------------
# Calculate properties of protein structure
#---------------------------------------------------------------------------------------

def best_hummer_q(traj, native):
    """Compute the fraction of native contacts according the definition from
    Best, Hummer and Eaton [1]

    Parameters
    ----------
    traj : md.Trajectory
        The trajectory to do the computation for
    native : md.Trajectory
        The 'native state'. This can be an entire trajecory, or just a single frame.
        Only the first conformation is used

    Returns
    -------
    q : np.array, shape=(len(traj),)
        The fraction of native contacts in each frame of `traj`

    References
    ----------
    ..[1] Best, Hummer, and Eaton, "Native contacts determine protein folding
          mechanisms in atomistic simulations" PNAS (2013)
    """

    BETA_CONST = 50  # 1/nm
    LAMBDA_CONST = 1.8
    NATIVE_CUTOFF = 0.45  # nanometers

    # get the indices of all of the heavy atoms
    heavy = native.topology.select_atom_indices('heavy')
    # get the pairs of heavy atoms which are farther than 3
    # residues apart
    heavy_pairs = np.array(
        [(i,j) for (i,j) in combinations(heavy, 2)
            if abs(native.topology.atom(i).residue.index - \
                   native.topology.atom(j).residue.index) > 3])

    # compute the distances between these pairs in the native state
    heavy_pairs_distances = md.compute_distances(native[0], heavy_pairs)[0]
    # and get the pairs s.t. the distance is less than NATIVE_CUTOFF
    native_contacts = heavy_pairs[heavy_pairs_distances < NATIVE_CUTOFF]
    print("Number of native contacts", len(native_contacts))

    # now compute these distances for the whole trajectory
    r = md.compute_distances(traj, native_contacts)
    # and recompute them for just the native state
    r0 = md.compute_distances(native[0], native_contacts)

    q = np.mean(1.0 / (1 + np.exp(BETA_CONST * (r - LAMBDA_CONST * r0))), axis=1)
    return q 
   
def frac_native_contacts(n, t):
    native = md.load(n)
    # t = md.load('../merged_pats.gro')
    traj = md.load(t)
    print('lenght of traj: ' + str( len(traj)))
    q = best_hummer_q(traj, native)
    return q

def calc_contacts(trj):
    C_index = [atom.index for atom in trj.topology.atoms 
               if ((atom.residue.name == 'GLY' and atom.name=='CA') or (atom.residue.name != 'GLY' and atom.name=='CB'))]
    C_comb = [(i,j) for (i,j) in combinations(C_index,2)]
    C_dist = md.compute_distances(trj, C_comb)
    contacts = C_dist[0] < 0.8
    return contacts

def calc_f1(n, t):
    native = md.load(n)
    traj = md.load(t)
    native_contacts = calc_contacts(native)
    trj_contacts    = [calc_contacts(t) for t in traj]
    f1_list =  [f1_score(native_contacts, tc) for tc in trj_contacts]
    return np.array(f1_list)

def calc_rg(t):
    traj = md.load(t)
    return md.compute_rg(traj)

