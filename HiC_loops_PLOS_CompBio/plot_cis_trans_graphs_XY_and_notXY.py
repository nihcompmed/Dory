import numpy as np
import cooler
import h5py
import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn as sns
from numba import njit
import pickle
import networkx as nx
import itertools as it
import math


def get_chrom(binn, chrom_start_bins_list):

    for chrom, start in chrom_start_bins_list:

        if start >= binn:
            break
        else:
            this_chrom = chrom

    return this_chrom

def get_info(cfile_name, exp_name, dim, G_XY, G_notXY, factor):

    cfile = exp_name+'.cool'
    
    cf = cooler.Cooler(cfile)
    
    chrom_start_bins = {}
    chrom_start_bins_list = []

    # Initialize graph

    for ch in cf.chromnames:
        beg, end = cf.extent(ch)
        chrom_start_bins[ch] = beg
        chrom_start_bins_list.append([ch, beg])
    
    chrom_start_bins_list.append([ 'last_bin', cf.extent(cf.chromnames[-1])[1] ])
    #print(chrom_start_bins_list)
    #exit()
    
    # Keep it simple
    # Initiate a chorm dict (except last bin)
    chrom_dict = dict()
    for [chrom, binn] in chrom_start_bins_list[:-1]:
        chrom_dict[chrom] = {'cis':dict(), 'trans':dict()}
    
    # Go over smoothened cycles
    smoothened = exp_name + '/smooth_H1_estimated_pers.txt' 
    
    
    ff = open(smoothened, 'r')
    
    
    for line in ff:
        if line[0] != 's':
            continue
        line = line.split(',')
        line = line[1:-1]
        line = [int(x) for x in line]

        # Smoothened loops are stored differently
    
        # Ignore degenerate cycles
        if len(line) < 3:
            continue
    
        nn = 0
    
        chroms_in_cycle = frozenset()
    
        while nn < len(line):
            bin1 = line[nn] 
            bin2 = line[(nn+1)%len(line)] 
    
            bin1_chrom = get_chrom(bin1, chrom_start_bins_list)
            bin2_chrom = get_chrom(bin2, chrom_start_bins_list)
    
            chroms_in_cycle = chroms_in_cycle.union(frozenset({bin1_chrom, bin2_chrom}))
    
            nn += 1
    

        if len(chroms_in_cycle) == 1:
            continue
        else:
            chroms = list(chroms_in_cycle)
            if 'chrX' in chroms or 'chrY' in chroms:
                for ch1, ch2 in it.combinations(chroms, 2):
                    if G_XY.has_edge(ch1, ch2):
                        G_XY[ch1][ch2]['weight'] += factor
                    else:
                        G_XY.add_edge(ch1, ch2, weight=factor)
            else:
                for ch1, ch2 in it.combinations(chroms, 2):
                    if G_notXY.has_edge(ch1, ch2):
                        G_notXY[ch1][ch2]['weight'] += factor
                    else:
                        G_notXY.add_edge(ch1, ch2, weight=factor)
        
    return G_XY, G_notXY



chrom_list = [\
'chr1',\
'chr2',\
'chr3',\
'chr4',\
'chr5',\
'chr6',\
'chr7',\
'chr8',\
'chr9',\
'chr10',\
'chr11',\
'chr12',\
'chr13',\
'chr14',\
'chr15',\
'chr16',\
'chr17',\
'chr18',\
'chr19',\
'chr20',\
'chr21',\
'chr22',\
'chrX',\
'chrY',\
]


dim = 1

G_XY = nx.Graph()
G_notXY = nx.Graph()

# control
cfile_name = '4DNFITH9KJF1'
exp_name = '4DNFITH9KJF1'
factor = 1
G_XY, G_notXY = get_info(cfile_name, exp_name, dim, G_XY, G_notXY, factor)

# with auxin
cfile_name = '4DNFI2X45Z5L'
exp_name = '4DNFI2X45Z5L'
factor = -1
G_XY, G_notXY = get_info(cfile_name, exp_name, dim, G_XY, G_notXY, factor)



edge_weight = list(np.array(list(nx.get_edge_attributes(G_XY,'weight').values())))

print(edge_weight)

#exit()

H = nx.Graph()

max_weight = 0

for u, v in G_XY.edges:
    this_weight = math.log2(abs(G_XY[u][v]['weight'])+1)
    H.add_edge(u, v, weight=this_weight)
    max_weight = max(this_weight, max_weight)
    #alpha=math.log2(auxin_G[u][v]['weight'])/max_edge_weight


pos = nx.spring_layout(H, seed=7)

for u, v in G_XY.edges:

    this_weight = G_XY[u][v]['weight']

    if this_weight > 0:
        color = '#005AB5'
    else:
        color = '#DC3220'

    this_weight = math.log2(abs(G_XY[u][v]['weight'])+1)

    nx.draw_networkx_edges(
        G_XY,
        pos,
        edgelist=[(u, v)],
        width=2.0,
        alpha=this_weight/max_weight,
        edge_color=color,
    )



labels = {}

for idx, node in enumerate(list(G_XY.nodes)):
    labels[node] = node[3:]

nx.draw_networkx_nodes(G_XY, pos, nodelist=list(G_XY.nodes)\
        , node_size = 600, node_color="black", alpha=0.75)


nx.draw_networkx_labels(G_XY, pos, labels, font_size=16, font_color='white')

plt.gca().axis('off')

plt.savefig('figures/HiC_trans_graph_diff_XY.pdf', dpi=600)
#plt.show()

plt.cla()
plt.clf()

# Plot notXY graph

max_weight = 0

for u, v in G_notXY.edges:
    this_weight = math.log2(abs(G_notXY[u][v]['weight'])+1)
    max_weight = max(this_weight, max_weight)


for u, v in G_notXY.edges:

    this_weight = G_notXY[u][v]['weight']

    if this_weight > 0:
        color = '#005AB5'
    else:
        color = '#DC3220'

    this_weight = math.log2(abs(G_notXY[u][v]['weight'])+1)

    nx.draw_networkx_edges(
        G_notXY,
        pos,
        edgelist=[(u, v)],
        width=2.0,
        alpha=this_weight/max_weight,
        edge_color=color,
    )

labels = {}

for idx, node in enumerate(list(G_notXY.nodes)):
    labels[node] = node[3:]

nx.draw_networkx_nodes(G_notXY, pos, nodelist=list(G_notXY.nodes)\
        , node_size = 600, node_color="black", alpha=0.75)

nx.draw_networkx_labels(G_notXY, pos, labels, font_size=16, font_color='white')

plt.gca().axis('off')

plt.savefig('figures/HiC_trans_graph_diff_notXY.pdf', dpi=600)

