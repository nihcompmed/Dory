import pydory as dory
import check_connected_helper as cch
import helper_functions as hf 
import helper_funcs_v2 as hf2
import networkx as nx
import matplotlib.pyplot as plt
import os
import numpy as np
import itertools
import seaborn as sns
from sklearn.neighbors import KernelDensity
import scipy
import matplotlib as mpl
import cooler 
import pickle
import math
import networkx as nx
import pandas as pd
from upsetplot import plot as upsetplot


reso = 5000

#reso = 10000

edge_file_name = [\
         '../fa_dsg_ddel_chr1_KR_reso'+str(reso)+'.csv'\
        ,'../fa_dsg_dpn_chr1_KR_reso'+str(reso)+'.csv'\
        ,'../fa_dsg_ddel_dpn_chr1_KR_reso'+str(reso)+'.csv'\
        ,'../fa_dsn_mnase_chr1_KR_reso'+str(reso)+'.csv'\
        ,'../fa_dpn_chr1_KR_reso'+str(reso)+'.csv'\
        ]


targets = [\
         '../fa_dsg_ddel_reso'+str(reso)+'_'\
        ,'../fa_dsg_dpn_reso'+str(reso)+'_'\
        ,'../fa_dsg_ddel_dpn_reso'+str(reso)+'_'\
        ,'../fa_dsg_mnase_reso'+str(reso)+'_'\
        ,'../fa_dpn_reso'+str(reso)+'_'\
        ]

labels = [\
         'fa_dsg_ddel_5kb'\
        ,'fa_dsg_dpn_5kb'\
        ,'fa_dsg_ddel_dpn_5kb'\
        ,'fa_dsg_mnase_5kb'\
        ,'fa_dpn_5kb'\
        ]

plt_labels = [\
         'FA+DSG, DdeI'\
        ,'FA+DSG, DpnII'\
        ,'FA+DSG, DdeI+DpnII'\
        ,'FA+DSG, Mnase'\
        ,'FA, Dpn'\
        ]

dim = 1
print('Processing dim', dim)
exp_count_info = pickle.load(open('exp_cis_trans_anchors_dim'+str(dim)+'.p', 'rb'))

# Initialize Pandas multi index
tuples = [\
           (0,0,0,0,0)\
          ,(0,0,0,0,1)\
          ,(0,0,0,1,0)\
          ,(0,0,0,1,1)\
          ,(0,0,1,0,0)\
          ,(0,0,1,0,1)\
          ,(0,0,1,1,0)\
          ,(0,0,1,1,1)\
          ,(0,1,0,0,0)\
          ,(0,1,0,0,1)\
          ,(0,1,0,1,0)\
          ,(0,1,0,1,1)\
          ,(0,1,1,0,0)\
          ,(0,1,1,0,1)\
          ,(0,1,1,1,0)\
          ,(0,1,1,1,1)\
          ,(1,0,0,0,0)\
          ,(1,0,0,0,1)\
          ,(1,0,0,1,0)\
          ,(1,0,0,1,1)\
          ,(1,0,1,0,0)\
          ,(1,0,1,0,1)\
          ,(1,0,1,1,0)\
          ,(1,0,1,1,1)\
          ,(1,1,0,0,0)\
          ,(1,1,0,0,1)\
          ,(1,1,0,1,0)\
          ,(1,1,0,1,1)\
          ,(1,1,1,0,0)\
          ,(1,1,1,0,1)\
          ,(1,1,1,1,0)\
          ,(1,1,1,1,1)\
        ]

index = pd.MultiIndex.from_tuples(tuples, names=plt_labels)


## FOR ONLY FIRST 3 EXPERIMENTS
#
#tuples = [\
#         (0,0,0)\
#        ,(0,0,1)\
#        ,(0,1,0)\
#        ,(0,1,1)\
#        ,(1,0,0)\
#        ,(1,0,1)\
#        ,(1,1,0)\
#        ,(1,1,1)\
#        ]
#index = pd.MultiIndex.from_tuples(tuples, names=plt_labels[:3])

s = pd.Series(np.zeros(len(tuples), dtype=int)\
        , index=index)

#print(index)

#s = pd.Series([1,2,3,4,5,6,7,8]\
#        , index=index)

#print(s[(0,0,0)])
#s[(0,1,0)] = 100
#print(s)
#exit()

#plot(s)
#plt.show()


for b_idx in range(2):

    s = pd.Series(np.zeros(len(tuples), dtype=int)\
        , index=index)

    common_anchors = dict()
    for tt in tuples:
        common_anchors[tt] = []

    # Create graph with node as (exp, chrom, anchor bin 1, anchor bin 2)
    G = nx.Graph()

    for e_idx in list(range(len(labels))):
    #for e_idx in range(3):

        exp = labels[e_idx]

        e1_info = exp_count_info[exp][b_idx]
        anchor_list = e1_info['trans']

        for anchor in anchor_list:

            chrom1 = anchor[0]
            bin1 = anchor[1]

            chrom2 = anchor[2]
            bin2 = anchor[3]

            G.add_node((e_idx, chrom1, bin1, chrom2, bin2))

    print(len(list(G.nodes)))

    for n1, n2 in itertools.combinations(list(G.nodes), 2):

        c11 = n1[1]
        c12 = n1[3]

        c21 = n2[1]
        c22 = n2[3]

        if c11 != c21 and c11 != c22:
            continue

        if c12 != c21 and c12 != c22:
            continue

        b11 = n1[2]
        b12 = n1[4]

        b21 = n2[2]
        b22 = n2[4]

        if c11 == c21 and c12 == c22:
            ddist = max(abs(b11-b21), abs(b12-b22))
        elif c11 == c22 and c12 == c21:
            ddist = max(abs(b11-b22), abs(b12-b21))

        if ddist < 3:
            G.add_edge(n1, n2)


    cliques = nx.find_cliques(G)

    for cl in cliques:

        this_tuple = np.zeros((5,), dtype=int)

        ## FOR ONLY FIRST 3
        #this_tuple = np.zeros((3,), dtype=int)

        for nn in cl:
            this_tuple[nn[0]] = 1

        this_tuple = tuple(this_tuple)
        s[this_tuple] += 1

        for nn in cl:
            common_anchors[this_tuple].append([nn[1], nn[2], nn[3], nn[4]])

    for tt in tuples:
        s[tt] = math.log(1 + s[tt])

    upsetplot(s)
    
    plt.ylabel('log(1+intersection size)')
    plt.savefig('trans_upset_dim'+str(dim)+'_bidx'+str(b_idx)+'.pdf', dpi=600)
    #plt.show()

    plt.cla()
    plt.clf()
    plt.close()

    pickle.dump(common_anchors, open('common_trans_anchors_dim'+str(dim)+'_bidx'+str(b_idx)+'.p', 'wb'))

    #continue

    #plt.savefig('trans_upset_dim'+str(dim)+'_bidx'+str(b_idx)+'.pdf', dpi=600)


    ##plt.show()
    #plt.cla()
    #plt.clf()
    #plt.close()





#    for chrom in graph_dict:
#
#        #count = 0
#        #print('\n', chrom, len(graph_dict[chrom].nodes))
#        for n1, n2 in itertools.combinations(list(graph_dict[chrom].nodes), 2):
#
#            b11 = n1[1]
#            b12 = n1[2]
#
#            b21 = n2[1]
#            b22 = n2[2]
#
#            ddist = max(abs(b11-b21), abs(b12-b22))
#
#            if ddist < 3:
#                graph_dict[chrom].add_edge(n1, n2)
#
#            #print(b11, b12, b21, b22, ddist)
#
#            #exit()
#
#
#            #count += 1
#            #if count % 100000 == 0:
#            #    print(count, end='\r')
#
#        for comp in nx.connected_components(graph_dict[chrom]):
#
#            #this_tuple = np.zeros((5,), dtype=int)
#            # FOR ONLY FIRST 3
#            this_tuple = np.zeros((3,), dtype=int)
#
#            for node in comp:
#                this_tuple[node[0]] = 1
#
#            this_tuple = tuple(this_tuple)
#            s[this_tuple] += 1
#
#    upsetplot(s)
#    #plt.show()
#    plt.savefig('trans_upset_dim'+str(dim)+'_bidx'+str(b_idx)+'.pdf', dpi=600)
#
#    plt.cla()
#    plt.clf()
#    plt.close()
#
#
#
