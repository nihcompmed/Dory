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


#plt_labels2 = [\
#         'FA+DSG, DdeI'\
#        ,'FA+DSG, DpnII'\
#        ,'FA+DSG, DdeI+DpnII'\
#        ]
#tuples2 = [\
#         (0,0,0)\
#        ,(0,0,1)\
#        ,(0,1,0)\
#        ,(0,1,1)\
#        ,(1,0,0)\
#        ,(1,0,1)\
#        ,(1,1,0)\
#        ,(1,1,1)\
#        ]
#index2 = pd.MultiIndex.from_tuples(tuples2, names=plt_labels2)

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
#
#s = pd.Series(np.zeros(len(tuples), dtype=int)\
#        , index=index)

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

    print('DOINT B_IDX', b_idx)

    anchor_thresh = b_idx + 1

    graph_dict = dict()

    common_anchors = dict()
    for tt in tuples:
        common_anchors[tt] = dict()

    # Initialize set intersection counts to 0
    s = pd.Series(np.zeros(len(tuples), dtype=int)\
        , index=index)

    # Create graph with node as (exp, chrom, anchor bin 1, anchor bin 2)
    for e_idx in list(range(len(labels))):
    #for e_idx in range(3):

        exp = labels[e_idx]

        e1_info = exp_count_info[exp][b_idx]
        anchor_list = e1_info['cis']

        for anchor in anchor_list:

            chrom = anchor[0]

            if chrom not in graph_dict:

                graph_dict[chrom] = nx.Graph()

            bin1 = anchor[1]
            bin2 = anchor[2]

            graph_dict[chrom].add_node((e_idx, bin1, bin2))

    for chrom in graph_dict:

        print(chrom)

        #count = 0
        #print('\n', chrom, len(graph_dict[chrom].nodes))
        for n1, n2 in itertools.combinations(list(graph_dict[chrom].nodes), 2):

            b11 = n1[1]
            b12 = n1[2]

            b21 = n2[1]
            b22 = n2[2]

            ddist = max(abs(b11-b21), abs(b12-b22))

            if ddist < anchor_thresh + 1:
                graph_dict[chrom].add_edge(n1, n2)

            #print(b11, b12, b21, b22, ddist)

            #exit()


            #count += 1
            #if count % 100000 == 0:
            #    print(count, end='\r')


        cliques = nx.find_cliques(graph_dict[chrom])

        for cl in cliques:
            #exps = []
            this_tuple = np.zeros((5,), dtype=int)
            for nn in cl:
                #exps.append(nn[0])
                this_tuple[nn[0]] = 1

            this_tuple = tuple(this_tuple)

            s[this_tuple] += 1

            # LEGACY: SAVE ANCHORS
            #for nn in cl:
            #    if chrom not in common_anchors[this_tuple]:
            #        common_anchors[this_tuple][chrom] = [[nn[1], nn[2]]]
            #    else:
            #        common_anchors[this_tuple][chrom].append([nn[1], nn[2]])

            # NEW: SAVE CLIQUES
            if chrom not in common_anchors[this_tuple]:
                common_anchors[this_tuple][chrom] = [cl]
            else:
                common_anchors[this_tuple][chrom].append(cl)

            #print(common_anchors)
            #input('w')

            #if this_tuple == (1,1,1,1,1):

            #    if chrom not in common_anchors:
            #        common_anchors[chrom] = []
            #    for nn in cl:
            #        common_anchors[chrom].append(nn[1])
            #        common_anchors[chrom].append(nn[2])

        ##if chrom in common_anchors:
        ##    common_anchors[chrom] = list(frozenset(common_anchors[chrom]))

    # LEGACY: SAVE ANCHORS
    #pickle.dump(common_anchors, open('common_cis_anchors_dim'+str(dim)+'_bidx'+str(b_idx)+'.p', 'wb'))

    # NEW: SAVE CLIQUES
    pickle.dump(common_anchors, open('common_cis_cliques_dim'+str(dim)+'_bidx'+str(b_idx)+'.p', 'wb'))

    #print(common_anchors)
    #exit()

    # log-scale for counts

    for tt in tuples:
        s[tt] = math.log(1 + s[tt])


    upsetplot(s)

    plt.ylabel('log(1+intersection size)')
    plt.savefig('cis_upset_dim'+str(dim)+'_bidx'+str(b_idx)+'.pdf', dpi=600)
    #plt.show()
    plt.cla()
    plt.clf()
    plt.close()

    continue

    #plt.cla()
    #plt.clf()
    #plt.close()



