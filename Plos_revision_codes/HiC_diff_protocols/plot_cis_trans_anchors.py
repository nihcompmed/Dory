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

reso = 5000

#reso = 10000

edge_file_name = [\
         '../fa_dsg_ddel_chr1_KR_reso'+str(reso)+'.csv'\
        ,'../fa_dsg_dpn_chr1_KR_reso'+str(reso)+'.csv'\
        ,'../fa_dsn_mnase_chr1_KR_reso'+str(reso)+'.csv'\
        ,'../fa_dsg_ddel_dpn_chr1_KR_reso'+str(reso)+'.csv'\
        ,'../fa_dpn_chr1_KR_reso'+str(reso)+'.csv'\
        ]


targets = [\
         '../fa_dsg_ddel_reso'+str(reso)+'_'\
        ,'../fa_dsg_dpn_reso'+str(reso)+'_'\
        ,'../fa_dsg_mnase_reso'+str(reso)+'_'\
        ,'../fa_dsg_ddel_dpn_reso'+str(reso)+'_'\
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
         'fa+dsg, ddel'\
        ,'fa+dsg, dpn'\
        ,'fa+dsg, ddel+dpn'\
        ,'fa+dsg, mnase'\
        ,'fa, dpn'\
        ]

dim = 1
print('Processing dim', dim)
exp_count_info = pickle.load(open('exp_cis_trans_anchors_dim'+str(dim)+'.p', 'rb'))


for b_idx in range(2):

    for e_idx1, e_idx2 in itertools.combinations(list(range(len(labels))), 2):

        exp1 = labels[e_idx1]
        exp2 = labels[e_idx2]

        plt_label1 = plt_labels[e_idx1]
        plt_label2 = plt_labels[e_idx2]

        print(exp1, exp2, b_idx)

        e1_info = exp_count_info[exp1][b_idx]
        e2_info = exp_count_info[exp2][b_idx]

        anchor1_list = e1_info['cis']
        anchor2_list = e2_info['cis']

        n1 = len(anchor1_list)
        n2 = len(anchor2_list)

        if n1 == 0 or n2 == 0:
            continue

        cis_dist_mat = np.ones((n1, n2))*(math.inf)

        # cis-anchors
        for a1_idx, anchor1 in enumerate(anchor1_list):
            
            a1_chrom = anchor1[0]
            a1_bin1 = anchor1[1]
            a1_bin2 = anchor1[2]


            for a2_idx, anchor2 in enumerate(anchor2_list):

                a2_chrom = anchor2[0]
                if a1_chrom != a2_chrom:
                    continue

                a2_bin1 = anchor2[1]
                a2_bin2 = anchor2[2]

                ddist = max(abs(a1_bin1-a2_bin1), abs(a1_bin2-a2_bin2)) 

                cis_dist_mat[a1_idx, a2_idx] = ddist


        minn_1 = np.min(cis_dist_mat, axis=0)
        minn_2 = np.min(cis_dist_mat, axis=1)

        minnn = list(minn_1) + list(minn_2)

        all_cis_dists = np.array(minnn)

        uni, counts = np.unique(all_cis_dists, return_counts=True)

        counts = counts/(n1+n2)*100

        plt.plot(uni, counts, label=plt_label1+' vs. '+plt_label2, alpha=0.8)

    plt.xscale('symlog', base=2)
    plt.yscale('symlog', base=2)

    plt.xlabel('closest matching anchor')

    plt.ylabel('counts percentage')

    plt.legend()

    plt.savefig('cis_anchor_distances_ALL_scaled_dim'+str(dim)+'_bidx'+str(b_idx)+'.pdf', dpi=600)
    #plt.show()

    plt.cla()
    plt.clf()
    plt.close()

        




