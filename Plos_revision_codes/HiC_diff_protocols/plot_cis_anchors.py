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

#plt_labels = [\
#         'FA+DSG, DdeI'\
#        ,'FA+DSG, DpnII'\
#        ,'FA+DSG, DdeI+dpnII'\
#        ,'FA+DSG, Mnase'\
#        ,'FA, DpnII'\
#        ]
plt_labels = [\
         '(2)'\
        ,'(3)'\
        ,'(4)'\
        ,'(5)'\
        ,'(1)'\
        ]

dim = 1
print('Processing dim', dim)
exp_count_info = pickle.load(open('exp_cis_trans_anchors_dim'+str(dim)+'.p', 'rb'))


for b_idx in range(2):

    for e_idx1, e_idx2 in itertools.combinations(list(range(len(labels))), 2):

        e1_g_dist_list = []
        e2_g_dist_list = []

        exp1 = labels[e_idx1]
        exp2 = labels[e_idx2]

        plt_label1 = plt_labels[e_idx1]
        plt_label2 = plt_labels[e_idx2]

        print(exp1, exp2, b_idx)

        e1_info = exp_count_info[exp1][b_idx]
        e2_info = exp_count_info[exp2][b_idx]

        anchor1_list = e1_info['cis']
        anchor2_list = e2_info['cis']

        for aa in anchor1_list:
            a1_bin1 = aa[1]
            a1_bin2 = aa[2]
            e1_g_dist_list.append(abs(a1_bin1 - a1_bin2))

        for aa in anchor2_list:
            a1_bin1 = aa[1]
            a1_bin2 = aa[2]
            e2_g_dist_list.append(abs(a1_bin1 - a1_bin2))

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

        #print(len(e1_g_dist_list), len(minn_1), len(e2_g_dist_list), len(minn_2))

        this_pair_g_dist_list = np.array(e1_g_dist_list + e2_g_dist_list)
        this_pair_anchor_dist_list = np.array(list(minn_2) + list(minn_1))


        non_inf_idxs = np.argwhere(this_pair_anchor_dist_list != math.inf).flatten()

        this_pair_g_dist_list = this_pair_g_dist_list[non_inf_idxs]
        this_pair_anchor_dist_list = this_pair_anchor_dist_list[non_inf_idxs]

        this_pair_g_dist_list = np.log1p(this_pair_g_dist_list)
        this_pair_anchor_dist_list = np.log1p(this_pair_anchor_dist_list)

        #print(np.amax(this_pair_anchor_dist_list))

        #exit()

        #plt.scatter(this_pair_g_dist_list, this_pair_anchor_dist_list)

        plt.hist2d(np.array(this_pair_g_dist_list), np.array(this_pair_anchor_dist_list), bins=50,norm=mpl.colors.LogNorm(),  cmap='gnuplot')

        plt.title('Experiment '+plt_labels[e_idx1]+' vs. '+plt_labels[e_idx2])

        plt.colorbar()
        

        #plt.xscale('symlog', base=2)
        #plt.yscale('symlog', base=2)
        plt.xlabel('genomic distance of anchor (log (1+x))')
        plt.ylabel('closest matching anchor distance (log(1+x))')
        plt.savefig('exp'+str(e_idx1)+'_exp'+str(e_idx2)+'_histogram_closest_anchors_g_dist_bidx'+str(b_idx)+'.pdf'\
                , dpi=600)

        plt.cla()
        plt.clf()
        plt.close()
        #plt.show()
        #exit()



        minnn = list(minn_1) + list(minn_2)

        all_cis_dists = np.array(minnn)

        uni, counts = np.unique(all_cis_dists, return_counts=True)

        counts = counts/(n1+n2)*100

        cum_counts = np.cumsum(counts)

        plt.plot(uni, cum_counts, label=plt_label1+' vs. '+plt_label2, alpha=0.8)

    plt.xscale('symlog', base=2)

    #plt.yscale('symlog', base=2)

    plt.xlabel('Distance to closest anchor', fontsize=12)

    plt.ylabel('cumulative percentage', fontsize=12)

    plt.legend()

    plt.savefig('cis_anchor_distances_ALL_scaled_dim'+str(dim)+'_bidx'+str(b_idx)+'.pdf', dpi=600)
    #plt.show()
    #exit()

    plt.cla()
    plt.clf()
    plt.close()

        
