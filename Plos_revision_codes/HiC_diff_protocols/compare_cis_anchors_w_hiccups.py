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

hiccup_files = [\
         '../fa_dsg_ddel_hiccups.bedpe'\
        ,'../fa_dsg_dpn_hiccups.bedpe'\
        ,'../fa_dsg_ddel_dpn_hiccups.bedpe'\
        ,'../fa_dsg_mnase_hiccups.bedpe'\
        ,'../fa_dpn_hiccups.bedpe'\
                ]


hiccups_data = []
for ff in hiccup_files:

    # Result is in base pairs
    # dict[chrom] = [cent1, cent2\
    # , start1, end1\
    # , start2, end2\
    # , c_size, c_count\
    # , c_donut, c_vertical\
    # , c_horizontal, c_lowleft]

    hiccups_data.append(hf2.get_hiccup_peaks(ff))


dim = 2

exp_count_info = pickle.load(open('exp_cis_trans_anchors_dim'+str(dim)+'.p', 'rb'))

reso = 5000
print(f'Doing dim {dim} reso {reso}')

for b_idx in range(2):

    for e_idx in list(range(len(labels))):

        exp1 = labels[e_idx]
        plt_label1 = plt_labels[e_idx]

        print(exp1, b_idx)

        e1_info = exp_count_info[exp1][b_idx]
        anchor1_list = e1_info['cis']

        all_matching = []


        # cis-anchors
        for a1_idx, anchor1 in enumerate(anchor1_list):
            
            a1_chrom = anchor1[0]

            # Scale to base pairs
            a1_bin1 = anchor1[1]*reso
            a1_bin2 = anchor1[2]*reso

            this_gdist = a1_bin2 - a1_bin1

            ddist = math.inf
            entry = None

            if 'chr'+a1_chrom not in hiccups_data[e_idx]:
                all_matching.append([this_gdist, math.inf])
                continue

            # Go through hiccup peaks and find closest anchor
            for hiccup in hiccups_data[e_idx]['chr'+a1_chrom]:

                c1 = min(hiccup[0], hiccup[1])
                c2 = max(hiccup[0], hiccup[1])

                val = max(abs(a1_bin1-c1), abs(a1_bin2-c2))

                if val < ddist:
                    ddist = val
                    entry = hiccup

            # Scale back to 5kb resolution
            all_matching.append([this_gdist, ddist] + entry)

        #print(all_matching)
        #exit()

        x_valls = []
        y_valls = []

        skipping = 0

        for matching in all_matching:
            if matching[1] != math.inf:
                # Scale to reso (5kb)
                x_valls.append(matching[0]/reso)
                y_valls.append(matching[1]/reso)
                # Scale to 5bk reso and take log2
                #x_valls.append(math.log2(matching[0]/reso))
            #else:
                #x_valls.append(matching[0])
            else:
                skipping += 1

        if len(x_valls) == 0:
            continue

        print(f'Not found {skipping} peaks')

        # log scale
        x_valls = np.log1p(np.array(x_valls))
        y_valls = np.log1p(np.array(y_valls))
        #plt.scatter(x_valls, y_valls, alpha=0.7, color='blue')
        plt.hist2d(x_valls, y_valls, bins=50, norm=mpl.colors.LogNorm(), cmap='gnuplot_r')
        plt.colorbar()

        plt.xlabel('genomic distance of cis-anchor (log 1+5kb)')
        plt.ylabel('nearest hiccup centroid (log 1+5kb)')

        plt.savefig('compare_hiccups_exp'+plt_labels[e_idx]\
                    +'_dim'+str(dim)\
                    +'_bidx'+str(b_idx)\
                    +'.pdf'\
                    , dpi=600)

        

        plt.cla()
        plt.clf()
        plt.close()
        
        #plt.show()
    
    #plt.hist(valls, bins=100)
    #plt.show()
    #exit()





#    for e_idx1, e_idx2 in itertools.combinations(list(range(len(labels))), 2):
#
#        exp1 = labels[e_idx1]
#        exp2 = labels[e_idx2]
#
#        plt_label1 = plt_labels[e_idx1]
#        plt_label2 = plt_labels[e_idx2]
#
#        print(exp1, exp2, b_idx)
#
#        e1_info = exp_count_info[exp1][b_idx]
#        e2_info = exp_count_info[exp2][b_idx]
#
#        anchor1_list = e1_info['cis']
#        anchor2_list = e2_info['cis']
#
#        n1 = len(anchor1_list)
#        n2 = len(anchor2_list)
#
#        cis_dist_mat = np.ones((n1, n2))*(math.inf)
#
#        # cis-anchors
#        for a1_idx, anchor1 in enumerate(anchor1_list):
#            
#            a1_chrom = anchor1[0]
#            a1_bin1 = anchor1[1]
#            a1_bin2 = anchor1[2]
#
#
#            for a2_idx, anchor2 in enumerate(anchor2_list):
#
#                a2_chrom = anchor2[0]
#                if a1_chrom != a2_chrom:
#                    continue
#
#                a2_bin1 = anchor2[1]
#                a2_bin2 = anchor2[2]
#
#                ddist = max(abs(a1_bin1-a2_bin1), abs(a1_bin2-a2_bin2)) 
#
#                cis_dist_mat[a1_idx, a2_idx] = ddist
#
#
#        minn_1 = np.min(cis_dist_mat, axis=0)
#        minn_2 = np.min(cis_dist_mat, axis=1)
#
#        minnn = list(minn_1) + list(minn_2)
#
#        all_cis_dists = np.array(minnn)
#
#        uni, counts = np.unique(all_cis_dists, return_counts=True)
#
#        counts = counts/(n1+n2)*100
#
#        plt.plot(uni, counts, label=plt_label1+' vs. '+plt_label2, alpha=0.8)
#
#    plt.xscale('symlog', base=2)
#    plt.yscale('symlog', base=2)
#
#    plt.xlabel('closest matching anchor')
#
#    plt.ylabel('counts percentage')
#
#    plt.legend()
#
#    plt.savefig('cis_anchor_distances_ALL_scaled_bidx_'+str(b_idx)+'.pdf', dpi=600)
#    #plt.show()
#
#    plt.cla()
#    plt.clf()
#    plt.close()
#
#        
#
#
#
#
#
#
