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

    info = pickle.load(open('common_cis_anchors_dim'+str(dim)+'_bidx'+str(b_idx)+'.p', 'rb'))

    all_matching = dict()

    # First get nearest ATAC peak
    for chrom in info:
        print('DOING CHROM', chrom)
        anchors = info[chrom]

        nn = 0
        while nn < len(anchors)-1:

            bin1 = anchors[nn]
            bin2 = anchors[nn+1]

            # Scale to bp
            bin1 = bin1*reso
            bin2 = bin2*reso
            
            nn += 2

            this_gdist = abs(bin2 - bin1)

            for e_idx in range(3):

                if e_idx not in all_matching:
                    all_matching[e_idx] = []

                if 'chr'+chrom not in hiccups_data[e_idx]:
                    all_matching[e_idx].append([this_gdist, math.inf])
                    continue

                ddist = math.inf

                # Go through hiccup peaks and find closest anchor
                for hiccup in hiccups_data[e_idx]['chr'+chrom]:

                    c1 = min(hiccup[0], hiccup[1])
                    c2 = max(hiccup[0], hiccup[1])

                    val = max(abs(bin1-c1), abs(bin2-c2))

                    if val < ddist:
                        ddist = val
                        entry = hiccup

                # Scale back to 5kb resolution
                all_matching[e_idx].append([this_gdist, ddist] + entry)

                #print(all_matching)
                #exit()


    for e_idx in range(3):

        x_valls = []
        y_valls = []

        skipping = 0

        for matching in all_matching[e_idx]:
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
        
        plt.savefig('compare_common_anchor_w_hiccups_exp'+plt_labels[e_idx]\
                    +'_dim'+str(dim)\
                    +'_bidx'+str(b_idx)\
                    +'.pdf'\
                    , dpi=600)
        
        
        
        plt.cla()
        plt.clf()
        plt.close()
    
        
