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

reso = 5000

dim = 1

for b_idx in range(1, 2):

    cliques = pickle.load(open('common_cis_cliques_dim'+str(dim)+'_bidx'+str(b_idx)+'.p', 'rb'))

    for tt in cliques:

        #print(info[tt])
        strtt = ''
        for t in tt:
            strtt = strtt+str(t)
        print(strtt)

        #not_found_hiccup_in_all = 0
        #found_hiccup_in_all = 0

        all_g_dist = []
        all_h_dist = []

        for chrom in cliques[tt]:

            print('CHROM', chrom)

            nearest_hiccup = []

            for cl in cliques[tt][chrom]:

                #print(cl)
                ddist = dict()

                hiccup_tuple = np.array([0,0,0,0,0], dtype=int)

                this_g_dist = 0

                for anchor in cl:
                    e_idx = anchor[0]
                    b1 = anchor[1]*reso
                    b2 = anchor[2]*reso
                    if e_idx not in ddist:
                        ddist[e_idx] = math.inf

                    if 'chr'+chrom not in hiccups_data[e_idx]:
                        ddist = None
                        break

                    for hiccup in hiccups_data[e_idx]['chr'+chrom]:

                        c1 = min(hiccup[0], hiccup[1])
                        c2 = max(hiccup[0], hiccup[1])

                        val = max(abs(b1-c1), abs(b2-c2))

                        if val < ddist[e_idx]:
                            hiccup_tuple[e_idx] = 1
                            ddist[e_idx] = val
                            this_g_dist = abs(b2-b1)/reso

                #if np.amin(hiccup_tuple) == 1:
                #    found_hiccup_in_all += 1

                if ddist is not None:
                    maxx = -math.inf
                    for e_idx in ddist:
                        ddist[e_idx] = ddist[e_idx]/reso
                        maxx = max(maxx, ddist[e_idx])

                    all_g_dist.append(math.log(1+this_g_dist))
                    all_h_dist.append(math.log(1+maxx))
                #    #print(ddist, maxx)
                #else:
                #    not_found_hiccup_in_all += 1
                #    #print('NOT COMMON TO ALL')

                ##exit()

        #print(not_found_hiccup_in_all, found_hiccup_in_all)
        #plt.hist2d(all_g_dist, all_h_dist, bins=50, norm=mpl.colors.LogNorm(), cmap='gnuplot_r')
        plt.hist2d(all_g_dist, all_h_dist, bins=50, cmap='gnuplot_r', cmin=1)
        
        #plt.gca().axhline(y=math.log(1+4), color='red', ls='--', lw=2)
        #plt.text(x=0.6*plt.gca().get_xlim()[1], y=1.05*math.log(1+4), s='bin distance = 4', color='red')

        plt.xlabel('log(1+genomic distance of anchor)', fontsize=14)
        plt.ylabel('log(1+closest hiccup anchor)', fontsize=14)
        
        plt.colorbar()

        #plt.show()

        out = 'hiccups_close_cis_combination_'+strtt+'_dim'+str(dim)+'_bidx'+str(b_idx)
        plt.savefig(out+'.pdf', dpi=600)

        plt.cla()
        plt.clf()
        plt.close()

        #all_h_dist = np.array(all_h_dist)
        #uni, counts = np.unique(all_h_dist, return_counts=True)

        #cum_counts = np.cumsum(counts)

        #percentages = cum_counts/cum_counts[-1]*100

        #plt.plot(uni, percentages)
        #plt.gca().axvline(x=math.log(1+3), color='red', ls='--', lw=2)
        ##plt.text(x=0.6*plt.gca().get_xlim()[1], y=1.05*math.log(1+4), s='bin distance = 2', color='red')

        #plt.show()







