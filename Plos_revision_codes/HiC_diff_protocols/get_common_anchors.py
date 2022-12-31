import pandas as pd
import pickle
import matplotlib.pyplot as plt
import numpy as np
import pandas
import math
import h5py
from numba import njit, prange
from numba import njit
import pickle
import cooler 
import seaborn as sns
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
import helper_funcs_v2 as hf2
import matplotlib.cm as cm
from upsetplot import plot as upsetplot
import networkx as nx
import itertools
import matplotlib as mpl


# Resolution 5kb first
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

plt_labels_arr = np.array(plt_labels)

exps = [\
         'fa_dsg_ddel_5kb.cool'\
        ,'fa_dsg_dpn_5kb.cool'\
        ,'fa_dsg_ddel_dpn_5kb.cool'\
        ,'fa_dsg_mnase_5kb.cool'\
        ,'fa_dpn_5kb.cool'\
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

labels_for_intersection = [0, 1, 2, 3, 4]

reso = 5000

#anchor_thresh = 2

#labels_for_intersection = [0, 1, 2]
log_scale = 1

tuples = hf2.get_binary_tuples_combination(len(labels_for_intersection))
index = pd.MultiIndex.from_tuples(tuples, names=plt_labels_arr[labels_for_intersection])

dim = 1
print('Processing dim', dim)
exp_count_info = pickle.load(open('exp_cis_trans_anchors_dim'+str(dim)+'.p', 'rb'))


nbins = 50

all_common_anchors = dict()

for b_idx in range(2):

    all_common_anchors[b_idx] = dict()

    # ANCHOR THRESH IS 1 if birth idx is for gene distance 1
    # ANCHOR THRESH IS 2 if birth idx is for gene distance 2
    anchor_thresh = b_idx + 1

    # FIRST FIND ALL MATCHING HICCUPS FOR THE DIFFERENT EXPERIMENTS
    all_matching = dict()
    all_matching_common = dict()
    for h_e_idx in range(len(hiccup_files)):
        all_matching[h_e_idx] = []
        all_matching_common[h_e_idx]= []

    minn_g_dist = anchor_thresh
    minn_h_dist = 0

    maxx_g_dist = -math.inf
    maxx_h_dist = -math.inf

    
    for e_idx in range(len(labels)):

        print(labels[e_idx])

        exp = labels[e_idx]

        e1_info = exp_count_info[exp][b_idx]

        anchor_list = e1_info['cis']
         
        for anchor in anchor_list:

            chrom = anchor[0]

            bin1 = anchor[1]
            bin2 = anchor[2]

            this_g_dist = abs(bin2 - bin1)

            minnn = math.inf
            entry = []

            a1_bin1 = bin1*reso
            a1_bin2 = bin2*reso

            if 'chr'+chrom not in hiccups_data[e_idx]:
                continue

            for hiccup in hiccups_data[e_idx]['chr'+chrom]:

                c1 = min(hiccup[0], hiccup[1])
                c2 = max(hiccup[0], hiccup[1])

                val = max(abs(a1_bin1-c1), abs(a1_bin2-c2))

                if val < minnn:
                    minnn = val
                    entry = hiccup

            all_matching[e_idx].append([this_g_dist, minnn] + entry)

            maxx_g_dist = max(maxx_g_dist, this_g_dist)
            maxx_h_dist = max(maxx_h_dist, minnn)

    maxx_h_dist = maxx_h_dist/reso

    # Binn uniformly in log-scale
    minn_g_dist = math.log(1 + minn_g_dist)
    minn_h_dist = math.log(1 + minn_h_dist)

    maxx_g_dist = math.log(1+maxx_g_dist)
    maxx_h_dist = math.log(1+maxx_h_dist)

    print(maxx_g_dist, maxx_h_dist)

    xbins = np.linspace(minn_g_dist, maxx_g_dist, num = nbins)
    ybins = np.linspace(minn_h_dist, maxx_h_dist, num = nbins)

    H_full = dict()

    for e_idx in range(len(labels)):

        xx_vals = []
        yy_vals = []

        skipping = 0

        for matching in all_matching[e_idx]:
            if matching[1] != math.inf:
                # Scale to reso (5kb)
                xx_vals.append(matching[0])
                yy_vals.append(matching[1]/reso)
                # Scale to 5bk reso and take log2
                #x_valls.append(math.log2(matching[0]/reso))
            #else:
                #x_valls.append(matching[0])
            else:
                skipping += 1

        if len(xx_vals) == 0:
            print('EMPTY, SKIPPING PLOT')
            continue


        print(f'Not found {skipping} peaks')
        # log scale
        xx_vals = np.log1p(np.array(xx_vals))
        yy_vals = np.log1p(np.array(yy_vals))

        H_full[e_idx], x_edges, y_edges =\
                np.histogram2d(xx_vals, yy_vals\
                                , bins=[xbins, ybins])

        #plt.imshow(H.T, cmap='gnuplot_r'\
        #            , norm = mpl.colors.LogNorm()\
        #            #, extent = (minn_g_dist, maxx_g_dist, minn_h_dist, maxx_h_dist)\
        #            , origin = 'lower'\
        #            )

        #plt.hist2d(xx_vals, yy_vals, bins=[xbins, ybins]\
        #            , cmap='gnuplot_r', cmin=1)

        #plt.colorbar()

        #plt.show()

        #exit()

        ##plt.scatter(x_valls, y_valls, alpha=0.7, color='blue')
        #plt.hist2d(x_valls, y_valls, bins=50, norm=mpl.colors.LogNorm(), cmap='gnuplot_r')
        #plt.colorbar()

        #plt.xlabel('genomic distance of cis-anchor (log 1+5kb)')
        #plt.ylabel('nearest hiccup centroid (log 1+5kb)')

        #plt.show()
        #plt.cla()
        #plt.clf()
        #plt.close()

    #print(maxx_g_dist, maxx_h_dist)


    graph_dict = dict()

    common_anchors = dict()



    g_dist = []

    # Initialize set intersection counts to 0
    s = pd.Series(np.zeros(len(tuples), dtype=int)\
        , index=index)

    for e_idx in labels_for_intersection:

        exp = labels[e_idx]

        e1_info = exp_count_info[exp][b_idx]

        anchor_list = e1_info['cis']

        for anchor in anchor_list:

            chrom = anchor[0]

            if chrom not in graph_dict:

                graph_dict[chrom] = nx.Graph()

            bin1 = anchor[1]
            bin2 = anchor[2]

            # Ignore anchor with genomic dist not more than anchor thresh
            # Rationale: birth threshold is 95%-tile of bin-dist 2
            if not abs(bin2 - bin1) > anchor_thresh:
                continue
            
            graph_dict[chrom].add_node((e_idx, bin1, bin2))

    for chrom in graph_dict:
        print('DOING CHROM', chrom)

        all_common_anchors[b_idx][chrom] = dict()

        for n1, n2 in itertools.combinations(list(graph_dict[chrom].nodes), 2):

            b11 = n1[1]
            b12 = n1[2]

            b21 = n2[1]
            b22 = n2[2]

            ddist = max(abs(b11-b21), abs(b12-b22))

            # Closeness criteria is within 2 bin distance b/c...
            # ...birth threshold is 95%-tile of bin distance 2
            if not ddist > anchor_thresh:
                graph_dict[chrom].add_edge(n1, n2)

        # Find maximal cliques
        #print('finding cliques...')
        cliques = nx.find_cliques(graph_dict[chrom])
        #print('found')
        for cl in cliques:
            this_tuple = np.zeros((len(labels_for_intersection,)), dtype=int)
            for node in cl:
                e_idx = node[0]
                this_tuple[labels_for_intersection.index(e_idx)] = 1

            # If tuple is all ones
            if np.amax(np.abs(this_tuple - 1)) == 0:
                for node in cl:

                    this_e_idx = node[0]

                    # Add this node to common anchors dictionary
                    this_exp = labels[this_e_idx]
                    if this_exp not in all_common_anchors[b_idx][chrom]:
                        all_common_anchors[b_idx][chrom][this_exp] = []

                    all_common_anchors[b_idx][chrom][this_exp].append([node[1], node[2]])

                    this_g_dist = abs(node[1] - node[2])

                    a1_bin1 = node[1]*reso
                    a1_bin2 = node[2]*reso

                    minnn = math.inf
                    entry = []
                    for hiccup in hiccups_data[this_e_idx]['chr'+chrom]:

                        c1 = min(hiccup[0], hiccup[1])
                        c2 = max(hiccup[0], hiccup[1])

                        val = max(abs(a1_bin1-c1), abs(a1_bin2-c2))

                        if val < minnn:
                            minnn = val
                            entry = hiccup

                    all_matching_common[this_e_idx].append([this_g_dist, minnn] + entry)

            s[tuple(this_tuple)] += 1

    H_common = dict()

    for h_e_idx in range(len(hiccup_files)):

        print('MATCHING WITH', hiccup_files[h_e_idx])

        x_valls = []
        y_valls = []

        skipping = 0

        for matching in all_matching_common[h_e_idx]:
            if matching[1] != math.inf:
                # Scale to reso (5kb)
                x_valls.append(matching[0])
                y_valls.append(matching[1]/reso)
                # Scale to 5bk reso and take log2
                #x_valls.append(math.log2(matching[0]/reso))
            #else:
                #x_valls.append(matching[0])
            else:
                skipping += 1

        if len(x_valls) == 0:
            print('EMPTY, SKIPPING PLOT')
            continue

        print(f'Not found {skipping} peaks')
        # log scale
        x_valls = np.log1p(np.array(x_valls))
        y_valls = np.log1p(np.array(y_valls))

        #plt.scatter(x_valls, y_valls, alpha=0.7, color='blue')

        H_common[h_e_idx], x_edges, y_edges =\
                np.histogram2d(x_valls, y_valls\
                                , bins=[xbins, ybins])

        #plt.hist2d(x_valls, y_valls, bins=50, norm=mpl.colors.LogNorm(), cmap='gnuplot_r')

        #plt.colorbar()

        #plt.xlabel('genomic distance of cis-anchor (log 1+5kb)')
        #plt.ylabel('nearest hiccup centroid (log 1+5kb)')

        #plt.show()
        #plt.cla()
        #plt.clf()
        #plt.close()

    #H_diff = dict()

    for e_idx in range(len(labels)):

        #fig, axs = plt.subplots(1, 3, figsize=(12,6), sharey=True, sharex=True)

        #H_diff[e_idx] = (H_full[e_idx] - H_common[e_idx])/H_full[e_idx]*100

        plt.imshow(H_full[e_idx].T, cmap='gnuplot_r'\
                    #, norm = mpl.colors.LogNorm()\
                    #, extent = (minn_g_dist, maxx_g_dist, minn_h_dist, maxx_h_dist)\
                    , origin = 'lower'\
                    )

        plt.colorbar()

        plt.xlabel('genomic distance log(1+5kb)')
        plt.ylabel('nearest huccup peak log(1+5kb)')

        plt.savefig('nearest_hiccup_ALL_exp'+labels[e_idx]+'_bidx'+str(b_idx)+'.pdf', dpi=600)
        plt.cla()
        plt.clf()
        plt.close()


        plt.imshow(H_common[e_idx].T, cmap='gnuplot_r'\
                    #, norm = mpl.colors.LogNorm()\
                    #, extent = (minn_g_dist, maxx_g_dist, minn_h_dist, maxx_h_dist)\
                    , origin = 'lower'\
                    )

        plt.colorbar()
        plt.xlabel('genomic distance log(1+5kb)')
        plt.ylabel('nearest huccup peak log(1+5kb)')

        plt.savefig('nearest_hiccup_COMMON_exp'+labels[e_idx]+'_bidx'+str(b_idx)+'.pdf', dpi=600)

        #axs[2].imshow(H_diff[e_idx].T, cmap='gnuplot_r'\
        #            #, norm = mpl.colors.LogNorm()\
        #            #, extent = (minn_g_dist, maxx_g_dist, minn_h_dist, maxx_h_dist)\
        #            , origin = 'lower'\
        #            )

        ##plt.colorbar()

        #plt.show()

        plt.cla()
        plt.clf()
        plt.close()

    #exit()


    #print(s[(1,1,1,0,0)])
    #print(s[(1,1,1)])
    for tt in tuples:
        print(tt, s[tt])
    
    if log_scale:
        for tt in tuples:
            s[tt] = math.log(1 + s[tt])
    
    upsetplot(s)

    axs = plt.gca()

    axs.set_ylabel('log(1+Intersection size)')
    
    #plt.show()

    plt.savefig('cis_upset_plot_5kb_bidx'+str(b_idx)+'.pdf', dpi=600)

    #exit()
    
    plt.cla()
    plt.clf()
    plt.close()


pickle.dump(all_common_anchors, open('cis_common_anchors.p', 'wb'))



## GET TRANS COMMON ANCHORS
#all_common_anchors = dict()





