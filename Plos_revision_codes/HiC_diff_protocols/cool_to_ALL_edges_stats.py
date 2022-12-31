import matplotlib.pyplot as plt
import numpy as np
import pandas
import math
import h5py
from numba import njit, prange
from numba import njit
import pickle

import cooler 

@njit
def search_edges(n_elements, bin1, bin2, pixels, weights, g_dist, edges\
                , count):
    
    for i in range(n_elements):

        bin1_e = bin1[i]
        bin2_e = bin2[i]

        if bin1_e == bin2_e:
            continue

        val = pixels[i]
        raw_bal_val = 1/val

        bal_val = val*weights[bin1_e]*weights[bin2_e]
        bal_val = 1/bal_val


        binn1 = min(bin1_e, bin2_e)
        binn2 = max(bin1_e, bin2_e)

        if binn2 - binn1 > g_dist:
            continue

        edges[count, 0] = binn1
        edges[count, 1] = binn2
        edges[count, 2] = raw_bal_val
        edges[count, 3] = bal_val
        count += 1

    return edges, count

@njit
def save_all_edges(jj, bin1, bin2, pixels, all_edges, count):

    for i in range(jj):
        val = pixels[i]
        bin1_e = bin1[i]
        bin2_e = bin2[i]

        # +1 for MaxHiC interactions file indexing
        all_edges[count, 0] = bin1_e+1
        all_edges[count, 1] = bin2_e+1
        all_edges[count, 2] = val

        count += 1

    return all_edges, count




# This creates only plots for genomic dist vs. distance estimate for ICE correction

exps_dir = '../cooler_files/'

exps = [\
         'fa_dsg_ddel_5kb.cool'\
        ,'fa_dsg_dpn_5kb.cool'\
        ,'fa_dsg_ddel_dpn_5kb.cool'\
        ,'fa_dsg_mnase_5kb.cool'\
        ,'fa_dpn_5kb.cool'\
        ,'fa_dsg_ddel_1kb.cool'\
        ,'fa_dsg_dpn_1kb.cool'\
        ,'fa_dsg_ddel_dpn_1kb.cool'\
        ,'fa_dsg_mnase_1kb.cool'\
        ,'fa_dpn_1kb.cool'\
       ]

labels = [\
         'fa_dsg_ddel_5kb'\
        ,'fa_dsg_dpn_5kb'\
        ,'fa_dsg_ddel_dpn_5kb'\
        ,'fa_dsg_mnase_5kb'\
        ,'fa_dpn_5kb'\
        ,'fa_dsg_ddel_1kb'\
        ,'fa_dsg_dpn_1kb'\
        ,'fa_dsg_ddel_dpn_1kb'\
        ,'fa_dsg_mnase_1kb'\
        ,'fa_dpn_1kb'\
        ]

reso_num = [\
          '5kb'\
        , '5kb'\
        , '5kb'\
        , '5kb'\
        , '5kb'\
        , '1kb'\
        , '1kb'\
        , '1kb'\
        , '1kb'\
        , '1kb'\
        ]

save_dir = '../cooler_results'


fig, axs = plt.subplots(1, 2)
all_ptiles = dict()

for idx, exp in enumerate(exps):

    exp = exps_dir + exp

    all_ptiles[labels[idx]] = dict()

    h5 = h5py.File(exp, 'r')

    n_elements = h5['pixels']['bin1_id'].shape[0]
    weights = h5['bins']['weight'][:]

    reso = reso_num[idx]

    #if reso == '1kb':
    #    break

    all_ptiles[labels[idx]][reso] = dict()

    cfile = exp
    
    cf = cooler.Cooler(cfile)

    print(f'Number of non-zero pixels for exp {exp} reso {reso} is {n_elements}')

    chunk_size = 50000

    edges = np.ones((int(n_elements/5), 4))
    count = 0
    
    i_start = 0
    i_end = chunk_size
    
    # threshold for bin distance
    g_thresh = 10

    while (i_start < n_elements):
        if (i_end > n_elements):
            i_end = n_elements
    
        if (i_start%10000==0):
            print('i_start/n_elements', i_start/n_elements, end='\r')
    
        bin1 = h5['pixels']['bin1_id'][i_start:i_end]
        bin2 = h5['pixels']['bin2_id'][i_start:i_end]
        pixels = h5['pixels']['count'][i_start:i_end]
    
        edges, count = search_edges(i_end-i_start, bin1, bin2, pixels, weights, g_thresh, edges\
                                    , count)

    
        i_start += chunk_size
        i_end += chunk_size

    edges = edges[:count]

    filename = '../cooler_results/edge_files/'+labels[idx]+'_raw_edges_upto_bindist'+str(g_thresh)+'.csv'
    ff = open(filename, 'w')

    filename = '../cooler_results/edge_files/'+labels[idx]+'_bal_edges_upto_bindist'+str(g_thresh)+'.csv'
    gg = open(filename, 'w')

    for ee in edges:
        ff.write(str(ee[0])+','+str(ee[1])+','+str(ee[2])+'\n')
        gg.write(str(ee[0])+','+str(ee[1])+','+str(ee[3])+'\n')

    ff.close()
    gg.close()




#    bin_dist = edges[:, 1] - edges[:, 0]
#
#
#    xx = []
#    raw_yy = []
#    bal_yy = []
#
#    for g_dist in list(range(1, g_thresh+1)):
#
#        this_idxs = np.argwhere(bin_dist == g_dist).flatten()
#
#        this_raw_vals = edges[this_idxs, 2]
#        this_bal_vals = edges[this_idxs, 3]
#
#        print(this_raw_vals, this_bal_vals)
#
#        # remove NaN from balanced
#        nan_idxs = np.isnan(this_bal_vals)
#        not_nan_idxs = np.argwhere(nan_idxs==False).flatten()
#        this_bal_vals = this_bal_vals[not_nan_idxs]
#
#        raw_95 = np.percentile(this_raw_vals, 95)
#        bal_95 = np.percentile(this_bal_vals, 95)
#
#        xx.append(g_dist)
#        raw_yy.append(raw_95)
#        bal_yy.append(bal_95)
#        
#
#    axs[0].plot(xx, raw_yy, label=labels[idx])
#    axs[1].plot(xx, bal_yy)
#
#
#axs[0].legend()
#plt.show()
#
#
#    
#    #if thresh == math.inf:
#    #    thresh_label = 'inf'
#    #else:
#    #    thresh_label = str(thresh)
#
#    ##print('EXITING AFTER MaxHiC interactions file for testing')
#
#    #g_dists = list(range(1, 11))
#
#    #this_color = 'black'
#
#
#
#    ## 1. Plot genomic distance vs. distance estimate
#    #for ch in range(n_chrom):
#
#
#    #    print("\nSaving edges for chrom", chrom_names[ch])
#    #    save_edges = edges[ch, :count[ch]]
#
#    #    if chrom_names[ch] == 'chrX':
#    #        this_color = 'blue'
#    #    elif chrom_names[ch] == 'chrY':
#    #        this_color = 'red'
#    #    else:
#    #        this_color = 'black'
#
#    #    genomic_dist = np.abs(save_edges[:,0] - save_edges[:,1])
#    #    dist = save_edges[:,2]
#
#    #    minns = []
#    #    maxxs = []
#    #    medians = []
#
#
#
#    #    all_ptiles[labels[idx]][reso][chrom_names[ch]] = dict()
#
#    #    percentiles = [5, 10, 25, 50, 75, 90, 95]
#
#
#    #    for g_dist in g_dists:
#
#    #        idxs = np.argwhere(genomic_dist == g_dist)
#    #        this_dist = dist[idxs]
#
#    #        this_ptiles = np.percentile(this_dist, percentiles)
#
#    #        all_ptiles[labels[idx]][reso][chrom_names[ch]][g_dist] = dict()
#
#    #        for p_idx, ptile in enumerate(this_ptiles):
#
#    #            all_ptiles[labels[idx]][reso][chrom_names[ch]][g_dist][percentiles[p_idx]] = ptile
#
#    #        #try:
#    #        #    all_p_tiles.append(np.percentile(this_dist, percentiles))
#    #        #except:
#    #        #    # quick fix for 3 percentiles
#    #        #    all_p_tiles.append(np.array([0, 0, 0]))
#
#    #        #minns.append(np.amin(this_dist))
#    #        #maxxs.append(np.amax(this_dist))
#    #        #medians.append(np.median(this_dist))
#
#
#    #    #plt.plot(g_dists, minns, color='blue', marker='o')
#    #    #plt.plot(g_dists, maxxs, color='red', marker='o')
#    #    #plt.plot(g_dists, medians, color='black', marker='o')
#
#    #    #plt.show()
#
#
#    #    #all_p_tiles = np.array(all_p_tiles)
#
#    #    #all_p_tiles_list.append(all_p_tiles)
#
#    #
#    ##print(all_ptiles[labels[idx]][reso])
#
#    ##exit()
#
#    ##fig, axs = plt.subplots(1, 3)
#
#    ##for ch in range(n_chrom):
#
#    ##    p_tiles_array = all_p_tiles_list[ch]
#
#    ##    for p_idx, ptile in enumerate(percentiles):
#    ##        axs[p_idx].plot(g_dists, p_tiles_array[:, p_idx], color=this_color\
#    ##                , alpha = (p_idx+1)/(len(percentiles)+1)\
#    ##                , label=str(ptile)+' percentile')
#
#    ###print(all_p_tiles)
#
#    ##axs[0].set_ylabel('distance estimate')
#
#    ##axs[0].set_xlabel('genomic distance')
#    ##axs[1].set_xlabel('genomic distance')
#    ##axs[2].set_xlabel('genomic distance')
#
#    ##axs[0].set_title('percentile '+str(percentiles[0]))
#    ##axs[1].set_title('percentile '+str(percentiles[1]))
#    ##axs[2].set_title('percentile '+str(percentiles[2]))
#
#    ###plt.xlabel('genomic distance')
#    ###plt.ylabel('distance estimate')
#
#    ##plt.tight_layout()
#
#    ###plt.legend()
#    ##plt_file = save_dir+'/'+labels[idx]+'_genomic_vs_dist_allchroms_reso'+reso+'.pdf'
#
#    ##plt.savefig(plt_file, dpi=600)
#    ###plt.show()
#    ##
#    ##plt.cla()
#    ##plt.clf()
#    #
#    ##exit()
#    
#print('saving stats')
#pickle.dump(all_ptiles, open('ice_norm_stats.p', 'wb'))
#
#
#
#
#
#
