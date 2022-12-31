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
def search_edges(n_elements, bin1, bin2, pixels, weights, thresh, edges\
                , count, chrom_start_bins):
    
    for i in range(n_elements):
        val = pixels[i]
        bin1_e = bin1[i]
        bin2_e = bin2[i]

        if not (bin1_e < bin2_e) :
            continue
        c1_idx = 0
        while (chrom_start_bins[c1_idx] < bin1_e):
            c1_idx += 1

        c2_idx = 0
        while (chrom_start_bins[c2_idx] < bin2_e):
            c2_idx += 1

        if (c1_idx != c2_idx):
            continue

        c_idx = c1_idx - 1


        bal_val = val*weights[bin1_e]*weights[bin2_e]
        bal_val = 1/bal_val

        if (bal_val < thresh):
            edges[c_idx, count[c_idx], 0] = bin1_e
            edges[c_idx, count[c_idx], 1] = bin2_e
            edges[c_idx, count[c_idx], 2] = bal_val
            count[c_idx] += 1

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

exps_dir = [\
        '../GM12878/'\
       ,'../GM12878/'\
       ,'../GM12878/'\
       ]

exps = [\
        '../GM12878/GM12878_combined_1kb.cool'\
       , '../GM12878/GM12878_combined_5kb.cool'\
       ,'../GM12878/GM12878_combined_10kb.cool'\
       ]

reso_num = ['1000', '10000']


thresh = 400

for idx, exp in enumerate(exps):

    if idx > 0:
        break

    exp_dir = exps_dir[idx]

    h5 = h5py.File(exp, 'r')

    n_elements = h5['pixels']['bin1_id'].shape[0]
    weights = h5['bins']['weight'][:]

    reso = reso_num[idx]

    cfile = exp
    
    cf = cooler.Cooler(cfile)

    chrom_start_bins_list = []
    chrom_names = []

    n_chrom = 0
    
    for ch in cf.chromnames:
        
        chrom_names.append(ch)
        beg, end = cf.extent(ch)
        chrom_start_bins_list.append(beg)
        n_chrom += 1

    chrom_start_bins_list.append(cf.extent(cf.chromnames[-1])[1])
    chrom_start_bins = np.array(chrom_start_bins_list)

    print(f'Number of non-zero pixels for exp {exp} reso {reso} is {n_elements}')

    chunk_size = 50000

    edges = np.ones((n_chrom, 20000000, 3))
    count = np.zeros((n_chrom, ), dtype=int)
    
    i_start = 0
    i_end = chunk_size
    
    while (i_start < n_elements):
        if (i_end > n_elements):
            i_end = n_elements
    
        if (i_start%10000==0):
            print('i_start/n_elements', i_start/n_elements, end='\r')
    
        bin1 = h5['pixels']['bin1_id'][i_start:i_end]
        bin2 = h5['pixels']['bin2_id'][i_start:i_end]
        pixels = h5['pixels']['count'][i_start:i_end]
    
        edges, count = search_edges(i_end-i_start, bin1, bin2, pixels, weights, thresh, edges\
                                    , count, chrom_start_bins)

    
        i_start += chunk_size
        i_end += chunk_size
    
    if thresh == math.inf:
        thresh_label = 'inf'
    else:
        thresh_label = str(thresh)



    #print('EXITING AFTER MaxHiC interactions file for testing')

    all_p_tiles_list = []

    g_dists = list(range(1, 11))

    this_color = 'black'

    # 1. Plot genomic distance vs. distance estimate
    for ch in range(n_chrom):


        print("\nSaving edges for chrom", chrom_names[ch])
        save_edges = edges[ch, :count[ch]]


        filename = exp_dir + 'chr'+chrom_names[ch] + '_edges_reso'+str(reso) + '_thresh' + str(thresh) + '.csv'



        ff = open(filename, 'w')
        for row in save_edges:
            row[0] = row[0] - chrom_start_bins_list[ch]
            row[1] = row[1] - chrom_start_bins_list[ch]
            ff.write(str(int(row[0])) + ',' + str(int(row[1])) + ',' + str(round(row[2], 2))+'\n')

        ff.close()






