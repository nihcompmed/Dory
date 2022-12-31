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

def update_cis_trans(c1, bin1, c2, bin2\
                    , cis_anchors, trans_anchors\
                    , anchor_thresh):

    # Edge c1-c2
    if c1 == c2:
        this_g_dist = abs(bin2 - bin1)
        if this_g_dist > anchor_thresh:
            cis_anchors.append((c1, min(bin1, bin2), max(bin1, bin2)))
    else:
        trans_anchors.append((c1, bin1, c2, bin2))

    return cis_anchors, trans_anchors

reso = 5000

#reso = 10000

labels = [\
         'fa_dsg_ddel_5kb'\
        ,'fa_dsg_dpn_5kb'\
        ,'fa_dsg_ddel_dpn_5kb'\
        ,'fa_dsg_mnase_5kb'\
        ,'fa_dpn_5kb'\
        ]

scale= [
        1\
        ,1\
        ,0.92\
        ,0.75\
        ,1.15\
        ]

shift = [
        0\
        ,0\
        ,0\
        ,0\
        ,0\
        ]

mu = [\
        14.029\
        ,15.822\
        ,10.418\
        ,5.948\
        ,45.507\
    ]

sigma = [\
        3.9248\
        ,4.4154\
        ,2.992\
        ,1.9266\
        ,12.3218\
        ]

posshift = [\
        3.342\
        ,3.331\
        ,3.223\
        ,2.644\
        ,3.3912\
            ]

# First do for H1
dim = 1
## Do for H2
#dim = 2
print('Processing dim', dim)

exp_counts = dict()

#anchor_thresh = 2


for e_idx, exp in enumerate(labels):

    # cfile
    cfile = '../cooler_files/'+exp + '.cool'

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

    if exp not in exp_counts:
        exp_counts[exp] = dict()

    for b_idx in range(2):

        print(exp, b_idx)

        # ANCHOR THRESH IS 1 if birth idx is for gene distance 1
        # ANCHOR THRESH IS 2 if birth idx is for gene distance 2
        anchor_thresh = b_idx + 1

        exp_counts[exp][b_idx] = dict()
        
        t1 = '../cooler_results/PH_results/'+exp+'_gdist'+str(b_idx+1)+'_ALL_scaled_'

        # Assuming they are checked for connectedness, hence
        # Cycles are stored as v0, v1, v2... where v_i, v_{i+1} is an edge
        cyc_file = t1+'minimal_V_birth_H'+str(dim)+'.txt'

        cis_anchors = []
        trans_anchors = []

        try:
            cycs = open(cyc_file, 'r')
        except:
            exp_counts[exp][b_idx]['cis'] = cis_anchors
            exp_counts[exp][b_idx]['trans'] = trans_anchors
            continue

        
        progress = 0

        if dim == 1:

            for line in cycs:

                print(progress, end='\r')
                progress += 1

                line = line.split(',')
                line = line[:-1]
                line = [int(x) for x in line]
                line.append(line[0])

                line = np.array(line)

                cyc_chroms = []

                for binn in line:
                    cyc_chroms.append(hf2.get_chrom(binn, chrom_start_bins_list, chrom_names))

                nn = 0
                while nn < len(line)-1:
                    c1 = cyc_chroms[nn][0]
                    bin1 = cyc_chroms[nn][1]

                    c2 = cyc_chroms[nn+1][0]
                    bin2 = cyc_chroms[nn+1][1]

                    cis_anchors, trans_anchors = \
                            update_cis_trans(c1, bin1, c2, bin2\
                                            , cis_anchors, trans_anchors\
                                            , anchor_thresh\
                                            )
                    nn += 1

        elif dim == 2:

            for line in cycs:

                print(progress, end='\r')
                progress += 1

                line = line.split(',')
                line = line[:-1]
                line = [int(x) for x in line]
                cyc_chroms = []

                for binn in line:
                    cyc_chroms.append(hf2.get_chrom(binn, chrom_start_bins_list, chrom_names))

                nn = 0
                while nn < len(line):
                    c1 = cyc_chroms[nn][0]
                    bin1 = cyc_chroms[nn][1]

                    c2 = cyc_chroms[nn+1][0]
                    bin2 = cyc_chroms[nn+1][1]

                    c3 = cyc_chroms[nn+2][0]
                    bin3 = cyc_chroms[nn+2][1]

                    # Edge c1-c2
                    cis_anchors, trans_anchors = \
                            update_cis_trans(c1, bin1, c2, bin2\
                                            , cis_anchors, trans_anchors\
                                            , anchor_thresh\
                                            )

                    # Edge c1-c3
                    cis_anchors, trans_anchors = \
                            update_cis_trans(c1, bin1, c3, bin3\
                                            , cis_anchors, trans_anchors\
                                            , anchor_thresh\
                                            )

                    # Edge c2-c3
                    cis_anchors, trans_anchors = \
                            update_cis_trans(c2, bin2, c3, bin3\
                                            , cis_anchors, trans_anchors\
                                            , anchor_thresh\
                                            )

                    nn += 3


                
        #print(cis_anchors, trans_anchors)
        exp_counts[exp][b_idx]['cis'] = cis_anchors
        exp_counts[exp][b_idx]['trans'] = trans_anchors


print('saving anchor valency for dim', dim)
pickle.dump(exp_counts, open('exp_cis_trans_anchors_dim'+str(dim)+'.p', 'wb'))

#
#            
#
#            # cis or trans?
#
#            cyc_g_dists = np.abs(line[1:] - line[:-1])
#
#            cyc_g_dists = cyc_g_dists[cyc_g_dists > 1]
#
#            uni, counts = np.unique(cyc_g_dists, return_counts=True)
#
#            #print(uni, counts)
#
#            #uni = np.flip(uni)
#            #counts = np.flip(counts)
#
#
#            cum_counts = np.cumsum(counts)
#            print(uni, cum_counts)
#            plt.plot(uni, cum_counts, marker='o', markersize=10)
#        
#        plt.ylabel('cumulative counts')
#        plt.xlabel('genomic distance')
#        plt.xscale('log', base=2)
#        plt.show()
#        exit()
#
#
#
#
#
