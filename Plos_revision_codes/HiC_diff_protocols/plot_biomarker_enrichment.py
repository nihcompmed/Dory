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


def get_chrom(binn, chrom_start_bins_list, chrom_names):

    for idx, start in enumerate(chrom_start_bins_list):

        if binn < start:
            return chrom_names[idx-1], binn-chrom_start_bins_list[idx-1]

#PH_pars = pickle.load(open('../cooler_results/edge_files/PH_pars_info_5kb.p', 'rb'))

# Resolution 5kb first

labels = [\
         'fa_dsg_ddel_5kb'\
        ,'fa_dsg_dpn_5kb'\
        ,'fa_dsg_ddel_dpn_5kb'\
        ,'fa_dsg_mnase_5kb'\
        ,'fa_dpn_5kb'\
        ]

exps = [\
         'fa_dsg_ddel_5kb.cool'\
        ,'fa_dsg_dpn_5kb.cool'\
        ,'fa_dsg_ddel_dpn_5kb.cool'\
        ,'fa_dsg_mnase_5kb.cool'\
        ,'fa_dpn_5kb.cool'\
       ]



data_dirr = '../Biomarkers/'


# BED PEAKS
fnames = [\
        data_dirr+'ENCFF294RSZ.bed'\
        ,data_dirr+'GSE165895_HFFc6_Atacseq.mRp.clN_peaks.narrowPeak'\
        ,data_dirr+'GSE165895_HFFc6_CTCF_CT.mRp.clN_peaks.narrowPeak'\
        ]

flabel = [\
        'ChIPseq_CTCF'\
        ,'Atacseq'\
        ,'CTCF'
        ]

atac_narrow_file = fnames[1]

# BIGWIG FILES
bigwig_fnames = [\
         data_dirr+'GSE165895_HFFc6_H3K27ac.bigWig'\
        ,data_dirr+'GSE165895_HFFc6_H3K4me3.bigWig'\
        ,data_dirr+'GSE165895_HFFc6_SMC1_CT.mRp.clN.bigWig'\
        ]

bigwig_flabel = [\
        'H3K27ac'\
        ,'H3K4me3'\
        ,'SMC1'\
        ]

#PH_pars = pickle.load(open('../cooler_results/edge_files/PH_pars_info_5kb.p', 'rb'))


reso = 5000

n_chroms = 24

counts = np.array([0, 0, 0])

for e_idx, exp in enumerate(labels):

    print(exp)
    if e_idx != 3:
        continue
    #exit()

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


    test_bed_file = open('test.bed', 'w')


    for idx in range(1, 2):

        print(exp, idx)

        cycle_info = []

        cis_cycles = []

        target = '../cooler_results/PH_results/'+exp+'_gdist'+str(idx+1)+'_ALL_'

        # H1 shortened cycles (and have been checked for connecedness)
        cyc_file = target + 'minimal_V_birth_H1.txt'

        cycs = open(cyc_file, 'r')

        progress = 0

        for line in cycs:

            print(progress, end='\r')
            progress += 1

            line = line.split(',')
            line = line[:-1]
            line = [int(x) for x in line]

            # Since this is after connectedness check
            # The cycle is stored as v_i, v_{i+1}
            # Append the first element to end to complete the cycle
            line.append(line[0])
            line = np.array(line)

            this_g_dist = np.abs(line[1:] - line[:-1])

            # multiple anchors?
            # Simply assume that any edge with more than one 1 genomic distance is an anchor
            # Later plot the statistics
            all_possible_anchors = np.argwhere(this_g_dist > 1).flatten()

            # OR pick the maximum one
            all_possible_anchors = np.argmax(this_g_dist).flatten()

            #print(all_possible_anchors)

            #print(len(all_possible_anchors))

            for anchor in all_possible_anchors:

                anchor_b1 = line[anchor]
                anchor_b2 = line[anchor+1]

                chrr_1, binn_1 = get_chrom(anchor_b1, chrom_start_bins_list, chrom_names)
                chrr_2, binn_2 = get_chrom(anchor_b2, chrom_start_bins_list, chrom_names)

                # Get ATAC peaks in the anchor bins
                atac_peaks_1 = hf2.get_narrow_range(atac_narrow_file, 'chr'+chrr_1\
                                                , (anchor_b1-1)*reso, (anchor_b1+2)*reso)

                for ll in atac_peaks_1:
                    for entry in ll[:-1]:
                        test_bed_file.write(entry+'\t')
                    test_bed_file.write(ll[-1])

                atac_peaks_2 = hf2.get_narrow_range(atac_narrow_file, 'chr'+chrr_2\
                                                , (anchor_b2-1)*reso, (anchor_b2+2)*reso)

                for ll in atac_peaks_2:
                    for entry in ll[:-1]:
                        test_bed_file.write(entry+'\t')
                    test_bed_file.write(ll[-1])
                #print(binn_1, binn_2)
                #print(atac_peaks_1, atac_peaks_2)

                if len(atac_peaks_1) and len(atac_peaks_2):
                    counts[2] += 1
                elif not len(atac_peaks_1) and not len(atac_peaks_2):
                    counts[0] += 1
                else:
                    counts[1] += 1

                print(counts)

        test_bed_file.close()
        exit()
        continue

#            iddx = np.argmax(this_g_dist)
#            maxx = np.max(this_g_dist)
#
#            #if maxx < 500:
#            #    continue
#
#            #if iddx > 1:
#            #    print('Multiple maximal anchors with g dist ', maxx)
#
#            anchor_b1 = line[iddx]
#            anchor_b2 = line[iddx+1]
#
#            chrr_1, binn_1 = get_chrom(anchor_b1, chrom_start_bins_list, chrom_names)
#            chrr_2, binn_2 = get_chrom(anchor_b2, chrom_start_bins_list, chrom_names)
#
#            print(chrr_1, binn_1, chrr_2, binn_2)
#
#            ###########################################
#            # Narrow peaks
#            ###########################################
#            fig, axs = plt.subplots(3, 2, figsize=(10, 8))
#
#            if chrr_1 != chrr_2:
#                fig.suptitle('trans anchor, genomic dist '+str(maxx))
#            else:
#                fig.suptitle('cis anchor, genomic dist '+str(maxx))
#
#            row = 0
#            col = 0
#            axs[row, col].set_title('anchor 1, chrom '+chrr_1+', bin '+str(binn_1))
#
#            for idx, bfile in enumerate(fnames):
#
#                hf2.plot_narrow_range(axs[row,col], bfile\
#                                     , 'chr'+chrr_1\
#                                     , (binn_1-2)*reso\
#                                     , (binn_1+2)*reso\
#                                     )
#
#                y_lims = axs[row, col].get_ylim()
#                axs[row, col].scatter(binn_1*reso/1000, y_lims[1], color='red')
#
#                axs[row, col].set_ylabel(flabel[idx])
#
#                row += 1
#
#            row = 0
#            col = 1
#            axs[row, col].set_title('anchor 2, chrom '+chrr_2+', bin '+str(binn_2))
#
#            for idx, bfile in enumerate(fnames):
#
#                hf2.plot_narrow_range(axs[row,col], bfile\
#                                     , 'chr'+chrr_2\
#                                     , (binn_2-2)*reso\
#                                     , (binn_2+2)*reso\
#                                     )
#
#                y_lims = axs[row, col].get_ylim()
#                axs[row, col].scatter(binn_2*reso/1000, y_lims[1], color='red')
#
#
#                row += 1
#
#            
#            plt.tight_layout()
#            plt.show()
#            plt.cla()
#            plt.clf()
#            plt.close()
#
#            #exit()
#
#
#            ###########################################
#            # BigWig peaks
#            ###########################################
#            fig, axs = plt.subplots(3, 2, figsize=(10, 8))
#
#            if chrr_1 != chrr_2:
#                fig.suptitle('trans anchor, genomic dist '+str(maxx))
#            else:
#                fig.suptitle('cis anchor, genomic dist '+str(maxx))
#
#            row = 0
#            col = 0
#            axs[row, col].set_title('anchor 1, chrom '+chrr_1+', bin '+str(binn_1))
#
#            for idx, bfile in enumerate(bigwig_fnames):
#
#                hf2.plot_bigwig_range(axs[row,col], bfile\
#                                     , 'chr'+chrr_1\
#                                     , (binn_1-2)*reso\
#                                     , (binn_1+2)*reso\
#                                     )
#
#                y_lims = axs[row, col].get_ylim()
#                axs[row, col].scatter(binn_1*reso/1000, y_lims[1], color='red')
#
#                axs[row, col].set_ylabel(bigwig_flabel[idx])
#
#                row += 1
#
#
#            row = 0
#            col = 1
#            axs[row, col].set_title('anchor 2, chrom '+chrr_2+', bin '+str(binn_2))
#
#            for idx, bfile in enumerate(bigwig_fnames):
#
#                hf2.plot_bigwig_range(axs[row,col], bfile\
#                                     , 'chr'+chrr_2\
#                                     , (binn_2-2)*reso\
#                                     , (binn_2+2)*reso\
#                                     )
#
#                y_lims = axs[row, col].get_ylim()
#                axs[row, col].scatter(binn_2*reso/1000, y_lims[1], color='red')
#
#                axs[row, col].set_ylabel(bigwig_flabel[idx])
#
#                row += 1
#
#
#            plt.tight_layout()
#
#            plt.show()
#
#
#
#        exit()
#    
#
#
#
#
#
#
#
#
#
