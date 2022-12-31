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

hiccup_files = [\
         '../fa_dsg_ddel_hiccups.bedpe'\
        ,'../fa_dsg_dpn_hiccups.bedpe'\
        ,'../fa_dsg_mnase_hiccups.bedpe'\
        ,'../fa_dsg_ddel_dpn_hiccups.bedpe'\
        ,'../fa_dpn_hiccups.bedpe'\
                ]

hiccups_data = []
for ff in hiccup_files:
    hiccups_data.append(hf2.get_hiccup_peaks(ff))


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

for h_idx, hiccup in enumerate(hiccup_files):

    print(hiccup)
    #exit()

    test_bed_file = open('test.bed', 'w')

    h_data = hiccups_data[h_idx]



    for chrom in h_data:

        for peak in h_data[chrom]:

            start1 = peak[2]
            end1 = peak[3]

            # Get ATAC peaks in the anchor bins
            atac_peaks_1 = hf2.get_narrow_range(atac_narrow_file, chrom\
                                            , start1, end1)
            start2 = peak[4]
            end2 = peak[5]

            for ll in atac_peaks_1:
                for entry in ll[:-1]:
                    test_bed_file.write(entry+'\t')
                test_bed_file.write(ll[-1])

            atac_peaks_2 = hf2.get_narrow_range(atac_narrow_file, chrom\
                                            , start2, end2)

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
