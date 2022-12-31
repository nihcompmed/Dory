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

dim = 1

for b_idx in range(2):

    info = pickle.load(open('common_cis_anchors_dim'+str(dim)+'_bidx'+str(b_idx)+'.p', 'rb'))

    #all_nearest = []

    for tt in info:

        #print(info[tt])
        strtt = ''
        for t in tt:
            strtt = strtt+str(t)
        print(strtt)


        bed_file = 'cis_combination_'+strtt+'_atac_dim'+str(dim)+'_bidx'+str(b_idx)+'.bed'
        #bed_file = 'test.bed'

        bed = open(bed_file, 'w')

        for chrom in info[tt]:

            print(chrom)

            unique_bed_lines = []

            #ccc = 0

            for anchor in info[tt][chrom]:

                #print(anchor)

                a1 = anchor[0]

                start = a1*reso
                end = (a1+1)*reso

                bed_lines = hf2.get_narrow_range(atac_narrow_file, 'chr'+chrom, start, end)
                #ccc += len(bed_lines)
                for ll in bed_lines:
                    if ll not in unique_bed_lines:
                        unique_bed_lines.append(ll)

                a2 = anchor[1]

                start = a2*reso
                end = (a2+1)*reso
                #start = max((a2-(b_idx+1))*reso, 0)
                #end = (a2+(b_idx+1))*reso

                bed_lines = hf2.get_narrow_range(atac_narrow_file, 'chr'+chrom, start, end)
                #ccc += len(bed_lines)
                for ll in bed_lines:
                    if ll not in unique_bed_lines:
                        unique_bed_lines.append(ll)

                #ccc2 = len(unique_bed_lines)

                #if ccc != ccc2:
                #    print('not equal', ccc, ccc2)


            for ll in unique_bed_lines:
                for entry in ll[:-1]:
                    bed.write(entry+'\t')
                bed.write(ll[-1])

        bed.close()

    #exit()
    

#    # First get nearest ATAC peak
#    for chrom in info:
#        print('DOING CHROM', chrom)
#        count += 1
#        anchors = info[chrom]
#
#        for anchor in anchors:
#
#            nearest = hf2.find_nearest_peak_bed('chr'+chrom, anchor, reso, atac_narrow_file)
#
#            all_nearest.append(nearest)
#
#            start = (anchor-1)*reso
#            end = (anchor+1)*reso
#            bed_lines_5kb = hf2.get_narrow_range(atac_narrow_file, 'chr'+chrom, start, end)
#
#            for ll in bed_lines_5kb:
#                for entry in ll[:-1]:
#                    bed_5kb.write(entry+'\t')
#                bed_5kb.write(ll[-1])
#            
#
#            start = (anchor-2)*reso
#            end = (anchor+2)*reso
#            bed_lines_10kb = hf2.get_narrow_range(atac_narrow_file, 'chr'+chrom, start, end)
#            for ll in bed_lines_10kb:
#                for entry in ll[:-1]:
#                    bed_10kb.write(entry+'\t')
#                bed_10kb.write(ll[-1])
#
#    bed_5kb.close()
#    bed_10kb.close()

        #if count == 2:
        #    break
#
#    # Scale to 1kb and log (1+x)
#    all_nearest = np.log1p(np.array(all_nearest)/1000)


#    sns.displot(x=all_nearest, kde=True)
#    plt.yscale('log', base=2)
#
#    # Mark ref_bin (in kb units) on xaxis
#    ref_bin = 5
#    xx_ref = math.log(1 + ref_bin)
#    plt.axvline(x=xx_ref, color='red')
#    plt.text(xx_ref, plt.gca().get_ylim()[1], str(ref_bin)+'kb')
#
#
#    # Mark ref_bin (in kb units) on xaxis
#    ref_bin = 10
#    xx_ref = math.log(1 + ref_bin)
#    plt.axvline(x=xx_ref, color='red')
#    plt.text(xx_ref, plt.gca().get_ylim()[1], str(ref_bin)+'kb')
#
#    plt.xlabel('nearest ATAC peak (log 1 + kb)')
#
#    plt.savefig('histogram_atac_nearness_cis_dim'+str(dim)+'_bidx'+str(b_idx)+'.pdf'\
#            , dpi=600)
#








