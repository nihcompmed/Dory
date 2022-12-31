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
import pandas as pd
from upsetplot import plot as upsetplot


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

ccre_types = [\
            'Low-DNase'\
            ,'DNase-only'\
            ,'DNase-H3K4me3'\
            ,'dELS'\
            ,'pELS'\
            ,'PLS'\
            ]

# Initialize Pandas multi index
tuples = [\
          (0,0,0,0,0,0)\
          ,(0,0,0,0,0,1)\
          ,(0,0,0,0,1,0)\
          ,(0,0,0,0,1,1)\
          ,(0,0,0,1,0,0)\
          ,(0,0,0,1,0,1)\
          ,(0,0,0,1,1,0)\
          ,(0,0,0,1,1,1)\
          ,(0,0,1,0,0,0)\
          ,(0,0,1,0,0,1)\
          ,(0,0,1,0,1,0)\
          ,(0,0,1,0,1,1)\
          ,(0,0,1,1,0,0)\
          ,(0,0,1,1,0,1)\
          ,(0,0,1,1,1,0)\
          ,(0,0,1,1,1,1)\
          ,(0,1,0,0,0,0)\
          ,(0,1,0,0,0,1)\
          ,(0,1,0,0,1,0)\
          ,(0,1,0,0,1,1)\
          ,(0,1,0,1,0,0)\
          ,(0,1,0,1,0,1)\
          ,(0,1,0,1,1,0)\
          ,(0,1,0,1,1,1)\
          ,(0,1,1,0,0,0)\
          ,(0,1,1,0,0,1)\
          ,(0,1,1,0,1,0)\
          ,(0,1,1,0,1,1)\
          ,(0,1,1,1,0,0)\
          ,(0,1,1,1,0,1)\
          ,(0,1,1,1,1,0)\
          ,(0,1,1,1,1,1)\
          ,(1,0,0,0,0,0)\
          ,(1,0,0,0,0,1)\
          ,(1,0,0,0,1,0)\
          ,(1,0,0,0,1,1)\
          ,(1,0,0,1,0,0)\
          ,(1,0,0,1,0,1)\
          ,(1,0,0,1,1,0)\
          ,(1,0,0,1,1,1)\
          ,(1,0,1,0,0,0)\
          ,(1,0,1,0,0,1)\
          ,(1,0,1,0,1,0)\
          ,(1,0,1,0,1,1)\
          ,(1,0,1,1,0,0)\
          ,(1,0,1,1,0,1)\
          ,(1,0,1,1,1,0)\
          ,(1,0,1,1,1,1)\
          ,(1,1,0,0,0,0)\
          ,(1,1,0,0,0,1)\
          ,(1,1,0,0,1,0)\
          ,(1,1,0,0,1,1)\
          ,(1,1,0,1,0,0)\
          ,(1,1,0,1,0,1)\
          ,(1,1,0,1,1,0)\
          ,(1,1,0,1,1,1)\
          ,(1,1,1,0,0,0)\
          ,(1,1,1,0,0,1)\
          ,(1,1,1,0,1,0)\
          ,(1,1,1,0,1,1)\
          ,(1,1,1,1,0,0)\
          ,(1,1,1,1,0,1)\
          ,(1,1,1,1,1,0)\
          ,(1,1,1,1,1,1)\
        ]

index = pd.MultiIndex.from_tuples(tuples, names=ccre_types)

reso = 5000

n_chroms = 24

counts = np.array([0, 0, 0])

ccre_dict = pickle.load(open('../Biomarkers/ccre_reso'+str(reso)+'.p', 'rb'))

# Load common anchors
dim = 1
print('Processing dim', dim)

for b_idx in range(2):

    # LOAD COMMON CLIQUES
    info = pickle.load(open('common_cis_cliques_dim'+str(dim)+'_bidx'+str(b_idx)+'.p', 'rb'))


    for exp_tt in info:

        s = pd.Series(np.zeros(len(tuples), dtype=int)\
                    , index=index)

        #print(info[tt])
        strtt = ''
        for t in exp_tt:
            strtt = strtt+str(t)
        print(exp_tt, strtt)

        count = 0

        pie_chart = np.array([0,0])

        for chrom in info[exp_tt]:

            print(chrom)

            this_ccre_dict = ccre_dict['chr'+chrom]

            for clique in info[exp_tt][chrom]:

                count += 1

                for anchor in clique:

                    binn1 = anchor[1]
                    start1 = binn1
                    end1 = binn1 + 1

                    binn2 = anchor[2]
                    start2 = binn2
                    end2 = binn2 + 1

                    #cre_set = []
                    this_tuple = np.zeros((6,), dtype=int)
                    for key in this_ccre_dict:

                        ccre_start = key[0]
                        ccre_end = key[1]

                        # If either ccre start or end is in the bin of the anchor
                        if\
                            (ccre_start >= start1 and ccre_start <= end1)\
                            or (ccre_end >= start1 and ccre_end <= end1)\
                            or (ccre_start >= start2 and ccre_start <= end2)\
                            or (ccre_end >= start2 and ccre_end <= end2)\
                            :

                            for cc in this_ccre_dict[key]:
                                c_idx = ccre_types.index(cc)
                                #print(c_idx)
                                this_tuple[c_idx] = 1
                                #cre_set.append(c_idx)


                #cre_set = np.array(list((frozenset(cre_set))), dtype=int)
                #this_tuple = np.zeros((6,), dtype=int)
                #this_tuple[cre_set] = 1
                this_tuple = tuple(this_tuple)
                s[this_tuple] += 1

                if \
                        this_tuple[1] == 0\
                    and this_tuple[2] == 0\
                    and this_tuple[3] == 0\
                    and this_tuple[4] == 0\
                    and this_tuple[5] == 0:
                        pie_chart[0] += 1
                else:
                        pie_chart[1] += 1



        #print(count)

        for t in tuples:
            s[t] = s[t]/count*100
            #s[t] = math.log1p(s[t])

        upsetplot(s)
        plt.title('Number of cliques '+str(count))

        plt.ylabel('Intersection percentage')

        #plt.show()

        #plt.ylabel('log(1+intersection size)')
        #plt.show()
        #exit()
        #plt.savefig('ccre_upset_dim'+str(dim)+'_bidx'+str(b_idx)+'.pdf', dpi=600)

        out = 'ccre_cis_combination_'+strtt+'_dim'+str(dim)+'_bidx'+str(b_idx)
        plt.savefig(out+'.pdf', dpi=600)
        plt.cla()
        plt.clf()
        plt.close()

        if np.sum(pie_chart) != 0:

            #print(pie_chart)

            #Pie chart
            plt.pie(pie_chart, labels=['no cCRE', 'has cCRE'], colors=['#D41159', '#1A85FF'])
            plt.title('Number of cliques '+str(count))
            out = 'ccre_cis_pie_chart_combination_'+strtt+'_dim'+str(dim)+'_bidx'+str(b_idx)

            #plt.show()

            #exit()

            plt.savefig(out+'.pdf', dpi=600)
            plt.cla()
            plt.clf()
            plt.close()


        ##plt.show()
        #
        #

    #exit()


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
#
#
#
#
#
#
#    print('DOINT B_IDX', b_idx)
#
#    anchor_thresh = b_idx + 1
#
#    graph_dict = dict()
#
#    common_anchors = dict()
#    for tt in tuples:
#        common_anchors[tt] = dict()
#
#    # Initialize set intersection counts to 0
#    s = pd.Series(np.zeros(len(tuples), dtype=int)\
#        , index=index)
#
#
#
#for e_idx, exp in enumerate(labels):
#
#
#
#
#
#
#
#    print(exp)
#    #if e_idx != 3:
#    #    continue
#    #exit()
#
#    ccre_type_count = dict()
#
#    # cfile
#    cfile = '../cooler_files/'+exp + '.cool'
#
#    cf = cooler.Cooler(cfile)
#
#    chrom_start_bins_list = []
#    chrom_names = []
#
#    n_chrom = 0
#    
#    for ch in cf.chromnames:
#        
#        chrom_names.append(ch)
#        beg, end = cf.extent(ch)
#        chrom_start_bins_list.append(beg)
#        n_chrom += 1
#
#    chrom_start_bins_list.append(cf.extent(cf.chromnames[-1])[1])
#    chrom_start_bins = np.array(chrom_start_bins_list)
#
#
#    for idx in range(1, 2):
#
#        print(exp, idx)
#
#        # Initialize set count
#        s = pd.Series(np.zeros(len(tuples), dtype=int)\
#                    , index=index)
#
#        cycle_info = []
#
#        cis_cycles = []
#
#        # Look at scaled
#        #target = '../cooler_results/PH_results/'+exp+'_gdist'+str(idx+1)+'_ALL_'
#        target = '../cooler_results/PH_results/'+exp+'_gdist'+str(idx+1)+'_ALL_scaled_'
#
#        # H1 shortened cycles (and have been checked for connectedness)
#        cyc_file = target + 'minimal_V_birth_H1.txt'
#
#        cycs = open(cyc_file, 'r')
#
#        progress = 0
#
#        has_cre = 0
#        no_cre = 0
#
#        for line in cycs:
#
#            print(progress, end='\r')
#            progress += 1
#
#            line = line.split(',')
#            line = line[:-1]
#            line = [int(x) for x in line]
#
#            # Since this is after connectedness check
#            # The cycle is stored as v_i, v_{i+1}
#            # Append the first element to end to complete the cycle
#            line.append(line[0])
#            line = np.array(line)
#
#            this_g_dist = np.abs(line[1:] - line[:-1])
#
#            # multiple anchors?
#            # Simply assume that any edge with more than one 2 genomic distance is an anchor
#            # Later plot the statistics
#            all_possible_anchors = np.argwhere(this_g_dist > 2).flatten()
#
#            ## OR pick the maximum one
#            #all_possible_anchors = np.argmax(this_g_dist).flatten()
#
#            #print(all_possible_anchors)
#
#            #print(len(all_possible_anchors))
#
#            flag = 0
#            for anchor in all_possible_anchors:
#
#                anchor_b1 = line[anchor]
#                anchor_b2 = line[anchor+1]
#
#                chrr_1, binn_1 = get_chrom(anchor_b1, chrom_start_bins_list, chrom_names)
#                chrr_2, binn_2 = get_chrom(anchor_b2, chrom_start_bins_list, chrom_names)
#                if chrr_1 != chrr_2:
#                    continue
#
#                this_ccre_dict = ccre_dict['chr'+chrr_1]
#                start = max(binn_1 - 1, 0)
#                end = binn_1 + 1
#
#                cre_set = []
#                for key in this_ccre_dict:
#
#                    ccre_start = key[0]
#                    ccre_end = key[1]
#
#                    if (ccre_start >= start and ccre_start <= end)\
#                        or (ccre_end >= start and ccre_end <= end):
#
#                        #print(this_ccre_dict[key])
#
#                        for cc in this_ccre_dict[key]:
#                            c_idx = ccre_types.index(cc)
#                            #print(c_idx)
#                            cre_set.append(c_idx)
#
#                cre_set = np.array(list((frozenset(cre_set))), dtype=int)
#                this_tuple = np.zeros((6,), dtype=int)
#                this_tuple[cre_set] = 1
#                this_tuple = tuple(this_tuple)
#                s[this_tuple] += 1
#
#
#                for bb in range(binn_1-1, binn_1+2):
#
#                    for rr in this_ccre_dict:
#                        start = rr[0]
#                        end = rr[1]
#
#                        if not bb < start and not bb > end:
#                            for cc in this_ccre_dict[rr]:
#                                if cc not in ccre_type_count:
#                                    ccre_type_count[cc] = 1
#                                else:
#                                    ccre_type_count[cc] += 1
#                            has_cre += 1
#                            flag = 1
#                            break
#
#                    if flag:
#                        break
#
#                if flag:
#                    break
#
#
#                this_ccre_dict = ccre_dict['chr'+chrr_2]
#                start = max(binn_2 - 1, 0)
#                end = binn_2 + 1
#
#                cre_set = []
#                for key in this_ccre_dict:
#
#                    ccre_start = key[0]
#                    ccre_end = key[1]
#
#                    if (ccre_start >= start and ccre_start <= end)\
#                        or (ccre_end >= start and ccre_end <= end):
#
#                        #print(this_ccre_dict[key])
#
#                        for cc in this_ccre_dict[key]:
#                            c_idx = ccre_types.index(cc)
#                            #print(c_idx)
#                            cre_set.append(c_idx)
#
#                cre_set = np.array(list((frozenset(cre_set))), dtype=int)
#                this_tuple = np.zeros((6,), dtype=int)
#                this_tuple[cre_set] = 1
#                this_tuple = tuple(this_tuple)
#                s[this_tuple] += 1
#
#                for bb in range(binn_2-1, binn_2+2):
#
#                    for rr in this_ccre_dict:
#                        start = rr[0]
#                        end = rr[1]
#
#                        if not bb < start and not bb > end:
#                            for cc in this_ccre_dict[rr]:
#                                if cc not in ccre_type_count:
#                                    ccre_type_count[cc] = 1
#                                else:
#                                    ccre_type_count[cc] += 1
#                            has_cre += 1
#                            flag = 1
#                            break
#
#                    if flag:
#                        break
#
#            
#
#            if not flag:
#                no_cre += 1
#
#
#        for tt in tuples:
#            s[tt] = math.log1p(s[tt])
#
#        upsetplot(s)
#        plt.show()
#        plt.cla()
#        plt.clf()
#        plt.close()
#
#        print(no_cre, has_cre)
#
#        print(ccre_type_count)
#
#        input('w')
#
#
#
##                # Get ATAC peaks in the anchor bins
##                atac_peaks_1 = hf2.get_narrow_range(atac_narrow_file, 'chr'+chrr_1\
##                                                , (anchor_b1-1)*reso, (anchor_b1+2)*reso)
##
##                for ll in atac_peaks_1:
##                    for entry in ll[:-1]:
##                        test_bed_file.write(entry+'\t')
##                    test_bed_file.write(ll[-1])
##
##                atac_peaks_2 = hf2.get_narrow_range(atac_narrow_file, 'chr'+chrr_2\
##                                                , (anchor_b2-1)*reso, (anchor_b2+2)*reso)
##
##                for ll in atac_peaks_2:
##                    for entry in ll[:-1]:
##                        test_bed_file.write(entry+'\t')
##                    test_bed_file.write(ll[-1])
##                #print(binn_1, binn_2)
##                #print(atac_peaks_1, atac_peaks_2)
##
##                if len(atac_peaks_1) and len(atac_peaks_2):
##                    counts[2] += 1
##                elif not len(atac_peaks_1) and not len(atac_peaks_2):
##                    counts[0] += 1
##                else:
##                    counts[1] += 1
##
##                print(counts)
##
##        test_bed_file.close()
##        exit()
##        continue
##
###            iddx = np.argmax(this_g_dist)
###            maxx = np.max(this_g_dist)
###
###            #if maxx < 500:
###            #    continue
###
###            #if iddx > 1:
###            #    print('Multiple maximal anchors with g dist ', maxx)
###
###            anchor_b1 = line[iddx]
###            anchor_b2 = line[iddx+1]
###
###            chrr_1, binn_1 = get_chrom(anchor_b1, chrom_start_bins_list, chrom_names)
###            chrr_2, binn_2 = get_chrom(anchor_b2, chrom_start_bins_list, chrom_names)
###
###            print(chrr_1, binn_1, chrr_2, binn_2)
###
###            ###########################################
###            # Narrow peaks
###            ###########################################
###            fig, axs = plt.subplots(3, 2, figsize=(10, 8))
###
###            if chrr_1 != chrr_2:
###                fig.suptitle('trans anchor, genomic dist '+str(maxx))
###            else:
###                fig.suptitle('cis anchor, genomic dist '+str(maxx))
###
###            row = 0
###            col = 0
###            axs[row, col].set_title('anchor 1, chrom '+chrr_1+', bin '+str(binn_1))
###
###            for idx, bfile in enumerate(fnames):
###
###                hf2.plot_narrow_range(axs[row,col], bfile\
###                                     , 'chr'+chrr_1\
###                                     , (binn_1-2)*reso\
###                                     , (binn_1+2)*reso\
###                                     )
###
###                y_lims = axs[row, col].get_ylim()
###                axs[row, col].scatter(binn_1*reso/1000, y_lims[1], color='red')
###
###                axs[row, col].set_ylabel(flabel[idx])
###
###                row += 1
###
###            row = 0
###            col = 1
###            axs[row, col].set_title('anchor 2, chrom '+chrr_2+', bin '+str(binn_2))
###
###            for idx, bfile in enumerate(fnames):
###
###                hf2.plot_narrow_range(axs[row,col], bfile\
###                                     , 'chr'+chrr_2\
###                                     , (binn_2-2)*reso\
###                                     , (binn_2+2)*reso\
###                                     )
###
###                y_lims = axs[row, col].get_ylim()
###                axs[row, col].scatter(binn_2*reso/1000, y_lims[1], color='red')
###
###
###                row += 1
###
###            
###            plt.tight_layout()
###            plt.show()
###            plt.cla()
###            plt.clf()
###            plt.close()
###
###            #exit()
###
###
###            ###########################################
###            # BigWig peaks
###            ###########################################
###            fig, axs = plt.subplots(3, 2, figsize=(10, 8))
###
###            if chrr_1 != chrr_2:
###                fig.suptitle('trans anchor, genomic dist '+str(maxx))
###            else:
###                fig.suptitle('cis anchor, genomic dist '+str(maxx))
###
###            row = 0
###            col = 0
###            axs[row, col].set_title('anchor 1, chrom '+chrr_1+', bin '+str(binn_1))
###
###            for idx, bfile in enumerate(bigwig_fnames):
###
###                hf2.plot_bigwig_range(axs[row,col], bfile\
###                                     , 'chr'+chrr_1\
###                                     , (binn_1-2)*reso\
###                                     , (binn_1+2)*reso\
###                                     )
###
###                y_lims = axs[row, col].get_ylim()
###                axs[row, col].scatter(binn_1*reso/1000, y_lims[1], color='red')
###
###                axs[row, col].set_ylabel(bigwig_flabel[idx])
###
###                row += 1
###
###
###            row = 0
###            col = 1
###            axs[row, col].set_title('anchor 2, chrom '+chrr_2+', bin '+str(binn_2))
###
###            for idx, bfile in enumerate(bigwig_fnames):
###
###                hf2.plot_bigwig_range(axs[row,col], bfile\
###                                     , 'chr'+chrr_2\
###                                     , (binn_2-2)*reso\
###                                     , (binn_2+2)*reso\
###                                     )
###
###                y_lims = axs[row, col].get_ylim()
###                axs[row, col].scatter(binn_2*reso/1000, y_lims[1], color='red')
###
###                axs[row, col].set_ylabel(bigwig_flabel[idx])
###
###                row += 1
###
###
###            plt.tight_layout()
###
###            plt.show()
###
###
###
###        exit()
###    
###
###
###
###
###
###
###
###
###
