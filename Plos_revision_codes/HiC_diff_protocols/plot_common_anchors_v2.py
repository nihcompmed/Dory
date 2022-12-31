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

dim = 1


ccre_dict = pickle.load(open('../Biomarkers/ccre_reso'+str(reso)+'.p', 'rb'))

# Dimension 1 cis common anchors
all_common_anchors = pickle.load(open('cis_common_anchors.p', 'rb'))


for exp in labels:

    for b_idx in range(2):

        all_chroms = []
        for key in all_common_anchors[b_idx]:
            print(key)

        exit()
    
        print(b_idx)
    
        # Initialize set count
        s = pd.Series(np.zeros(len(tuples), dtype=int)\
                    , index=index)
    
        info = pickle.load(open('common_cis_anchors_dim'+str(dim)+'_bidx'+str(b_idx)+'.p', 'rb'))
    
        all_nearest = []
    
        count = 0
    
        bed_file_5kb = 'cis_atac_5kb_dim'+str(dim)+'_bidx'+str(b_idx)+'.bed'
        bed_file_10kb = 'cis_atac_10kb_dim'+str(dim)+'_bidx'+str(b_idx)+'.bed'
    
        bed_5kb = open(bed_file_5kb, 'w')
        bed_10kb = open(bed_file_10kb, 'w')
    
        ccre_type_count = dict()
    
        #count_no_cre = 0
        #count_has_cre = 0
    
        count_anchors = 0
    
        # First get nearest ATAC peak
        for chrom in info:
            print('DOING CHROM', chrom)
            count += 1
            anchors = info[chrom]
    
            print(anchors)
            exit()
    
            this_ccre_dict = ccre_dict['chr'+chrom]
    
            for anchor in anchors:
    
                count_anchors += 1
    
                start = max(anchor - 1, 0)
                end = anchor + 1
    
                has_cre = 0
    
                cre_set = []
    
                for key in this_ccre_dict:
    
                    ccre_start = key[0]
                    ccre_end = key[1]
    
                    if (ccre_start >= start and ccre_start <= end)\
                        or (ccre_end >= start and ccre_end <= end):
    
                        #print(this_ccre_dict[key])
    
                        for cc in this_ccre_dict[key]:
                            c_idx = ccre_types.index(cc)
                            #print(c_idx)
                            cre_set.append(c_idx)
                            #if cc not in ccre_type_count:
                            #    ccre_type_count[cc] = 1
                            #else:
                            #    ccre_type_count[cc] += 1
    
                cre_set = np.array(list((frozenset(cre_set))), dtype=int)
    
                this_tuple = np.zeros((6,), dtype=int)
    
                #print(cre_set)
    
                this_tuple[cre_set] = 1
    
                this_tuple = tuple(this_tuple)
    
                s[this_tuple] += 1
    
    
    
        #print(count_no_cre, count_has_cre)
        #print(ccre_type_count)
    
        for tt in tuples:
            s[tt] = s[tt]/count_anchors*100
    
        upsetplot(s)
        plt.show()
        plt.cla()
        plt.clf()
        plt.close()
        #exit()
    
    
    
    
    
    
    
    
    
