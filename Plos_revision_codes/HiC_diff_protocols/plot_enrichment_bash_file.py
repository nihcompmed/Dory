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

ff = open('plot_enrichment.sh', 'w')

ff.write('#!/bin/sh\n')

dim = 1

tuples = [\
          (0,0,0,0,0)\
          ,(0,0,0,0,1)\
          ,(0,0,0,1,0)\
          ,(0,0,0,1,1)\
          ,(0,0,1,0,0)\
          ,(0,0,1,0,1)\
          ,(0,0,1,1,0)\
          ,(0,0,1,1,1)\
          ,(0,1,0,0,0)\
          ,(0,1,0,0,1)\
          ,(0,1,0,1,0)\
          ,(0,1,0,1,1)\
          ,(0,1,1,0,0)\
          ,(0,1,1,0,1)\
          ,(0,1,1,1,0)\
          ,(0,1,1,1,1)\
          ,(1,0,0,0,0)\
          ,(1,0,0,0,1)\
          ,(1,0,0,1,0)\
          ,(1,0,0,1,1)\
          ,(1,0,1,0,0)\
          ,(1,0,1,0,1)\
          ,(1,0,1,1,0)\
          ,(1,0,1,1,1)\
          ,(1,1,0,0,0)\
          ,(1,1,0,0,1)\
          ,(1,1,0,1,0)\
          ,(1,1,0,1,1)\
          ,(1,1,1,0,0)\
          ,(1,1,1,0,1)\
          ,(1,1,1,1,0)\
          ,(1,1,1,1,1)\
        ]


for b_idx in range(1,2):

    count = 0

    for tt in tuples:

        count += 1
        if count < 3:
            continue


        strtt = ''
        for t in tt:
            strtt = strtt+str(t)


        bed_file = 'cis_combination_'+strtt+'_atac_dim'+str(dim)+'_bidx'+str(b_idx)+'.bed'
        out = 'cis_combination_'+strtt+'_atac_dim'+str(dim)+'_bidx'+str(b_idx)
        

        for m_idx, marker_file in enumerate(bigwig_fnames):

            ff.write('computeMatrix reference-point --referencePoint center -b 1000 -a 1000 --verbose -p 8 -R '\
                        +bed_file+' -S '+marker_file+' -o '+out+'_'+bigwig_flabel[m_idx]+'.gz')

            ff.write('\n')

            ff.write('plotProfile -m '+out+'_'+\
                        bigwig_flabel[m_idx]+'.gz -out '+out+'_'+\
                        bigwig_flabel[m_idx]+'_MEAN.pdf --dpi 600 --averageType mean')

            ff.write('\n')

            ff.write('plotProfile -m '+out+'_'+\
                        bigwig_flabel[m_idx]+'.gz -out '+out+'_'+\
                        bigwig_flabel[m_idx]+'_MEDIAN.pdf --dpi 600 --averageType median')

            ff.write('\n')

            ff.write('plotProfile -m '+out+'_'+\
                        bigwig_flabel[m_idx]+'.gz -out '+out+'_'+\
                        bigwig_flabel[m_idx]+'_MIN.pdf --dpi 600 --averageType min')

            ff.write('\n')

            ff.write('plotProfile -m '+out+'_'+\
                        bigwig_flabel[m_idx]+'.gz -out '+out+'_'+\
                        bigwig_flabel[m_idx]+'_MAX.pdf --dpi 600 --averageType max')

            ff.write('\n')

            ff.write('plotProfile -m '+out+'_'+\
                        bigwig_flabel[m_idx]+'.gz -out '+out+'_'+\
                        bigwig_flabel[m_idx]+'_STD.pdf --dpi 600 --averageType std')

            ff.write('\n')

ff.close()
