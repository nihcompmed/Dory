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


reso = 5000

#reso = 10000

edge_file_name = [\
         '../fa_dsg_ddel_chr1_KR_reso'+str(reso)+'.csv'\
        ,'../fa_dsg_dpn_chr1_KR_reso'+str(reso)+'.csv'\
        ,'../fa_dsn_mnase_chr1_KR_reso'+str(reso)+'.csv'\
        ,'../fa_dsg_ddel_dpn_chr1_KR_reso'+str(reso)+'.csv'\
        ,'../fa_dpn_chr1_KR_reso'+str(reso)+'.csv'\
        ]


targets = [\
         '../fa_dsg_ddel_reso'+str(reso)+'_'\
        ,'../fa_dsg_dpn_reso'+str(reso)+'_'\
        ,'../fa_dsg_mnase_reso'+str(reso)+'_'\
        ,'../fa_dsg_ddel_dpn_reso'+str(reso)+'_'\
        ,'../fa_dpn_reso'+str(reso)+'_'\
        ]

labels = [\
        'FA+DSG, Ddel'\
        ,'FA+DSG, Dpn'\
        ,'FA+DSG, MNase'\
        ,'FA+DSG, Ddel-Dpn'\
        ,'FA, Dpn'\
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
    chrom = 'chr1'
    hiccups_data.append(hf2.get_hiccup_peaks(ff, reso, chrom))



b_idx = 1

dim = 1

for idx, source in enumerate(edge_file_name):

    main_target = targets[idx]

    target = main_target + 'threshidx' + str(b_idx) +'_'


    ### using greedy-shortening set
    #cyc_file = target + 'minimal_V_birth_H'+str(dim)+'.txt'
    
    # using smoothed
    cyc_file = target+'smoothen_H'+str(dim)+'.txt'
    
    # Get cycs as list of edges
    cycs = hf2.get_cycs_as_edge_lists(cyc_file)

    # Hiccup anchors
    hiccup_x = hiccups_data[idx][:, 0]
    hiccup_y = hiccups_data[idx][:, 1]

    #print(hiccups_data[idx][:, :2])

    #exit()

    hf2.get_cycs_hiccups_dist(cycs, hiccups_data[idx][:, :2], reso)



    exit()

    # For every loop, find nearest hiccup peak






    






