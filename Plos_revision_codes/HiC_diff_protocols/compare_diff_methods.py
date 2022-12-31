import pydory as dory
import check_connected_helper as cch
import helper_functions as hf 
import helper_funcs_v2 as hf2
import networkx as nx
import matplotlib.pyplot as plt
import os
import numpy as np


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
    hiccups_data.append(hf2.get_hiccup_peaks(ff, reso))



dim = 1

PH_reso = reso


colors = [\
        '#E69F00'\
        ,'#56B4E9'\
        ,'#009E73'\
        ,'#D55E00'\
        ,'#CC79A7'\
        ]

row_size = 3
col_size = 2

if reso == 10000:
    threshs = [
            0.00914071\
            ,0.0083232\
            ,0.01087386\
            ,0.0090625\
            ,0.01120925\
            ]
    
    birth_threshs = [
            [0.00119952 ,0.00303552 ,0.00497419]\
            ,[0.00117106 ,0.00283449 ,0.00457273]\
            ,[0.00059415 ,0.00314951 ,0.00558935]\
            ,[0.00091857 ,0.00278725 ,0.00474511]\
            ,[0.00249192 ,0.00455864 ,0.00669773]\
            ]
elif reso == 5000:
    threshs = [
            0.01622359\
            ,0.01561109\
            ,0.01703869\
            ,0.01615373\
            ,0.02395455\
            ]
    
    birth_threshs = [
            [0.00240757, 0.00534985, 0.00879996]\
            ,[0.00276559 ,0.00545233 ,0.00868856]\
            ,[0.00073935 ,0.00449164 ,0.00836386]\
            ,[0.00185776 ,0.00480261 ,0.00835398]\
            ,[0.00710244 ,0.01101183 ,0.01509577]
            ]


graph_files = [\
             '../fa_dsg_ddel_graph.gml'\
            ,'../fa_dsg_dpn_graph.gml'\
            ,'../fa_dsg_mnase_graph.gml'\
            ,'../fa_dsg_ddel_dpn_graph.gml'\
            ,'../fa_dpn_graph.gml'\
                ]

all_graphs = []

for idx, gfile in enumerate(graph_files):

    if os.path.isfile(gfile):
        print('Loading graph', idx)
        all_graphs.append(nx.read_gml(gfile))
        continue


    edge_file = edge_file_name[idx]

    print('Making graph', idx)

    # Make a graph on source up to 200
    global_G = nx.Graph()
    
    edges_data = open(edge_file, 'r')

    this_thresh = threshs[idx]
    
    for line in edges_data:
    
        line = line.strip('\n')
        line = line.split(',')
    
        v1 = int(line[0])*reso
        v2 = int(line[1])*reso
        weight = float(line[2])

        if weight > this_thresh:
            continue
    
        global_G.add_edge(v1, v2, weight=weight)
    
    edges_data.close()

    # save graph
    nx.write_gml(global_G, gfile)

    all_graphs.append(global_G)






for b_idx in range(3):

    #ax = plt.gca()


    fig, axs = plt.subplots(row_size\
                    , col_size\
                    , figsize=(12, 8)\
                    , sharex=True\
                    , sharey=True\
                    )

    row = 0
    col = 0

    for idx, source in enumerate(edge_file_name):

        main_target = targets[idx]

        target = main_target + 'threshidx' + str(b_idx) +'_'


        ### using greedy-shortening set
        #cyc_file = target + 'minimal_V_birth_H'+str(dim)+'.txt'
        
        # using smoothed
        cyc_file = target+'smoothen_H'+str(dim)+'.txt'
        
        # Get cycs as list of edges
        cycs = hf2.get_cycs_as_edge_lists(cyc_file)
        
        hf2.make_line_graph_plot(axs[row, col]\
                        , cycs\
                        , colors[idx]\
                        , 0\
                        , scale_reso=PH_reso)
        
        axs[row, col].set_title(labels[idx])

        hiccup_x = hiccups_data[idx][:, 0]
        hiccup_y = hiccups_data[idx][:, 1]

        axs[row, col].scatter(hiccup_x, hiccup_y, color='red', s=40, marker='+', zorder=10000)


        col += 1
        if col == col_size:
            col = 0
            row += 1
        
    plt.show()

    plt.cla()
    plt.clf()
    plt.close()










