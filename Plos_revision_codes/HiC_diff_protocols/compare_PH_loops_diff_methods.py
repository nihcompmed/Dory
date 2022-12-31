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



for b_idx in range(3):

    #ax = plt.gca()


    #fig, axs = plt.subplots(row_size\
    #                , col_size\
    #                , figsize=(12, 8)\
    #                , sharex=True\
    #                , sharey=True\
    #                )

    #row = 0
    #col = 0

    cycs_all_exps = []
    genomic_dists_all_exps = []

    for idx, source in enumerate(edge_file_name):

        main_target = targets[idx]

        target = main_target + 'threshidx' + str(b_idx) +'_'


        ### using greedy-shortening set
        #cyc_file = target + 'minimal_V_birth_H'+str(dim)+'.txt'
        
        # using smoothed
        cyc_file = target+'smoothen_H'+str(dim)+'.txt'
        
        # Get cycs as list of edges
        # These are on the scale of original resolution
        cycs = hf2.get_cycs_as_edge_lists(cyc_file)


        max_anchors = []

        # Get max genomic dists along boundary
        cycs_genomic_dists = []
        for cyc in cycs:
            #print(cyc)
            cyc = np.array(cyc)
            ddist = np.abs(cyc[:,0] - cyc[:,1])

            maxx = np.argmax(ddist).flatten()
            this_max_anchors = cyc[maxx]
            max_anchors.append(this_max_anchors)

            vval = np.amax(ddist)
            cycs_genomic_dists.append(vval)


        cycs_genomic_dists = np.array(cycs_genomic_dists)

        genomic_dists_all_exps.append(cycs_genomic_dists)

        cycs_all_exps.append(max_anchors)

        for max_anchor in max_anchors:

            print(max_anchor)
            exit()

        #print(cycs_all_exps)
        #exit()


        
        #hf2.make_line_graph_plot(axs[row, col]\
        #                , cycs\
        #                , colors[idx]\
        #                , 0\
        #                , scale_reso=PH_reso)
        #
        #axs[row, col].set_title(labels[idx])

        #col += 1
        #if col == col_size:
        #    col = 0
        #    row += 1


    n_exp = len(cycs_all_exps)


    for (exp1, exp2) in itertools.combinations(list(range(n_exp)), 2):

        print('doing', labels[exp1], ' and ', labels[exp2])

        cycs1 = cycs_all_exps[exp1]
        cycs2 = cycs_all_exps[exp2]

        n_cycs1 = len(cycs1)
        n_cycs2 = len(cycs2)

        dist_mat = np.ones((n_cycs1, n_cycs2))*(-1)

        llist = [list(range(n_cycs1)), list(range(n_cycs2))]

        for c1_idx in range(n_cycs1):

            for c2_idx in range(n_cycs2):

                print(c1_idx, c2_idx, end='\r')

                c1 = cycs1[c1_idx]
                c2 = cycs2[c2_idx]
                pair_cyc_val = 0

                lllist = [c1, c2]

                for e1, e2 in itertools.product(*lllist):
                    x1, x2 = e1
                    y1, y2 = e2
                    val = max(abs(x1 - y1), abs(x2 - y2))
                    pair_cyc_val = max(pair_cyc_val, val)
                
                dist_mat[c1_idx, c2_idx] = pair_cyc_val

            
        minns1 = np.min(dist_mat, axis=0)
        plt.scatter(genomic_dists_all_exps[exp2]*reso/1000, minns1*reso/1000, alpha=0.5, color='blue', label=labels[exp2])

        minns2 = np.min(dist_mat, axis=1)
        plt.scatter(genomic_dists_all_exps[exp1]*reso/1000, minns2*reso/1000, alpha=0.5, color='red', label=labels[exp1])

        plt.legend()

        #minns = list(minns1) + list(minns2)

        #sns.distplot(minns)

        plt.xlabel('genomic dist')
        plt.ylabel('nearest anchor dist')

        plt.xscale('log')
        plt.yscale('log')

        #plt.show()
        plt.savefig(labels[exp1]+'_'+labels[exp2]+'_threshidx'+str(b_idx)+'_scatter.pdf', dpi=600)

        plt.cla()
        plt.clf()
        plt.close()

        plt.scatter(genomic_dists_all_exps[exp2]*reso/1000, minns1/genomic_dists_all_exps[exp2]\
                    , alpha=0.5, color='blue', label=labels[exp2])
        plt.scatter(genomic_dists_all_exps[exp1]*reso/1000, minns2/genomic_dists_all_exps[exp1]\
                    , alpha=0.5, color='red', label=labels[exp1])
        plt.legend()
        plt.xlabel('genomic dist')
        plt.ylabel('nearest anchor dist/genomic dist')

        plt.xscale('log')
        plt.yscale('log')

        plt.savefig(labels[exp1]+'_'+labels[exp2]+'_threshidx'+str(b_idx)+'_scatter_ratio.pdf', dpi=600)

        plt.cla()
        plt.clf()
        plt.close()
        #plt.show()

        #exit()

        #break



        #minns = np.array(list(minns1) + list(minns2))

        #print(np.median(minns))

        #sns.distplot(minns, label=labels[exp1]+' vs. '+labels[exp2])

    #plt.show()
    #plt.cla()
    #plt.clf()

    
    #print(np.min(dist_mat, axis=0), np.min(dist_mat, axis=1))
    #plt.hist()
    #exit()
            
    
    #plt.show()

    #plt.cla()
    #plt.clf()
    #plt.close()










