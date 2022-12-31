import pydory as dory
import check_connected_helper as cch
import helper_functions as hf 
import networkx as nx
import matplotlib.pyplot as plt
import os

reso = 5000
reso = 10000

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



filetype = 2
dim = 1
threads = 4

for idx, source in enumerate(edge_file_name):

    main_target = targets[idx]
    thresh = threshs[idx]

    for b_idx, birth_thresh in enumerate(birth_threshs[idx]):

        target = main_target + 'threshidx' + str(b_idx) +'_'

        print('Computing PH...')

        # Birth thresh is eps
        dory.compute_PH(source, 0, thresh, filetype, threads, target, dim, 1, 1, birth_thresh, 1)

        print('Done. Modifying for connectedness...')

        # Modify for connectedness
        dim = 1
        cyc_file = target + 'minimal_V_birth_H'+str(dim)+'.txt'

        if not os.path.isfile(cyc_file):
            continue

        cch.modify_minimal(cyc_file, target, dim)

        print('Done. Smoothing...')

        # Smoothen 
        # 1. Get graph
        filetype = 2
        G = hf.get_edge_graph(source, filetype, birth_thresh)
        
        # 2. Smoothen with trivial triangles
        dim = 1
        # H1 boundary
        ff = open(target+'minimal_V_birth_H'+str(dim)+'.txt', 'r')
        
        # New file that will have smoothed boundary
        gg = open(target+'smoothen_H'+str(dim)+'.txt', 'w')
        
        count = 0
        for line in ff:
            print(count, end='\r')
            count += 1
        
            line = line.split(',')
            line = line[:-1]
            line = [int(x) for x in line]
        
            line = hf.smoothen_trivial_triangles(line, G)
        
            for x in line:
                gg.write(str(x)+',')
            gg.write('\n')
        
        ff.close()
        gg.close()


        print('Done.')

