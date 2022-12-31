import pydory as dory
import pickle
import numpy as np
import time
import check_connected_helper as cch
import os


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

chroms = [
         '1'\
        ,'2'\
        ,'3'\
        ,'4'\
        ,'5'\
        ,'6'\
        ,'7'\
        ,'8'\
        ,'9'\
        ,'10'\
        ,'11'\
        ,'12'\
        ,'13'\
        ,'14'\
        ,'15'\
        ,'16'\
        ,'17'\
        ,'18'\
        ,'19'\
        ,'20'\
        ,'21'\
        ,'22'\
        ,'X'\
        ,'Y'\
        ]


save_dir = '../cooler_results/PH_results/'

#PH_pars = pickle.load(open('../cooler_results/edge_files/PH_pars_info_5kb.p', 'rb'))

PH_threshs = [\
        32.47655\
        ,29.70\
        ,32.264\
        ,31.597\
        ,29.25137
        ]

births = [\
         [5.0989,10.689]\
        ,[5.07815,9.9911]\
        ,[4.3887,10.425]\
        ,[2.9765,11.479]\
        ,[6.6947,11.299]\
        ]

#print(PH_pars)

#exit()


#reso_num = dict()
#reso_num['1kb']=1000
#reso_num['5kb']=5000

#ptile = 95

PH_info_file = open('PH_computation_info.txt', 'w')

# PH parameters
filetype = 2

threads = 4


reso = 5000

for e_idx, exp in enumerate(labels):

    # source edge file for all

    source = '../cooler_results/edge_files/'+exp+'_edges_ALL_scaled.csv'

    #PH_thresh = PH_pars[exp]['PH_thresh']
    PH_thresh = PH_threshs[e_idx]

    for idx, birth_thresh in enumerate(births[e_idx]):

        # gdist is birth thresh idx
        target = '../cooler_results/PH_results/'+exp+'_gdist'+str(idx+1)+'_ALL_scaled_'

        print(exp, idx, birth_thresh, PH_thresh)

        dim = 2

        print(f'Computing PH dim {dim}...')
        start = time.time()

        # Birth thresh is eps
        dory.compute_PH(source, 0, PH_thresh, filetype, threads, target, dim, 1, 1, birth_thresh, 1)

        print('Done. Modifying for connectedness...')

        # Modify for connectedness
        dim = 1
        cyc_file = target + 'minimal_V_birth_H'+str(dim)+'.txt'

        if not os.path.isfile(cyc_file):
            continue

        cch.modify_minimal(cyc_file, target, dim)

        end = time.time()

        time_taken = round(end-start, 2)

        PH_info_file.write(f'{exp}, {reso}, genomic dist {idx+1}, PH thresh {PH_thresh}, birth thresh {birth_thresh}, time taken {time_taken} seconds\n')


PH_info_file.close()








