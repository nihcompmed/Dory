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

PH_pars = pickle.load(open('../cooler_results/edge_files/PH_pars_info_5kb.p', 'rb'))

print(PH_pars)

#exit()


#reso_num = dict()
#reso_num['1kb']=1000
#reso_num['5kb']=5000

#ptile = 95

PH_info_file = open('PH_computation_info.txt', 'w')

# PH parameters
filetype = 2

dim = 2
threads = 4

reso = 5000

for exp in PH_pars:

    # source edge file for all

    source = '../cooler_results/edge_files/'+exp+'_edges_ALL.csv'

    PH_thresh = PH_pars[exp]['PH_thresh']

    for idx, birth_thresh in enumerate(PH_pars[exp]['births']):

        # gdist is birth thresh idx
        target = '../cooler_results/PH_results/'+exp+'_gdist'+str(idx+1)+'_ALL_'

        print(exp, idx, birth_thresh, PH_thresh)

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







