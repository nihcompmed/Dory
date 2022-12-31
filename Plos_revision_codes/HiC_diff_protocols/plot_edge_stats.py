import matplotlib.pyplot as plt
import numpy as np
import math
import pickle


labels = [\
         'fa_dsg_ddel_5kb'\
        ,'fa_dsg_dpn_5kb'\
        ,'fa_dsg_ddel_dpn_5kb'\
        ,'fa_dsg_mnase_5kb'\
        ,'fa_dpn_5kb'\
        ,'fa_dsg_ddel_1kb'\
        ,'fa_dsg_dpn_1kb'\
        ,'fa_dsg_ddel_dpn_1kb'\
        ,'fa_dsg_mnase_1kb'\
        ,'fa_dpn_1kb'\
        ]

exps = [
         'fa_dsg_ddel'\
        ,'fa_dsg_dpn'\
        ,'fa_dsg_ddel_dpn'\
        ,'fa_dsg_mnase'\
        ,'fa_dpn'\
        ]



reso_num = [\
          '5kb'\
        , '5kb'\
        , '5kb'\
        , '5kb'\
        , '5kb'\
        , '1kb'\
        , '1kb'\
        , '1kb'\
        , '1kb'\
        , '1kb'\
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

g_dists = list(range(1, 11))

percentiles = [5, 50, 95]

save_dir = '../cooler_results'

all_stats = pickle.load(open('ice_norm_stats.p', 'rb'))

reso_num = dict()
reso_num['1kb']=1000
reso_num['5kb']=5000

for reso in ['1kb', '5kb']:

    for chrom in chroms:

        fig, axs = plt.subplots(1, 3, sharey=True, figsize=(10, 8))

        for p_idx, ptile in enumerate(percentiles):

        
            for exp in exps:

                label = exp+'_'+reso
        
                info = all_stats[label][reso]

                xx = []
                yy = []

                for g_dist in g_dists:

                    xx.append(g_dist*reso_num[reso]/1000)
                    yy.append(info[chrom][g_dist][ptile])

                axs[p_idx].plot(xx, yy, alpha=0.9, label=exp)

        axs[0].legend()
        axs[0].set_xlabel('distance along chrom in kb')
        axs[0].set_ylabel('spatial estimate from ICE norm counts')

        axs[0].set_title('chrom chr'+chrom+', reso '+reso+', ptile '+str(percentiles[0]))
        axs[1].set_title('ptile '+str(percentiles[1]))
        axs[2].set_title('ptile '+str(percentiles[2]))

        plt.show()
        exit()












