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

xticks = [
         'fa+dsg, ddel'\
        ,'fa+dsg, dpn'\
        ,'fa+dsg, ddel_dpn'\
        ,'fa+dsg, mnase'\
        ,'fa, dpn'\
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

ptile = 95

for reso in ['1kb', '5kb']:

    for g_dist in g_dists:
         
        xx = []
        yy = []
        vall = g_dist*reso_num[reso]/1000
        
        for e_idx, exp in enumerate(exps):

            label = exp+'_'+reso
        
            info = all_stats[label][reso]

            this_exp_yy = []

            for chrom in chroms:

                xx.append(e_idx)
                yy.append(info[chrom][g_dist][ptile])

                this_exp_yy.append(info[chrom][g_dist][ptile])

            this_exp_yy = np.array(this_exp_yy)

            print(exp, g_dist, reso, np.max(this_exp_yy))

        plt.scatter(xx, yy, color='blue', alpha=0.6)

        plt.xticks(ticks=[0,1,2,3,4], labels=xticks, rotation=10)

        plt.title('reso '+reso+', separation along chromosome in kb '+ str(vall))

        plt.ylabel('spatial estimate from cooler default balance')


        plt.savefig('../cooler_results/reso'+reso+'_bindist'+str(vall)+'_ptile'+str(ptile)+'.pdf', dpi=600)

        plt.show()

        plt.cla()
        plt.clf()
        plt.close()



#    for chrom in chroms:
#
#        fig, axs = plt.subplots(1, 3, sharey=True, figsize=(10, 8))
#
#        for p_idx, ptile in enumerate(percentiles):
#
#        
#
#                label = exp+'_'+reso
#        
#                info = all_stats[label][reso]
#
#
#
#
#                axs[p_idx].plot(xx, yy, alpha=0.9, label=exp)
#
#        axs[0].legend()
#        axs[0].set_xlabel('distance along chrom in kb')
#        axs[0].set_ylabel('spatial estimate from ICE norm counts')
#
#        axs[0].set_title('chrom chr'+chrom+', reso '+reso+', ptile '+str(percentiles[0]))
#        axs[1].set_title('ptile '+str(percentiles[1]))
#        axs[2].set_title('ptile '+str(percentiles[2]))
#
#        plt.show()
#        exit()
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
