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
from sklearn.neighbors import KernelDensity
import scipy
import matplotlib as mpl


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
         'fa_dsg_ddel_5kb'\
        ,'fa_dsg_dpn_5kb'\
        ,'fa_dsg_ddel_dpn_5kb'\
        ,'fa_dsg_mnase_5kb'\
        ,'fa_dpn_5kb'\
        ]

scale= [
        1\
        ,1\
        ,0.92\
        ,0.75\
        ,1.15\
        ]

shift = [
        0\
        ,0\
        ,0\
        ,0\
        ,0\
        ]

mu = [\
        14.029\
        ,15.822\
        ,10.418\
        ,5.948\
        ,45.507\
    ]

sigma = [\
        3.9248\
        ,4.4154\
        ,2.992\
        ,1.9266\
        ,12.3218\
        ]

posshift = [\
        3.342\
        ,3.331\
        ,3.223\
        ,2.644\
        ,3.3912\
            ]

ptiles = np.arange(5, 100, 5)

for b_idx in range(1, 2):

    dim = 1

    
    #for t1, t2 in itertools.combinations(labels, 2):
    for e_idx, exp in enumerate(labels):


        if e_idx != 3:
            continue
    
        t1 = '../cooler_results/PH_results/'+exp+'_gdist'+str(b_idx+1)+'_ALL_scaled_'
        t2 = '../cooler_results/PH_results/'+exp+'_gdist'+str(b_idx+1)+'_ALL_'

        pers1_file = t1+'H'+str(dim)+'_pers_data.txt'
        pers2_file = t2+'H'+str(dim)+'_pers_data.txt'
    
        pers1 = np.loadtxt(pers1_file, delimiter=',')
        pers2 = np.loadtxt(pers2_file, delimiter=',')

        # Compare 2D distribution of dead features first
        dead1 = np.argwhere(pers1[:,1]!=-1).flatten()
        pers1 = pers1[dead1]
        diff1 = pers1[:,1] - pers1[:,0]

        dead2 = np.argwhere(pers2[:,1]!=-1).flatten()
        pers2 = pers2[dead2]


        #################################
        ## Plotting log norm binning of PD
        #################################
        #fig, axs = plt.subplots(1, 2)
        #print('plotting')

        #axs[0].hist2d(pers1[:,0], pers1[:,1], bins=100\
        #            , norm=mpl.colors.LogNorm(), cmap='gnuplot')

        #axs[1].hist2d(pers2[:,0], pers2[:,1], bins=100\
        #            , norm=mpl.colors.LogNorm(), cmap='gnuplot')

        #plt.show()

        #exit()


        ################################
        # Scale pers2 by the same scaling as done for experiment
        ################################
        pers2 = hf2.scale_hic_estimate_array(pers2\
                                            , mu[e_idx], sigma[e_idx], posshift[e_idx]\
                                            , scale[e_idx], shift[e_idx]\
                                            )


        diff2 = pers2[:,1] - pers2[:,0]

        print('number of dead features in data 1', pers1.shape\
              ,'number of dead features in data 2', pers2.shape)
        
        #plt.scatter(pers1[:, 0], pers1[:, 1], color='blue', alpha=0.4)
        #plt.scatter(pers2[:, 0], pers2[:, 1], color='red', alpha=0.4)
        #plt.show()
        #exit()

        #plt.hist(diff1, bins=100, color='red', alpha=0.5)
        #plt.hist(diff2, bins=100, color='blue', alpha=0.5)
        #plt.show()

        #print(diff1, diff2)

        #exit()

        ##KDE APPROACH

        #xx_min = min(np.amin(pers1[:,0]), np.amin(pers2[:,0]))
        #xx_max = max(np.amax(pers1[:,0]), np.amax(pers2[:,0]))

        #yy_min = min(np.amin(pers1[:,1]), np.amin(pers2[:,1]))
        #yy_max = max(np.amax(pers1[:,1]), np.amax(pers2[:,1]))
        #print(xx_min, xx_max, yy_min, yy_max)
        ##exit()

        ## make grid, for now 100 by 100
        #xbins = 100j
        #ybins = 100j
        #bandwidth = max(np.amax(diff1), np.amax(diff2))/50

        #xx, yy = np.mgrid[xx_min:xx_max:xbins, 
        #              yy_min:yy_max:ybins]

        #xy_sample = np.vstack([yy.ravel(), xx.ravel()]).T

        #x = pers1[:, 0]
        #y = pers1[:, 1]
        #xy_train1  = np.vstack([y, x]).T

        #kde_skl = KernelDensity(bandwidth=bandwidth)
        #kde_skl.fit(xy_train1)
        #z1 = np.exp(kde_skl.score_samples(xy_sample))

        #z1 = np.reshape(z1, xx.shape)

        #print(z1.shape)

        #plt.imshow(z1, alpha=0.3, origin='lower')
        ##plt.scatter(x, y)
        #plt.show()

        #x = pers2[:, 0]
        #y = pers2[:, 1]
        #xy_train2  = np.vstack([y, x]).T

        #kde_skl = KernelDensity(bandwidth=bandwidth)
        #kde_skl.fit(xy_train2)
        #z2 = np.exp(kde_skl.score_samples(xy_sample))

        #print(np.sum(z1.flatten()), np.sum(z1.flatten()))

        ##print(z1.shape, z2.shape)

        ##exit()
        #print(scipy.stats.entropy(z1.flatten(), z2.flatten()))



        ### Compute KDE
        ##xx, yy, zz = hf2.kde2D(pers1[:,0], pers1[:,1], np.amax(diff1)/100)
        ##xx, yy, zz = hf2.kde2D(pers1[:,0], pers1[:,1], np.amax(diff1)/100)
        ##print(xx.shape, yy.shape, zz.shape)
        #exit()



        # Gradation by persistence
        # Gradation by PD that has higher persistence
        maxx1 = np.amax(diff1)
        maxx2 = np.amax(diff2)
        if maxx1 >= maxx2:
            pers_thresholds = np.percentile(diff1, ptiles)
        else:
            pers_thresholds = np.percentile(diff2, ptiles)

        for p_idx, pers_thresh in enumerate(pers_thresholds[:-1]):
            idxs1 = np.argwhere(diff1 >= pers_thresh).flatten()
            this_pers1 = pers1[idxs1]

            this_diff = this_pers1[:, 1] - this_pers1[:, 0]
            idxs1 = np.argwhere(this_diff < pers_thresholds[p_idx+1]).flatten()
            this_pers1 = this_pers1[idxs1]



            idxs2 = np.argwhere(diff2 >= pers_thresh).flatten()
            this_pers2 = pers2[idxs2]
        
            this_diff = this_pers2[:, 1] - this_pers2[:, 0]
            idxs2 = np.argwhere(this_diff < pers_thresholds[p_idx+1]).flatten()
            this_pers2 = this_pers2[idxs2]


            print(hf2.KL_symm(this_pers1, this_pers2))

            input('w')

    
        exit()









