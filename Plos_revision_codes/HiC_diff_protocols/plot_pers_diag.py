import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import helper_funcs_v2 as hf2

def plot_pers_hist2d(pers, ax, title):

    births = pers[:,0]
    deaths = pers[:,1]

    undead_idxs = np.argwhere(deaths==-1).flatten()
    dead_idxs = np.argwhere(deaths!=-1).flatten()

    #deaths[undead_idxs] = undead_maxx

    plot_pairs = pers[dead_idxs]

    plot_pairs = plot_pairs/np.amax(plot_pairs)

    maxx = np.amax(plot_pairs)
    undead_maxx = 1.05*maxx

    #ax.hist2d(births, deaths, bins=100, norm=mpl.colors.LogNorm(), cmap='gnuplot2_r')
    #ax.hist2d(births[dead_idxs], deaths[dead_idxs], bins=100, norm=mpl.colors.LogNorm(), cmap='gnuplot')
    ax.hist2d(plot_pairs[:, 0], plot_pairs[:, 1], bins=100, norm=mpl.colors.LogNorm(), cmap='gnuplot')

    #plt.colorbar()

    # y=x
    ax.plot([0,maxx], [0, maxx], ls='--', color='black', lw=2, alpha=0.8)
    ax.set_xlabel('birth', fontsize=14)
    ax.set_ylabel('death', fontsize=14)

    # undead
    ax.plot([0,maxx], [undead_maxx, undead_maxx], ls='--', color='black', lw=2, alpha=0.8)

    ax.set_title(title)

    ax.set_xticks([])
    ax.set_yticks([])

    # Persistence bands
    delta_percent = 0.1
    diff = plot_pairs[:, 1] - plot_pairs[:, 0]

    maxx_pers = np.amax(diff)
    delta = delta_percent*maxx_pers

    up_diff = maxx_pers
    low_diff = up_diff - delta

    bands = []
    
    for ii in range(8):

        idxs = np.argwhere(diff <= up_diff).flatten()
        short_list = plot_pairs[idxs]

        short_list_diff = short_list[:, 1] - short_list[:, 0]
        idxs = np.argwhere(short_list_diff > low_diff).flatten()

        short_list = short_list[idxs]

        bands.append(short_list)

        #print(np.amax(short_list[:,1] - short_list[:,0]))

        #input('w')

        up_diff -= delta
        low_diff -= delta

    #plt.show()

    return bands



exps = [\
         'fa_dsg_ddel_5kb.cool'\
        ,'fa_dsg_dpn_5kb.cool'\
        ,'fa_dsg_ddel_dpn_5kb.cool'\
        ,'fa_dsg_mnase_5kb.cool'\
        ,'fa_dpn_5kb.cool'\
#        ,'fa_dsg_ddel_1kb.cool'\
#        ,'fa_dsg_dpn_1kb.cool'\
#        ,'fa_dsg_ddel_dpn_1kb.cool'\
#        ,'fa_dsg_mnase_1kb.cool'\
        #,'fa_dpn_1kb.cool'\
       ]


labels = [\
         'fa_dsg_ddel_5kb'\
        ,'fa_dsg_dpn_5kb'\
        ,'fa_dsg_ddel_dpn_5kb'\
        ,'fa_dsg_mnase_5kb'\
        ,'fa_dpn_5kb'\
#        ,'fa_dsg_ddel_1kb'\
#        ,'fa_dsg_dpn_1kb'\
#        ,'fa_dsg_ddel_dpn_1kb'\
#        ,'fa_dsg_mnase_1kb'\
#        ,'fa_dpn_1kb'\
        ]


reso_num = [\
          '5kb'\
        , '5kb'\
        , '5kb'\
        , '5kb'\
        , '5kb'\
#        , '1kb'\
#        , '1kb'\
#        , '1kb'\
#        , '1kb'\
#        , '1kb'\
        ]


#plt_labels = [\
#         'fa+dsg, ddel, 5kb'\
#        ,'fa+dsg, dpn, 5kb'\
#        ,'fa+dsg, ddel+dpn, 5kb'\
#        ,'fa+dsg, mnase, 5kb'\
#        ,'fa+dsg, mnase, 5kb'\
#        ]

plt_labels = [\
         '(2) FA+DSG, DdeI, 5kb'\
        ,'(3) FA+DSG, DpnII, 5kb'\
        ,'(4) FA+DSG, DdeI+DpnII, 5kb'\
        ,'(5) FA+DSG, Mnase, 5kb'\
        ,'(1) FA, DpnII, 5kb'\
        ,'FA+DSG, DdeI, 1kb'\
        ,'FA+DSG, DpnII, 1kb'\
        ,'FA+DSG, DdeI+DpnII, 1kb'\
        ,'FA+DSG, Mnase, 1kb'\
        ,'FA, Dpn, 1kb'\
        ]

g_dist = 1


row = 0
col = 0


for e_idx, exp in enumerate(exps):

    #if e_idx != 3:
    #    continue

    fig, axs = plt.subplots(1, 2)

    reso = reso_num[e_idx]

    # source edge file for all

    target = '../cooler_results/PH_results/'+labels[e_idx]+'_gdist'+str(g_dist)+'_ALL_'

    ## H2 pers
    #pers_file = target+'H2_pers_data.txt'
    #try:
    #    pers = np.loadtxt(pers_file, delimiter=',')
    #    plot_pers_hist2d(pers, axs[row, col], plt_labels[e_idx])
    #except:
    #    print('empty')

    # H1 pers
    pers_file = target+'H1_pers_data.txt'
    #try:
    pers = np.loadtxt(pers_file, delimiter=',')

    bands1 = plot_pers_hist2d(pers, axs[0], 'Not scaled')

    #except:
    #    print('empty')


    target = '../cooler_results/PH_results/'+labels[e_idx]+'_gdist'+str(g_dist)+'_ALL_scaled_'
    pers_file = target+'H1_pers_data.txt'
    #try:
    pers = np.loadtxt(pers_file, delimiter=',')

    bands2 = plot_pers_hist2d(pers, axs[1], 'Scaled')

    plt.tight_layout()

    plt.show()
    plt.cla()
    plt.clf()
    plt.close()


    for idx, band1 in enumerate(bands1):

        band2 = bands2[idx]

        print(band1.shape, band2.shape)

        #print(hf2.KL_symm(band1, band2))

        #input('w')


    exit()

    #except:
    #    print('empty')



    #if col == 2:
    #    row += 1
    #    col = 0







        

        




