import matplotlib.pyplot as plt
import numpy as np
import pandas
import math
import h5py
from numba import njit, prange
from numba import njit
import pickle
import cooler 
import seaborn as sns
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

# This creates only plots for genomic dist vs. distance estimate for ICE correction

exps_dir = '../cooler_files/'

exps = [\
         'fa_dsg_ddel_5kb.cool'\
        ,'fa_dsg_dpn_5kb.cool'\
        ,'fa_dsg_ddel_dpn_5kb.cool'\
        ,'fa_dsg_mnase_5kb.cool'\
        ,'fa_dpn_5kb.cool'\
        ,'fa_dsg_ddel_1kb.cool'\
        ,'fa_dsg_dpn_1kb.cool'\
        ,'fa_dsg_ddel_dpn_1kb.cool'\
        ,'fa_dsg_mnase_1kb.cool'\
        ,'fa_dpn_1kb.cool'\
       ]

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


plt_labels = [\
         'fa+dsg, ddel, 5kb'\
        ,'fa+dsg, dpn, 5kb'\
        ,'fa+dsg, ddel+dpn, 5kb'\
        ,'fa+dsg, mnase, 5kb'\
        ,'fa, dpn, 5kb'\
        ,'fa+dsg, ddel, 1kb'\
        ,'fa+dsg, dpn, 1kb'\
        ,'fa+dsg, ddel+dpn, 1kb'\
        ,'fa+dsg, mnase, 1kb'\
        ,'fa, dpn, 1kb'\
        ]

def get_ptiles(raw_edges, bal_edges, g_thresh):

    bin_dist = raw_edges[:, 1] - raw_edges[:, 0]

    xx = []
    raw_yy = []
    bal_yy = []

    percentiles = list(range(85,100))

    fig, axs = plt.subplots(1, 2)

    for g_dist in list(range(1, g_thresh+1)):

        this_idxs = np.argwhere(bin_dist == g_dist).flatten()

        this_raw_vals = raw_edges[this_idxs, 2]
        this_bal_vals = bal_edges[this_idxs, 2]

        print(this_raw_vals, this_bal_vals)

        # remove NaN from balanced
        nan_idxs = np.isnan(this_bal_vals)
        not_nan_idxs = np.argwhere(nan_idxs==False).flatten()
        this_bal_vals = this_bal_vals[not_nan_idxs]

        raw_ptiles = np.percentile(this_raw_vals, percentiles)
        bal_ptiles = np.percentile(this_bal_vals, percentiles)

        axs[0].plot(percentiles, raw_ptiles, ls = '-')
        axs[1].plot(percentiles, bal_ptiles, ls = '--')

    axs[0].plot()

    plt.show()

    exit()

def get_box_plots(edges, ax, reso):


    # remove NaN from balanced
    nan_idxs = np.isnan(edges[:,2])
    not_nan_idxs = np.argwhere(nan_idxs==False).flatten()
    edges = edges[not_nan_idxs]

    bin_dist = edges[:, 1] - edges[:, 0]

    sns.boxplot(x=bin_dist*reso/1000, y=edges[:,2], ax=ax)

    percentiles = [5, 95]

    for g_dist in list(range(1, g_thresh+1)):

        this_idxs = np.argwhere(bin_dist == g_dist).flatten()

        this_bal_vals = edges[this_idxs, 2]

        # remove NaN from balanced
        nan_idxs = np.isnan(this_bal_vals)
        not_nan_idxs = np.argwhere(nan_idxs==False).flatten()
        this_bal_vals = this_bal_vals[not_nan_idxs]

        bal_ptiles = np.percentile(this_bal_vals, percentiles)

        ax.scatter(g_dist-1, bal_ptiles[0]\
                    , marker='x', color='red', s=100, zorder=1000)
        ax.scatter(g_dist-1, bal_ptiles[1]\
                    , marker='o', color='red', s=50, zorder=1000)



    legend_elements = [\
                   Line2D([0], [0], marker='x', color='r', label='5%-tile',
                          markerfacecolor='r', markersize=15),
                   Line2D([0], [0], marker='o', color='r', label='95%-tile',
                          markerfacecolor='r', markersize=10),
                                 ]
    


    ax.legend(handles=legend_elements)


## PLOT EXAMPLE WITH BOX PLOTS
#idx = 0
g_thresh = 10
#
## Resolution 5kb
#raw_edges_file = '../cooler_results/edge_files/'+labels[idx]+'_raw_edges_upto_bindist'+str(g_thresh)+'.csv'
#bal_edges_file = '../cooler_results/edge_files/'+labels[idx]+'_bal_edges_upto_bindist'+str(g_thresh)+'.csv'
#
##raw_edges= np.loadtxt(raw_edges_file, delimiter=',')
#bal_edges= np.loadtxt(bal_edges_file, delimiter=',')
#
##get_ptiles(raw_edges, bal_edges, g_thresh)
#
#ax = plt.gca()
#
#reso = 5000
#
#get_box_plots(bal_edges, ax, reso)
#
#plt.yscale('log', base=2)
#
#plt.xlabel('bin distance in kb')
#plt.ylabel('spatial estimate = 1/(balanced count)')
#
#plt.title(plt_labels[idx])
#
#plt.savefig('../cooler_results/figures/example_boxplot_estimate_vs_genomic_distance.pdf', dpi=600)


plt.cla()
plt.clf()
plt.close()

# Plot 5 and 95%tiles for all experiments at different resolutions

def plot_ptiles(edges, ax, color, reso, scale, shift):

    # remove NaN from balanced
    nan_idxs = np.isnan(edges[:,2])
    not_nan_idxs = np.argwhere(nan_idxs==False).flatten()
    edges = edges[not_nan_idxs]

    bin_dist = edges[:, 1] - edges[:, 0]

    percentiles = [5, 95]

    xx = []
    yy_5 = []
    yy_95 = []

    for g_dist in list(range(1, g_thresh+1)):

        this_idxs = np.argwhere(bin_dist == g_dist).flatten()

        this_bal_vals = edges[this_idxs, 2]

        # remove NaN from balanced
        nan_idxs = np.isnan(this_bal_vals)
        not_nan_idxs = np.argwhere(nan_idxs==False).flatten()
        this_bal_vals = this_bal_vals[not_nan_idxs]

        bal_ptiles = np.percentile(this_bal_vals, percentiles)

        xx.append(g_dist*reso/1000)
        yy_5.append(shift+bal_ptiles[0]**scale)
        yy_95.append(shift+bal_ptiles[1]**scale)

    ax.plot(xx, yy_5, ls='--', color=color)
    ax.plot(xx, yy_95, color=color)

    #legend_elements = [\
    #               Line2D([0], [0], marker='x', color='r', label='5%-tile',
    #                      markerfacecolor='r', markersize=15),
    #               Line2D([0], [0], marker='o', color='r', label='95%-tile',
    #                      markerfacecolor='r', markersize=10),
    #                             ]
    #


    #ax.legend(handles=legend_elements)


# Resolution = 5kb
reso = 5000
labels = [\
         'fa_dsg_ddel_5kb'\
        ,'fa_dsg_dpn_5kb'\
        ,'fa_dsg_ddel_dpn_5kb'\
        ,'fa_dsg_mnase_5kb'\
        ,'fa_dpn_5kb'\
        ]

colors = [\
        '#FFC20A'\
        ,'#0C7BDC'\
        ,'#994F00'\
        ,'#DC3220'\
        ,'black'\
        ]

scaling_factor = [
        1\
        ,1\
        ,1\
        ,0.915\
        ,1\
        ]

shift_factor = [
        0\
        ,0\
        ,0\
        ,7\
        ,0\
        ]

ignore_idxs = [4]

for idx, label in enumerate(labels):

    if idx in ignore_idxs:
        continue

    print(idx, label)

    bal_edges_file = '../cooler_results/edge_files/'+labels[idx]+'_bal_edges_upto_bindist'+str(g_thresh)+'.csv'
    bal_edges= np.loadtxt(bal_edges_file, delimiter=',')

    ax = plt.gca()

    plot_ptiles(bal_edges, ax, colors[idx], reso, scaling_factor[idx], shift_factor[idx])


# Legend
ax = plt.gca()

legend_elements = [\
                   Patch(facecolor=colors[0], label=plt_labels[0]),
                   Patch(facecolor=colors[1], label=plt_labels[1]),
                   Patch(facecolor=colors[2], label=plt_labels[2]),
                   Patch(facecolor=colors[3], label=plt_labels[3]),
                   Patch(facecolor=colors[4], label=plt_labels[4]),
                   Line2D([0], [0], color='k', lw=2, ls='--', label='5%-tile'),
                   Line2D([0], [0], color='k', lw=2, label='95%-tile'),
                    ]

ax.legend(handles=legend_elements)
plt.xlabel('bin distance in kb')
plt.ylabel('spatial estimate = 1/(balanced counts)')

plt.show()

exit()

#plt.savefig('../')
plt.savefig('../cooler_results/figures/ptiles_estimate_vs_genomic_distance_5kb.pdf', dpi=600)

plt.cla()
plt.clf()
plt.close()


# Resolution = 5kb
reso = 1000
labels = [\
         'fa_dsg_ddel_1kb'\
        ,'fa_dsg_dpn_1kb'\
        ,'fa_dsg_ddel_dpn_1kb'\
        ,'fa_dsg_mnase_1kb'\
        ,'fa_dpn_1kb'\
        ]

plt_labels = [\
        'fa+dsg, ddel, 1kb'\
        ,'fa+dsg, dpn, 1kb'\
        ,'fa+dsg, ddel+dpn, 1kb'\
        ,'fa+dsg, mnase, 1kb'\
        ,'fa, dpn, 1kb'\
        ]

colors = [\
        '#FFC20A'\
        ,'#0C7BDC'\
        ,'#994F00'\
        ,'#DC3220'\
        ,'black'\
        ]

for idx, label in enumerate(labels):

    print(idx, label)

    bal_edges_file = '../cooler_results/edge_files/'+labels[idx]+'_bal_edges_upto_bindist'+str(g_thresh)+'.csv'
    bal_edges= np.loadtxt(bal_edges_file, delimiter=',')

    ax = plt.gca()

    plot_ptiles(bal_edges, ax, colors[idx], reso, scale)


# Legend
ax = plt.gca()

legend_elements = [\
                   Patch(facecolor=colors[0], label=plt_labels[0]),
                   Patch(facecolor=colors[1], label=plt_labels[1]),
                   Patch(facecolor=colors[2], label=plt_labels[2]),
                   Patch(facecolor=colors[3], label=plt_labels[3]),
                   Patch(facecolor=colors[4], label=plt_labels[4]),
                   Line2D([0], [0], color='k', lw=2, ls='--', label='5%-tile'),
                   Line2D([0], [0], color='k', lw=2, label='95%-tile'),
                    ]

ax.legend(handles=legend_elements)
plt.xlabel('bin distance in kb')
plt.ylabel('spatial estimate = 1/(balanced counts)')

#plt.savefig('../')
plt.savefig('../cooler_results/figures/ptiles_estimate_vs_genomic_distance_1kb.pdf', dpi=600)

plt.cla()
plt.clf()
plt.close()












