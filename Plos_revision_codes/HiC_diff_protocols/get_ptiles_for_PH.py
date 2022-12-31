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



# 5kb
# PH_thresh is 95%tile of 25kb
# birth thresh of 95%tile of 5kb, 10kb



# Resolution = 5kb
reso = 5000
labels = [\
         'fa_dsg_ddel_5kb'\
        ,'fa_dsg_dpn_5kb'\
        ,'fa_dsg_ddel_dpn_5kb'\
        ,'fa_dsg_mnase_5kb'\
        ,'fa_dpn_5kb'\
        ]

PH_pars_info = dict()

g_thresh = 10

for idx, label in enumerate(labels):

    print(idx, label)

    edges_file = '../cooler_results/edge_files/'+labels[idx]+'_bal_edges_upto_bindist'+str(g_thresh)+'.csv'
    edges= np.loadtxt(edges_file, delimiter=',')

    # remove NaN from balanced
    nan_idxs = np.isnan(edges[:,2])
    not_nan_idxs = np.argwhere(nan_idxs==False).flatten()
    edges = edges[not_nan_idxs]

    bin_dist = edges[:, 1] - edges[:, 0]


    # g_dist = 1 (5kb)
    this_idxs = np.argwhere(bin_dist == 1).flatten()
    this_edges = edges[this_idxs, 2]
    birth1 = np.percentile(this_edges, 95)

    # g_dist = 2 (10kb)
    this_idxs = np.argwhere(bin_dist == 2).flatten()
    this_edges = edges[this_idxs, 2]
    birth2 = np.percentile(this_edges, 95)

    # g_dist = 3 (15kb)
    this_idxs = np.argwhere(bin_dist == 3).flatten()
    this_edges = edges[this_idxs, 2]
    birth3 = np.percentile(this_edges, 95)

    # g_dist = 5 (25kb)
    this_idxs = np.argwhere(bin_dist == 5).flatten()
    this_edges = edges[this_idxs, 2]
    PH_thresh = np.percentile(this_edges, 95)

    PH_pars_info[label] = {\
                            'PH_thresh':PH_thresh\
                            , 'births':[birth1, birth2, birth3]\
                            }

pickle.dump(PH_pars_info, open('../cooler_results/edge_files/PH_pars_info_5kb.p', 'wb'))








