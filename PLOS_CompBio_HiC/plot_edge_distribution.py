import numpy as np
import cooler
import h5py
import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn as sns
from numba import njit
import pickle
import math


fig, axs = plt.subplots(1, 2, sharey=True)

exp_name = '4DNFITH9KJF1.cool'
data = pickle.load(open(exp_name+"bin_dist_par.p", 'rb'))
bin_dist = data['bin_dist']
spatial = data['spatial_dist']
spatial_corrected = data['spatial_dist_corrected']
sns.boxplot(x=bin_dist, y=spatial_corrected, ax=axs[0])

exp_name = '4DNFI2X45Z5L.cool'
data = pickle.load(open(exp_name+"bin_dist_par.p", 'rb'))
bin_dist = data['bin_dist']
spatial = data['spatial_dist']
spatial_corrected = data['spatial_dist_corrected']
sns.boxplot(x=bin_dist, y=spatial_corrected, ax=axs[1])


axs[0].hlines(150, 0, 4, ls='--', color='red', lw=2, alpha=0.5)
axs[1].hlines(150, 0, 4, ls='--', color='red', lw=2, alpha=0.5)

axs[0].hlines(100, 0, 4, ls='-.', color='black', lw=2, alpha=0.5)
axs[1].hlines(100, 0, 4, ls='-.', color='black', lw=2, alpha=0.5)

axs[0].set_yscale("log", base=10)
axs[1].set_yscale("log", base=10)

axs[0].set_xlabel('bin distance', fontsize=14)
axs[1].set_xlabel('bin distance', fontsize=14)

axs[0].set_ylabel('pairwise distance estimate', fontsize=14)

axs[0].set_title('Control')
axs[1].set_title('Auxin-treated')


plt.tight_layout()

plt.savefig('figures/edge_distribution.pdf', dpi=600)

plt.cla()
plt.clf()

