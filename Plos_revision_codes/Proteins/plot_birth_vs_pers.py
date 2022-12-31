import pickle
import persim
import itertools
import networkx as nx
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl




#ff = 'all_final_results_v2.p'
#ff = 'all_final_results.p'
ff = 'pers_dict_homologs.p'
#ff = 'shortlisted_sets.p'
#ff = 'sig_pdb.p'

thresh = 10
epsil = 3.5

count = 0

pers_data = pickle.load(open(ff, 'rb'))

all_births = []
all_pers = []

for pdb in pers_data:

    pers1 = pers_data[pdb.lower()]
    pers1[pers1[:,1]==-1, 1] = thresh+epsil

    diff = pers1[:,1] - pers1[:,0]

    all_births += list(pers1[:,0])
    all_pers += list(diff)

plt.hist2d(all_births, all_pers, bins=100, norm=mpl.colors.LogNorm(), cmap='gnuplot')

plt.xlabel('births', fontsize=14)
plt.ylabel('persistence', fontsize=14)
plt.axhline(y=3.5, ls='--', color='red', lw=2, alpha=0.8)
#plt.axhline(y=3.25, ls='--', color='red', lw=2, alpha=0.8)
#plt.axhline(y=3, ls='--', color='red', lw=2, alpha=0.8)
plt.axhline(y=2.75, ls='--', color='red', lw=2, alpha=0.8)
plt.axvline(x=10, ls='--', color='cyan', lw=2, alpha=0.8)
plt.colorbar()

plt.savefig('proteins_birth_vs_persistence.pdf', dpi=600)

#plt.show()


