import pickle
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import seaborn as sns


info_dict = pickle.load(open('PD_distances_pairwise.p', 'rb'))


all_L0 = info_dict['vals']['L0']
all_bn = info_dict['vals']['bn']
all_ws = info_dict['vals']['ws']
all_sw = info_dict['vals']['sw']

pickle.dump(info_dict, open('PD_distances_pairwise.p', 'wb'))

fig, axs = plt.subplots(1, 3)

axs[0].hist2d(all_L0, all_bn, bins=50, norm=mpl.colors.LogNorm(), cmap='gnuplot')
axs[0].set_xlabel('L0 norm')
axs[0].set_ylabel('bottleneck')


axs[1].hist2d(all_L0, all_ws, bins=50, norm=mpl.colors.LogNorm(), cmap='gnuplot')
axs[1].set_xlabel('L0 norm')
axs[1].set_ylabel('wasserstein')

axs[2].hist2d(all_L0, all_sw, bins=50, norm=mpl.colors.LogNorm(), cmap='gnuplot')
axs[2].set_xlabel('L0 norm')
axs[2].set_ylabel('sliced wasserstein')

plt.tight_layout()

plt.savefig('proteins_compare_PD_distances.pdf', dpi=600)


#plt.show()

plt.cla()
plt.clf()
plt.close()

all_L0 = np.array(info_dict['times']['L0'])
all_bn = np.array(info_dict['times']['bn'])
all_ws = np.array(info_dict['times']['ws'])
all_sw = np.array(info_dict['times']['sw'])

bn_ratios = np.log10(all_bn/all_L0)
ws_ratios = np.log10(all_ws/all_L0)
sw_ratios = np.log10(all_sw/all_L0)

ax = plt.gca()

sns.histplot(bn_ratios, kde=True, alpha=0.5, label='Bottleneck',ax=ax)
sns.histplot(ws_ratios, kde=True, alpha=0.5, label='Wasserstein',ax=ax)
sns.histplot(sw_ratios, kde=True, alpha=0.5, label='Sliced Wasserstein',ax=ax)

plt.ylabel('Count', fontsize=16)
plt.xlabel(r'log$_{10}$(time/(L$_0$ time))', fontsize=16)

plt.legend()

#plt.cla()
#plt.clg()
#plt.close()

#plt.show()
#exit()

#pickle.dump(info_dict, open('PD_distances_pairwise.p', 'wb'))

#fig, axs = plt.subplots(1, 3)
#
#axs[0].hist2d(all_L0, all_bn, bins=50, norm=mpl.colors.LogNorm(), cmap='gnuplot')
#axs[0].set_xlabel('L0 norm')
#axs[0].set_ylabel('bottleneck')
#
#
#axs[1].hist2d(all_L0, all_ws, bins=50, norm=mpl.colors.LogNorm(), cmap='gnuplot')
#axs[1].set_xlabel('L0 norm')
#axs[1].set_ylabel('wasserstein')
#
#axs[2].hist2d(all_L0, all_sw, bins=50, norm=mpl.colors.LogNorm(), cmap='gnuplot')
#axs[2].set_xlabel('L0 norm')
#axs[2].set_ylabel('sliced wasserstein')
#
#plt.tight_layout()


plt.savefig('proteins_compare_PD_distances_computation_times.pdf', dpi=600)
#plt.show()

#plt.savefig('proteins_compare_PD_distances.pdf', dpi=600)


#plt.show()

plt.cla()
plt.clf()
plt.close()


