import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from scipy import interpolate
mpl.rcParams['text.usetex'] = True

def get_info(pers_pairs):



    b_c = pers_pairs[:,0]
    
    d_c = pers_pairs[:,1]
    
    b_c_uni, b_counts = np.unique(b_c, return_counts = True)
    b_c_cumcounts = np.cumsum(b_counts)
    
    d_c = d_c[d_c > -1]
    
    d_c_uni, d_counts = np.unique(d_c, return_counts = True)
    d_c_cumcounts = np.cumsum(d_counts)
    d_counts = -d_counts
    
    
    all_c_counts = np.array(list(b_counts) + list(d_counts))
    
    all_c_idxs = np.array(list(b_c_uni) + list(d_c_uni))
    
    sort_idxs = np.argsort(all_c_idxs)
    
    all_c_idxs = all_c_idxs[sort_idxs]
    all_c_counts = all_c_counts[sort_idxs]
    
    c_sum = np.cumsum(all_c_counts)

    return all_c_idxs, c_sum





f_c = 'Datasets/HiC/control400H1_pers_data.txt'
pers_pairs_c = np.loadtxt(f_c, delimiter=',')

#EDIT
f_a = 'Datasets/HiC/auxin400H1_pers_data.txt'
pers_pairs_a = np.loadtxt(f_a, delimiter=',')

all1, c1 = get_info(pers_pairs_c)

all2, c2 = get_info(pers_pairs_a)

min1 = np.amin(all1)
min2 = np.amin(all2)
xmin = max(min1, min2)

fcontrol = interpolate.interp1d(all1, c1)

fauxin = interpolate.interp1d(all2, c2)

xnew = np.arange(xmin, 400, 1)

controlnew = fcontrol(xnew)
auxinnew = fauxin(xnew)

diff = (auxinnew - controlnew)/controlnew*100

plt.plot(xnew, diff, lw=2, alpha=0.85, label = r'\#loops', color = 'tab:red')



#EDIT
f_c = 'Datasets/HiC/control400H2_pers_data.txt'
pers_pairs_c = np.loadtxt(f_c, delimiter=',')

#EDIT
f_a = 'Datasets/HiC/auxin400H2_pers_data.txt'
pers_pairs_a = np.loadtxt(f_a, delimiter=',')

all1, c1 = get_info(pers_pairs_c)

all2, c2 = get_info(pers_pairs_a)

min1 = np.amin(all1)
min2 = np.amin(all2)
xmin = max(min1, min2)

fcontrol = interpolate.interp1d(all1, c1)

fauxin = interpolate.interp1d(all2, c2)

xnew = np.arange(xmin, 400, 1)

controlnew = fcontrol(xnew)
auxinnew = fauxin(xnew)

diff = (auxinnew - controlnew)/controlnew*100


#plt.ticklabel_format(axis="y", style="sci", scilimits=(0,0))
#plt.scatter(b_c_uni, b_c_cumcounts)
#plt.scatter(d_c_uni, d_c_cumcounts)

plt.rc('xtick',labelsize=20)
plt.rc('ytick',labelsize=20)

plt.xticks(fontsize=18)
plt.yticks(fontsize=18)

plt.plot(xnew, diff, lw=2, alpha=0.85, label = r'\#voids', color = 'tab:blue')



#plt.plot(all1, c1, lw=4, alpha = 0.85, label = 'control')
#plt.plot(all2, c2, lw=4, ls='--', alpha = 0.85, label = 'with auxin')

#EDIT

#plt.ylabel(r'$\mathbf{\beta_1}$ (\#holes)', fontsize=16)
#plt.ylabel(r'$\mathbf{\beta_2}$ (\#voids)', fontsize=16)
plt.ylabel(r'$\%$ change', fontsize=16)

plt.xlabel('Parameter', fontsize=16)
plt.tight_layout()
plt.legend(prop={'size':14})

#plt.show()
#EDIT
plt.savefig('../HiC_num_diff.pdf', format = 'pdf')
#plt.savefig('../HiC_diff_b1.pdf', format = 'pdf')
#plt.savefig('../HiC_diff_b2.pdf', format = 'pdf')


