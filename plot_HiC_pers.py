import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
mpl.rcParams['text.usetex'] = True

def plot_pers(ff, fname):

    print('plotting pers pairs', ff)

    pers_pairs = np.loadtxt(ff, delimiter=',')

    plt.hist2d(pers_pairs[:,0], pers_pairs[:,1], bins=100,norm=mpl.colors.LogNorm(),cmap=mpl.cm.jet)
    
    cbar = plt.colorbar()
    
    cbar.set_label(r'\#pairs')
    
    plt.rc('xtick',labelsize=14)
    plt.rc('ytick',labelsize=14)
    
    plt.xlabel('birth', fontsize=16)
    plt.ylabel('death', fontsize=16)

    plt.tight_layout()
    
    plt.savefig('figures/'+fname+'.png', format = 'png')

    plt.cla()
    plt.clf()


ff = 'Datasets/HiC/Dory_control_H1_pers_data.txt'
plot_pers(ff, 'HiC_control_H1')
pers_pairs = np.loadtxt(ff, delimiter=',')

ff = 'Datasets/HiC/Dory_control_H2_pers_data.txt'
plot_pers(ff, 'HiC_control_H2')
pers_pairs = np.loadtxt(ff, delimiter=',')


ff = 'Datasets/HiC/Dory_auxin_H1_pers_data.txt'
plot_pers(ff, 'HiC_auxin_H1')
pers_pairs = np.loadtxt(ff, delimiter=',')

ff = 'Datasets/HiC/Dory_auxin_H2_pers_data.txt'
plot_pers(ff, 'HiC_auxin_H2')
pers_pairs = np.loadtxt(ff, delimiter=',')
