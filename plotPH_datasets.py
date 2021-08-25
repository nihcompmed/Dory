import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import math
mpl.rcParams['text.usetex'] = True


def plotPD(data1, label1, data2, label2, fname):

        
        max1 = np.amax(data1)
        max2 = np.amax(data2)
        maxx = 1.1*max(max1, max2)
        
        data1[data1[:,1] == -1, 1] = maxx
        data2[data2[:,1] == -1, 1] = maxx
        
        plt.scatter(data1[:,0], data1[:,1], marker = '+', alpha=0.85, color='tab:blue', label=label1)
        plt.scatter(data2[:,0], data2[:,1], marker = 'o' , alpha=0.6, label = label2\
                                , facecolors='none', edgecolors='tab:red', linewidths=2)
        ax = plt.gca()
        #ax.set_yticks([maxx])
        #ax.set_yticklabels([r'$\infty$'], fontsize=20)
        plt.hlines(maxx, xmin=0, xmax=maxx, ls='--', alpha=0.3, color='black')
        plt.plot([0, maxx], [0, maxx], ls='--', alpha=0.5, color='black')
        
        plt.xlabel('birth', fontsize=24)
        plt.ylabel('death', fontsize=24)
        
        plt.tight_layout()
        plt.legend()
        
        plt.legend(prop={'size': 20})

        
        #plt.show()
        plt.savefig('figures/'+fname+'.png', format='png')
        plt.savefig('figures/'+fname+'.pdf', format='pdf')
        
        plt.cla()
        plt.clf()



dirr = 'Datasets/ripser_PD_results/'

datasets = ['Dragon', 'Fract', 'o3', 'torus4']

label1 = 'Ripser'
label2 = 'Dory'

data1 = np.loadtxt(dirr+'ripsdragon.txt', delimiter=',')
data2 = np.loadtxt('Datasets/Dragon/DoryH1_pers_data.txt', delimiter=',')
plotPD(data1, label1, data2, label2, fname='ripsdragonH1')


data1 = np.loadtxt(dirr+'ripsfract.txt', delimiter=',')
data2 = np.loadtxt('Datasets/Fract/DoryH1_pers_data.txt', delimiter=',')
plotPD(data1, label1, data2, label2, fname='ripsfractH1')


data1 = np.loadtxt(dirr+'ripsfractH2.txt', delimiter=',')
data2 = np.loadtxt('Datasets/Fract/DoryH2_pers_data.txt', delimiter=',')
plotPD(data1, label1, data2, label2, fname='ripsfractH2')


data1 = np.loadtxt(dirr+'ripso3.txt', delimiter=',')
data2 = np.loadtxt('Datasets/o3/DoryH1_pers_data.txt', delimiter=',')
plotPD(data1, label1, data2, label2, fname='ripso3H1')


data1 = np.loadtxt(dirr+'ripso3H2.txt', delimiter=',')
data2 = np.loadtxt('Datasets/o3/DoryH2_pers_data.txt', delimiter=',')
plotPD(data1, label1, data2, label2, fname='ripso3H2')


data1 = np.loadtxt(dirr+'ripstorus4.txt', delimiter=',')
data2 = np.loadtxt('Datasets/torus4/DoryH1_pers_data.txt', delimiter=',')
plotPD(data1, label1, data2, label2, fname='ripstorusH1')

data1 = np.loadtxt(dirr+'ripstorus4H2.txt', delimiter=',')
data2 = np.loadtxt('Datasets/torus4/DoryH2_pers_data.txt', delimiter=',')
plotPD(data1, label1, data2, label2, fname='ripstorusH2')

label1 = 'Gudhi'
label2 = 'Dory'

data1 = np.loadtxt('Datasets/o3/Gudhi_H1.csv', delimiter=',')
data2 = np.loadtxt('Datasets/o3/DoryH1_pers_data.txt', delimiter=',')
plotPD(data1, label1, data2, label2, fname='Gudhio3H1')


data1 = np.loadtxt('Datasets/o3/Gudhi_H2.csv', delimiter=',')
data2 = np.loadtxt('Datasets/o3/DoryH2_pers_data.txt', delimiter=',')
plotPD(data1, label1, data2, label2, fname='Gudhio3H2')


data1 = np.loadtxt('Datasets/torus4/Gudhi_H1.csv', delimiter=',')
data2 = np.loadtxt('Datasets/torus4/DoryH1_pers_data.txt', delimiter=',')
plotPD(data1, label1, data2, label2, fname='GudhitorusH1')

data1 = np.loadtxt('Datasets/torus4/Gudhi_H2.csv', delimiter=',')
data2 = np.loadtxt('Datasets/torus4/DoryH2_pers_data.txt', delimiter=',')
plotPD(data1, label1, data2, label2, fname='GudhitorusH2')


label1 = 'Eirene'
label2 = 'Dory'

data1 = np.loadtxt('Datasets/Dragon/Eirene_H1.csv', delimiter=',')
data2 = np.loadtxt('Datasets/Dragon/DoryH1_pers_data.txt', delimiter=',')
plotPD(data1, label1, data2, label2, fname='EirenedragonH1')

data1 = np.loadtxt('Datasets/Fract/Eirene_H1.csv', delimiter=',')
data2 = np.loadtxt('Datasets/Fract/DoryH1_pers_data.txt', delimiter=',')
plotPD(data1, label1, data2, label2, fname='EirenefractH1')

data1 = np.loadtxt('Datasets/Fract/Eirene_H2.csv', delimiter=',')
data2 = np.loadtxt('Datasets/Fract/DoryH2_pers_data.txt', delimiter=',')
plotPD(data1, label1, data2, label2, fname='EirenefractH2')

data1 = np.loadtxt('Datasets/o3/Eirene_H1.csv', delimiter=',')
data2 = np.loadtxt('Datasets/o3/DoryH1_pers_data.txt', delimiter=',')
plotPD(data1, label1, data2, label2, fname='Eireneo3H1')


data1 = np.loadtxt('Datasets/o3/Eirene_H2.csv', delimiter=',')
data2 = np.loadtxt('Datasets/o3/DoryH2_pers_data.txt', delimiter=',')
plotPD(data1, label1, data2, label2, fname='Eireneo3H2')


data1 = np.loadtxt('Datasets/torus4/Eirene_H1.csv', delimiter=',')
data2 = np.loadtxt('Datasets/torus4/DoryH1_pers_data.txt', delimiter=',')
plotPD(data1, label1, data2, label2, fname='EirenetorusH1')



