import numpy as np
import timeit
import time
from ripser import Rips, ripser

dirr = 'Datasets/'

#################################
## dragon
#################################
start_time = time.time()

dataset = 'dragon'
d1 = dirr+'dragon2000_locs.csv'
data = np.loadtxt(d1, delimiter=',')

print('Processing', dataset)
dgms = ripser(data, maxdim=1, thresh=1000, distance_matrix=False)['dgms']

np.savetxt(dirr+'/'+dataset+'_Ripser_H1.csv', dgms[1], delimiter=',')

print('time taken', time.time() - start_time)

#################################

#################################
## fract
#################################
start = time.time()
dataset = 'fract'
d1 = dirr+'fractal_r_distmat.csv'
data = np.loadtxt(d1, delimiter=',')

print('Processing', dataset)
dgms = ripser(data, maxdim=2, thresh=1000, distance_matrix=True)['dgms']

np.savetxt(dirr+'/'+dataset+'_Ripser_H1.csv', dgms[1], delimiter=',')
np.savetxt(dirr+'/'+dataset+'_Ripser_H2.csv', dgms[2], delimiter=',')

print('time taken', time.time() - start_time)


#################################
## o3
#################################

start_time = time.time()
dataset = 'o3'
d1 = dirr+'o3_8192_locs.csv'
data = np.loadtxt(d1, delimiter=',')

print('Processing', dataset)
dgms = ripser(data, maxdim=2, thresh=1, distance_matrix=False)['dgms']


np.savetxt(dirr+'/'+dataset+'_Ripser_H1.csv', dgms[1], delimiter=',')
np.savetxt(dirr+'/'+dataset+'_Ripser_H2.csv', dgms[2], delimiter=',')

print('time taken', time.time() - start_time)


#################################
## torus4 (H1)
#################################
start = time.time()
dataset = 'torus4'
d1 = dirr+'torus4_locs.csv'
data = np.loadtxt(d1, delimiter=',')

print('Processing', dataset)

dgms = ripser(data, maxdim=1, thresh=0.15, distance_matrix=False)['dgms']


np.savetxt(dirr+'/'+dataset+'_Ripser_H1.csv', dgms[1], delimiter=',')

print('time taken', time.time() - start_time)

#################################
## torus4 (H2)
#################################
start = time.time()
dataset = 'torus4'
d1 = dirr+'torus4_locs.csv'
data = np.loadtxt(d1, delimiter=',')

print('Processing', dataset)

dgms = ripser(data, maxdim=2, thresh=0.15, distance_matrix=False)['dgms']


np.savetxt(dirr+'/'+dataset+'_Ripser_H1.csv', dgms[1], delimiter=',')
np.savetxt(dirr+'/'+dataset+'_Ripser_H2.csv', dgms[2], delimiter=',')

print('time taken', time.time() - start_time)



