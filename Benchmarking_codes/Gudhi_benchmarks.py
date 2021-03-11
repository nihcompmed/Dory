import numpy as np
import gudhi as gd  

import time

dirr = '../Datasets/'

#######################################
### dragon
#######################################
#start = time.time()
#dataset = 'dragon'
#d1 = dirr+'dragon2000_locs.csv'
#
#data = np.loadtxt(d1, delimiter=',')
#
#print('Processing ', dataset)
#
#skeleton = gd.RipsComplex(
#    points = data, 
#    max_edge_length = 1000
#)
#
#Rips_simplex_tree = skeleton.create_simplex_tree(max_dimension = 2)
#
#BarCodes_Rips0 = Rips_simplex_tree.persistence()
#
#vals = Rips_simplex_tree.persistence_intervals_in_dimension(1)
#np.savetxt(dirr+'/dragon_Gudhi_H1.csv', vals, delimiter=',')
#
#print('Time taken ', time.time() - start)
#
#######################################
#
#
#######################################
### fract
#######################################
#start = time.time()
#dataset = 'fract'
#d1 = dirr+'fractal_r_distmat.csv'
#
#data = np.loadtxt(d1, delimiter=',')
#
#print('Processing ', dataset)
#
#skeleton = gd.RipsComplex(
#    distance_matrix = data, 
#    max_edge_length = 1000
#)
#
#Rips_simplex_tree = skeleton.create_simplex_tree(max_dimension = 3)
#
#BarCodes_Rips0 = Rips_simplex_tree.persistence()
#
#vals = Rips_simplex_tree.persistence_intervals_in_dimension(1)
#np.savetxt(dirr+'/fract_Gudhi_H1.csv', vals, delimiter=',')
#
#vals = Rips_simplex_tree.persistence_intervals_in_dimension(2)
#np.savetxt(dirr+'/fract_Gudhi_H2.csv', vals, delimiter=',')
#
#print('Time taken ', time.time() - start)
#
#######################################

######################################
## o3
######################################
start = time.time()
dataset = 'o3'
d1 = dirr+'o3_8192_locs.csv'

print('Processing ', dataset)

data = np.loadtxt(d1, delimiter=',')

skeleton = gd.RipsComplex(
    points = data, 
    max_edge_length = 1
)


Rips_simplex_tree = skeleton.create_simplex_tree(max_dimension = 3)

BarCodes_Rips0 = Rips_simplex_tree.persistence()

vals = Rips_simplex_tree.persistence_intervals_in_dimension(1)
np.savetxt(dirr+'/o3_Gudhi_H1.csv', vals, delimiter=',')

vals = Rips_simplex_tree.persistence_intervals_in_dimension(2)
np.savetxt(dirr+'/o3_Gudhi_H2.csv', vals, delimiter=',')


print('Time taken ',time.time()-start)


######################################
## torus4
######################################
start = time.time()
dataset = 'torus4'
d1 = dirr+'torus4_locs.csv'

print('Processing ', dataset)

data = np.loadtxt(d1, delimiter=',')

skeleton = gd.RipsComplex(
    points = data, 
    max_edge_length = 0.15
)


Rips_simplex_tree = skeleton.create_simplex_tree(max_dimension = 3)

BarCodes_Rips0 = Rips_simplex_tree.persistence()

vals = Rips_simplex_tree.persistence_intervals_in_dimension(1)
np.savetxt(dirr+'/torus4_Gudhi_H1.csv', vals, delimiter=',')

vals = Rips_simplex_tree.persistence_intervals_in_dimension(2)
np.savetxt(dirr+'/torus4_Gudhi_H2.csv', vals, delimiter=',')

finish = time.time()

print('Time taken ',time.time()-start)


