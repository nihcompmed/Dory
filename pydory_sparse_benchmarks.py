import time
import pydory as dory

dirr = 'Datasets/'

#################################
## dragon
#################################

dataset = 'dragon'
source = dirr+'dragon2000_locs.csv'
target = dirr + 'DoryS_' + dataset + '_'

print('#########################')
print('\nProcessing', dataset)
print('\n#########################')

thresh = 1000
threads = 4
filetype = 1
dim = 1
dory.compute_PH(source, thresh, filetype, threads, target, dim)


#################################

#################################
## fract
#################################
dataset = 'fract'
source = dirr+'fractal_r_distmat.csv'
target = dirr + 'DoryS_' + dataset + '_'

print('#########################')
print('\nProcessing', dataset)
print('\n#########################')
thresh = 1000
threads = 4
filetype = 0
dim = 2
dory.compute_PH(source, thresh, filetype, threads, target, dim)



#################################
## o3
#################################

dataset = 'o3'
source = dirr+'o3_8192_locs.csv'
target = dirr + 'DoryS_' + dataset + '_'
print('#########################')
print('\nProcessing', dataset)
print('\n#########################')

thresh = 1
threads = 4
filetype = 1
dim = 2
dory.compute_PH(source, thresh, filetype, threads, target, dim)



#################################
## torus4 (H2)
#################################

dataset = 'torus4'
source = dirr+'torus4_locs.csv'
target = dirr + 'DoryS_' + dataset + '_'

print('#########################')
print('\nProcessing', dataset)
print('\n#########################')

thresh = 0.15
threads = 4
filetype = 1
dim = 2
dory.compute_PH(source, thresh, filetype, threads, target, dim)




