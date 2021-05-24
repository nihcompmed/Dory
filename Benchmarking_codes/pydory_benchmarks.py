import pydory as dory

dirr = '../Datasets/'

#################################
## dragon
#################################

dataset = 'dragon'
source = dirr+'dragon2000_locs.csv'
target = dirr + 'Dory_' + dataset + '_'

print('#########################')
print('\nProcessing', dataset)
print('\n#########################')

thresh = 1000
threads = 4
filetype = 1
dim = 1
compute_cycles = 0
reduce_cyc_lengths = 0
dory.compute_PH(source, thresh, filetype, threads, target, dim, compute_cycles, reduce_cyc_lengths)


#################################

#################################
## fract
#################################
dataset = 'fract'
source = dirr+'fractal_r_distmat.csv'
target = dirr + 'Dory_' + dataset + '_'

print('#########################')
print('\nProcessing', dataset)
print('\n#########################')
thresh = 1000
threads = 4
filetype = 0
dim = 2
dory.compute_PH(source, thresh, filetype, threads, target, dim, compute_cycles, reduce_cyc_lengths)



#################################
## o3
#################################

dataset = 'o3'
source = dirr+'o3_8192_locs.csv'
target = dirr + 'Dory_' + dataset + '_'
print('#########################')
print('\nProcessing', dataset)
print('\n#########################')

thresh = 1
threads = 4
filetype = 1
dim = 2
dory.compute_PH(source, thresh, filetype, threads, target, dim, compute_cycles, reduce_cyc_lengths)



#################################
## torus4 (H2)
#################################

dataset = 'torus4'
source = dirr+'torus4_locs.csv'
target = dirr + 'Dory_' + dataset + '_'

print('#########################')
print('\nProcessing', dataset)
print('\n#########################')

thresh = 0.15
threads = 4
filetype = 1
dim = 2
dory.compute_PH(source, thresh, filetype, threads, target, dim, compute_cycles, reduce_cyc_lengths)

#################################
## HiC control
#################################

dataset = 'HiC_control'
source = dirr+'HiC/control_400.csv'
target = dirr + 'Dory_' + dataset + '_'

print('#########################')
print('\nProcessing', dataset)
print('\n#########################')

thresh = 400
threads = 4
filetype = 2
dim = 2
dory.compute_PH(source, thresh, filetype, threads, target, dim, compute_cycles, reduce_cyc_lengths)

#################################
## HiC auxin
#################################

dataset = 'HiC_auxin'
source = dirr+'HiC/auxin_400.csv'
target = dirr + 'Dory_' + dataset + '_'

print('#########################')
print('\nProcessing', dataset)
print('\n#########################')

thresh = 400
threads = 4
filetype = 2
dim = 2
dory.compute_PH(source, thresh, filetype, threads, target, dim, compute_cycles, reduce_cyc_lengths)

