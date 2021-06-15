import pydory as dory

dirr = '../Datasets/'

compute_cycles = 0
reduce_cyc_lengths = 0
suppress_output = 0

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
dory.compute_PH(source, thresh, filetype, threads, target, dim, compute_cycles, reduce_cyc_lengths, suppress_output)


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
dory.compute_PH(source, thresh, filetype, threads, target, dim, compute_cycles, reduce_cyc_lengths, suppress_output)



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
dory.compute_PH(source, thresh, filetype, threads, target, dim, compute_cycles, reduce_cyc_lengths, suppress_output)



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
dory.compute_PH(source, thresh, filetype, threads, target, dim, compute_cycles, reduce_cyc_lengths, suppress_output)

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
dory.compute_PH(source, thresh, filetype, threads, target, dim, compute_cycles, reduce_cyc_lengths, suppress_output)

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
dory.compute_PH(source, thresh, filetype, threads, target, dim, compute_cycles, reduce_cyc_lengths, suppress_output)

