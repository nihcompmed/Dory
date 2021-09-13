import pydory as dory
import sys
import time


main_dirr = 'Datasets/'


compute_cycles = 0
reduce_cyc_lengths = 0
suppress_output = 1
threads = 4
lower_thresh = 0

#################################
## dragon
#################################

dataset = 'Dragon'
dirr = main_dirr+dataset+'/'
source = dirr+'dragon2000_locs.csv'
target = dirr + 'Dory'

print('#########################')
print('Processing', dataset)
print('#########################')

thresh = 1000
filetype = 1
dim = 1



start = time.time()

dory.compute_PH(source, lower_thresh, thresh, filetype, threads, target, dim, compute_cycles, reduce_cyc_lengths, suppress_output)
print('Time taken for', dataset, ':', time.time() - start)


#################################

#################################
## fract
#################################
dataset = 'fract'
dirr = main_dirr+dataset+'/'
source = dirr+'fractal_r_distmat.csv'
target = dirr + 'Dory'

print('#########################')
print('Processing', dataset)
print('#########################')
thresh = 1000
filetype = 0
dim = 2

start = time.time()
dory.compute_PH(source, lower_thresh, thresh, filetype, threads, target, dim, compute_cycles, reduce_cyc_lengths, suppress_output)
print('Time taken for', dataset, ':', time.time() - start)



#################################
## o3
#################################

dataset = 'o3'
dirr = main_dirr+dataset+'/'
source = dirr+'o3_8192.csv'
target = dirr + 'Dory'

print('#########################')
print('Processing', dataset)
print('#########################')

thresh = 1
filetype = 1
dim = 2

start = time.time()
dory.compute_PH(source, lower_thresh, thresh, filetype, threads, target, dim, compute_cycles, reduce_cyc_lengths, suppress_output)
print('Time taken for', dataset, ':', time.time() - start)

#################################
## torus4 (H1)
#################################

dataset = 'torus4'
dirr = main_dirr+dataset+'/'
source = dirr+'torus4_locs.csv'
target = dirr + 'Dory'

print('#########################')
print('Processing', dataset)
print('#########################')

thresh = 0.15
filetype = 1
dim = 1

start = time.time()
dory.compute_PH(source, lower_thresh, thresh, filetype, threads, target, dim, compute_cycles, reduce_cyc_lengths, suppress_output)
print('Time taken for', dataset, ' (H1):', time.time() - start)


#################################
## torus4 (H2)
#################################

dataset = 'torus4'
dirr = main_dirr+dataset+'/'
source = dirr+'torus4_locs.csv'
target = dirr + 'Dory'

print('#########################')
print('Processing', dataset)
print('#########################')

thresh = 0.15
filetype = 1
dim = 2

start = time.time()
dory.compute_PH(source, lower_thresh, thresh, filetype, threads, target, dim, compute_cycles, reduce_cyc_lengths, suppress_output)
print('Time taken for', dataset, ' (H2):', time.time() - start)

#################################
## HiC control
#################################

dataset = 'HiC'
dirr = main_dirr+dataset+'/'
source = dirr+'control_500.csv'
target = dirr + 'Dory_control_'

print('#########################')
print('Processing', dataset)
print('#########################')

thresh = 400
filetype = 2
dim = 2

start = time.time()
dory.compute_PH(source, lower_thresh, thresh, filetype, threads, target, dim, compute_cycles, reduce_cyc_lengths, suppress_output)
print('Time taken for', dataset, ':', time.time() - start)

#################################
## HiC auxin
#################################

dataset = 'HiC'
dirr = main_dirr+dataset+'/'
source = dirr+'auxin_500.csv'
target = dirr + 'Dory_auxin_'

print('#########################')
print('Processing', dataset)
print('#########################')

thresh = 400
filetype = 2
dim = 2

start = time.time()
dory.compute_PH(source, lower_thresh, thresh, filetype, threads, target, dim, compute_cycles, reduce_cyc_lengths, suppress_output)
print('Time taken for', dataset, ':', time.time() - start)

