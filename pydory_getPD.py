import pydory as dory
import sys
import time


main_dirr = 'Datasets/'

# Main command:
# dory.compute_PH(source, lower_thresh, thresh\
#            , filetype, threads, target, dim\
#            , compute_cycles, reduce_cyc_lengths, cyc_thresh, suppress_output\
#            , hom_batch_size, cohom_batch_size)

# source: string, Input file in csv format. Details in filetype below.

# lower_thresh: float. Experimental. Best to set at 0.

# thresh: float. Threshold for PH computation. Set to -1 to consider all possible simplices on the point-cloud.

# filetype --- 0, 1, 2, 3
# Filetype details
# Accepts comma-separated files
# 0: Distance square matrix
# 1: Point-cloud (locations in N-dimensional space) in the format --- x1, x2, ..., xN
# 2: List of edges in the format --- v1, v2, edge length
# 3: Sorted list of edges in the format --- v1, v2, edge length, sorted in increasing order of edge length

# threads --- integer, number of threads for serial-parallel computation, generally 4.

# target --- string, prefix for location of saving the results
# Output files for PD
# H1 persistence diagram is stored as a txt file: target+'H1_pers_data.txt'
# H2 persistence diagram is stored as a txt file: target+'H2_pers_data.txt'

# dim --- 1, 2.
# dim = 1 for up to and including H1
# dim = 2 for up to and including H2

# compute_cycles --- 0, 1.
# compute_cycles = 0 if only want PD.
# compute_cycles = 1 if want representative boundaries.

# reduce_cyc_lengths --- 0, 1.
# reduce_cyc_lengths = 0 to ignore greedy shortening.
# reduce_cyc_lengths = 1 for greedy shortening.
# Note: Set compute_cycles and reduce_cyc_lengths to same value.

# cyc_thresh: float. Threshold for cycle computation.
# Only cycles born at parameter <= cyc_thresh will be considered for shortening
# and saved as target+'minimal_V_birth_H1.txt', target+'minimal_V_birth_H2.txt'
# Each line of the saved file is a cycle.
# Each H1 cycle is stored as v0, v1, v2, v3, ..... where (v_2i, v_(2i+1)) is an edge in the H1-cycle.
# Each H2 cycle is stored as v0, v1, v2, v3, ..... where (v_3i, v_(3i+1), v_(3i+2)) is a triangular face in the H2-cycle.

# suppress_output: 0, 1
# suppress_output = 0 to show progress.
# suppress_output = 1 to hide progress.
# Note: Progress is to give a rough estimate of how long it might take to finish computation, and does not imply linear in time progression.

# hom_batch_size: integer
# Batch size for serial-parallel computation of H0 computation. Suggested value 1000.

# cohom_batch_size: integer
# Batch size for serial-parallel computation of coH1 and coH2 computation. Suggested value 100.


compute_cycles = 0
reduce_cyc_lengths = 0
suppress_output = 1
threads = 4
lower_thresh = 0

hom_ws = 1000
cohom_ws = 100




#################################
## dragon (point-cloud, H1)
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

dory.compute_PH(source, lower_thresh, thresh, filetype, threads, target, dim, compute_cycles, reduce_cyc_lengths, thresh, suppress_output, hom_ws, cohom_ws)
print('Time taken for', dataset, ':', time.time() - start)


#################################

#################################
## fract (distance matrix, H2)
#################################
dataset = 'Fract'
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
dory.compute_PH(source, lower_thresh, thresh, filetype, threads, target, dim, compute_cycles, reduce_cyc_lengths, thresh, suppress_output, hom_ws, cohom_ws)
print('Time taken for', dataset, ':', time.time() - start)



#################################
## o3 (point-cloud, H2)
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
dory.compute_PH(source, lower_thresh, thresh, filetype, threads, target, dim, compute_cycles, reduce_cyc_lengths, thresh, suppress_output, hom_ws, cohom_ws)
print('Time taken for', dataset, ':', time.time() - start)

#################################
## torus4 (point-cloud, H1)
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
dory.compute_PH(source, lower_thresh, thresh, filetype, threads, target, dim, compute_cycles, reduce_cyc_lengths, thresh, suppress_output, hom_ws, cohom_ws)
print('Time taken for', dataset, ' (H1):', time.time() - start)


#################################
## torus4 (point-cloud, H2)
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
dory.compute_PH(source, lower_thresh, thresh, filetype, threads, target, dim, compute_cycles, reduce_cyc_lengths, thresh, suppress_output, hom_ws, cohom_ws)
print('Time taken for', dataset, ' (H2):', time.time() - start)

#################################
## HiC control (list of edges, H2)
#################################

dataset = 'HiC'
dirr = main_dirr+dataset+'/'
source = dirr+'control_400.csv'
target = dirr + 'Dory_control_'

print('#########################')
print('Processing', dataset)
print('#########################')

thresh = 400
filetype = 2
dim = 2

start = time.time()
dory.compute_PH(source, lower_thresh, thresh, filetype, threads, target, dim, compute_cycles, reduce_cyc_lengths, thresh, suppress_output, hom_ws, cohom_ws)
print('Time taken for', dataset, ':', time.time() - start)

#################################
## HiC auxin (list of edges, H2)
#################################

dataset = 'HiC'
dirr = main_dirr+dataset+'/'
source = dirr+'auxin_400.csv'
target = dirr + 'Dory_auxin_'

print('#########################')
print('Processing', dataset)
print('#########################')

thresh = 400
filetype = 2
dim = 2

start = time.time()
dory.compute_PH(source, lower_thresh, thresh, filetype, threads, target, dim, compute_cycles, reduce_cyc_lengths, thresh, suppress_output, hom_ws, cohom_ws)
print('Time taken for', dataset, ':', time.time() - start)

