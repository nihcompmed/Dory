import pydory as dory
import check_connected_helper as cch
import helper_functions as hf 
import networkx as nx
import matplotlib.pyplot as plt


# Set experiment name
exp_name = 'control'
#exp_name = 'auxin'

source = exp_name+'/edges.csv'
target = exp_name+'/'

filetype = 2

dim = 1
threads = 4
thresh = 100
epsil = 50

# Can pick all thresh as 150 because edge file has 150 thresh
dory.compute_PH(source, 0, thresh+epsil, filetype, threads, target, dim, 1, 1, thresh, 1)

# Modify for connectedness
dim = 1
cyc_file = target + 'minimal_V_birth_H'+str(dim)+'.txt'
cch.modify_minimal(cyc_file, target, dim)


# Smoothen 
# 1. Get graph
filetype = 2
G = hf.get_edge_graph(source, filetype, thresh)

# 2. Smoothen with trivial triangles
dim = 1
# H1 boundary
ff = open(target+'minimal_V_birth_H'+str(dim)+'.txt', 'r')

# New file that will have smoothed boundary
gg = open(target+'smoothen_H'+str(dim)+'.txt', 'w')

count = 0
for line in ff:
    print(count, end='\r')
    count += 1

    line = line.split(',')
    line = line[:-1]
    line = [int(x) for x in line]

    line = hf.smoothen_trivial_triangles(line, G)

    for x in line:
        gg.write(str(x)+',')
    gg.write('\n')

ff.close()
gg.close()
