import numpy as np
import matplotlib.pyplot as plt
import math
import networkx as nx
import pickle
import pyBigWig
from matplotlib.collections import LineCollection
import helper_plot_chromosome as hpc


chrom_sizes = open('4DNFI823LSII.chrom.sizes')
chrom_list = []
chrom_lenn = dict()
for line in chrom_sizes:
    line = line.strip('\n')
    line = line.split('\t')
    chrom = line[0]
    lenn = int(int(line[1])/1000)+1
    chrom_lenn[chrom] = lenn
    chrom_list.append([chrom, lenn])

chrom_start_bin = []

start_bin_counter = 0

for idx in range(len(chrom_list)):
    chrom = chrom_list[idx][0]

    chrom_start_bin.append([chrom, start_bin_counter])

    start_bin_counter += chrom_list[idx][1]


chrom_start_bin_dict = dict()

for entry in chrom_start_bin:
    chrom = entry[0]
    start_bin = entry[1]

    chrom_start_bin_dict[chrom] = start_bin


#Test chrom 1
chrom_start_bin = chrom_start_bin_dict['chr1']
chrom_end_bin = chrom_start_bin_dict['chr2']

plt.figure(figsize=(10, 10))

ax = plt.gca()

ax, bin_coords = hpc.plot_chrom_border(ax, chrom_start_bin, chrom_end_bin)


# Plot the cycles
short_cycs = open('control/diff_thresholds/thresh_150_minimal_V_birth_H1.txt', 'r')

cyc_paths = []

count = 0

for line in short_cycs:

    print('cyc number', count, end='\r')
    count += 1

    cyc = line.split(',')
    cyc = cyc[:-1]


    cyc = [int(x) for x in cyc]

    lenn = len(cyc)

    nn = 0

    # Check cyc in chr1
    skip_trans_cycle = 0
    while nn < lenn:

        x1 = cyc[nn]
        if x1 >= chrom_start_bin and x1 <= chrom_end_bin:
            nn += 1
            continue
        skip_trans_cycle = 1
        break
    
    if skip_trans_cycle:
        #print('SKIPPING')
        continue
    

    G = nx.Graph()

    nn = 0

    max_bin_dist = 0

    while nn < lenn:

        # Short cycles info is of the form (v1, v2), (v3, v4), ...

        x1 = cyc[nn]
        y1 = cyc[(nn+1)%lenn]

        max_bin_dist = max(max_bin_dist, abs(x1-y1))

        G.add_edge(x1, y1)

        nn += 2

    # PLOTTING LOOPS WITH LONG-RANGE INTERACTION 
    if max_bin_dist < 10000:
        continue

    cyc_basis = nx.cycle_basis(G)

    for cyc in cyc_basis:

        nn = 0

        while (nn < len(cyc)):

            bin1 = cyc[nn]
            bin2 = cyc[(nn+1)%len(cyc)]

            path = hpc.bez_curve_polarized(bin1, bin2, bin_coords, 'green', '--', 0.8, ax)

            cyc_paths.append(path)

            nn += 1

        

bezier_collection = LineCollection(cyc_paths, pickradius=10, colors='red', linestyle='-', alpha=0.2)

bezier_collection.set_picker(True)

ax.add_collection(bezier_collection)


plt.show()

