import numpy as np
import matplotlib.pyplot as plt
import math
import networkx as nx
import pickle


def rot_45(mat):

    ang = -math.pi/4

    rot_mat = np.array([[math.cos(ang), -math.sin(ang)]\
                        ,[math.sin(ang), math.cos(ang)]])

    for idx, row in enumerate(mat):
        mat[idx,:] = np.matmul(rot_mat, row.T)

    return mat


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
chrom1_start_bin = chrom_start_bin_dict['chr1']
chrom1_end_bin = chrom_start_bin_dict['chr2']

print(chrom1_start_bin, chrom1_end_bin)


# Plot the edges
edge_file = open('control/edges.csv', 'r')
edge_pts = []

count = 0
for line in edge_file:
    print(count, end='\r')
    count += 1
    line = line.split(',')

    bin1 = int(line[0])
    bin2 = int(line[1])

    #if abs(bin1 - bin2) < 5:
    #    continue

    if bin1 >= chrom1_start_bin and bin1 <= chrom1_end_bin:
        if bin2 >= chrom1_start_bin and bin2 <= chrom1_end_bin:


            dist = float(line[2])
            if dist >= 150:
                continue

            edge_pts.append([bin1, bin2])

edge_pts = np.array(edge_pts)
edge_pts = rot_45(edge_pts)

plt.scatter(edge_pts[:,0], edge_pts[:,1], alpha=0.2, color='blue')

short_cycs = open('control/diff_thresholds/thresh_150_minimal_V_birth_H1.txt', 'r')

count = 0

all_pts = []

bin_cyc_counts = dict()

for line in short_cycs:

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
        if x1 >= chrom1_start_bin and x1 <= chrom1_end_bin:
            nn += 1
            continue
        skip_trans_cycle = 1
        break
    
    if skip_trans_cycle:
        #print('SKIPPING')
        continue
    
    #print(cyc)

    G = nx.Graph()

    nn = 0

    while nn < lenn:

        x1 = cyc[nn]
        y1 = cyc[(nn+1)%lenn]

        G.add_edge(x1, y1)

        if x1 not in bin_cyc_counts:
            bin_cyc_counts[x1] = 1
            all_measures[x1-chrom1_start_bin, 0] = 1
        else:
            bin_cyc_counts[x1] += 1
            all_measures[x1-chrom1_start_bin, 0] += 1

        if y1 not in bin_cyc_counts:
            bin_cyc_counts[y1] = 1
            all_measures[y1-chrom1_start_bin, 0] = 1
        else:
            bin_cyc_counts[y1] += 1
            all_measures[y1-chrom1_start_bin, 0] += 1

        nn += 2

    cyc_basis = nx.cycle_basis(G)


    for cyc in cyc_basis:

        pts = []

        nn = 0

        while (nn < len(cyc)):

            x1 = cyc[nn]
            y1 = cyc[(nn+1)%len(cyc)]


            if y1 < x1:
                c = y1
                y1 = x1
                x1 = c

            pts.append([x1, y1])

            nn += 1

        
        pts = np.array(pts)

        pts = rot_45(pts)

        plt.plot(pts[:,0], pts[:,1], color='red', alpha = 0.2)



plt.show()

plt.cla()
plt.clf()
plt.close()

