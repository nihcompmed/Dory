import numpy as np
import networkx as nx
import itertools





# 0. Select experiment
exp = 'control'
exp = 'auxin'

edge_file = exp+'/edges.csv'

ff = open(edge_file, 'r')

# 1. Make graph for control
G = nx.Graph()

count = 0

for line in ff:

    line = line.strip('\n')
    line = line.split(',')

    bin1 = int(line[0])
    bin2 = int(line[1])
    dist = float(line[2])

    G.add_edge(bin1, bin2, weight=dist)

ff.close()

# Read smoothed cycles
ff = open(exp + '/smoothen_H1.txt', 'r')

gg = open(exp + '/smooth_H1_estimated_pers.txt', 'w')

all_cycs = 0
non_degenerate = 0
non_trivial = 0

for line in ff:
    line = line.strip('\n')
    line = line.split(',')
    line = line[:-1]
    line = [int(x) for x in line]

    all_cycs += 1

    if len(line) < 3:
        #print('skipping degenerate')
        continue

    non_degenerate += 1

    nn = 0
    maxx = 0
    while nn < len(line):

        bin1 = line[nn]
        bin2 = line[(nn+1)%len(line)]
        maxx = max(maxx, G[bin1][bin2]['weight'])


        nn += 1

    birth = maxx


    maxx = 0
    for bin1, bin2 in itertools.combinations(line, 2):
        if not G.has_edge(bin1, bin2):
            maxx = -1
            break

        maxx = max(maxx, G[bin1][bin2]['weight'])


    death_estimate = maxx

    if death_estimate == -1 or death_estimate - birth >= 50:
        non_trivial += 1
        gg.write(str(birth) + ',' + str(death_estimate) + '\n')

        gg.write('smooth cycle,')

        nn = 0
        while nn < len(line):
            gg.write(str(line[nn]) + ',')
            nn += 1
        gg.write('\n')

print('all', all_cycs, 'non degenerate', non_degenerate, 'possibly significant', non_trivial)

gg.close()

ff.close()



