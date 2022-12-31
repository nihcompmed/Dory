import numpy as np
import networkx as nx
import itertools as it
import subprocess
import os



# Check if boundary of three simplices is connected, accept list of 3-simplices [[v1, v2, v3], ...]
# Return components
def connected_2simplices(boundary):

    G = nx.Graph()
    nn = 0
    while nn < len(boundary):

        G.add_edge(int(boundary[nn]), int(boundary[nn+1]))

        nn += 2


    #print(boundary)
    # simply return the cycle basis

    cycs = nx.cycle_basis(G)

    return cycs

    #print(cycs)

    #exit()




    #
    #L = nx.line_graph(G)

    #comps = nx.connected_components(L)

    #comps_list = []

    #for comp in comps:
    #    this_comp = []
    #    for c in comp:
    #        this_comp += list(c)

    #    comps_list.append(this_comp)
    #
    #return comps_list


# Check if boundary of three simplices is connected, accept list of 3-simplices [[v1, v2, v3], ...]
# Return components
def connected_3simplices(boundary):

    G = nx.Graph()

    nn = 0
    while nn < len(boundary):

        G.add_node(frozenset({boundary[nn], boundary[nn+1], boundary[nn+2]}))

        nn += 3

    print('Number of nodes', len(list(G.nodes)))

    
    # Computational costly
    # Iterating over all pairs of triangular faces
    for t1, t2 in it.combinations(list(G.nodes), 2):

        if len(t1.intersection(t2)) == 2:
            G.add_edge(t1, t2)

    comps = nx.connected_components(G)

    comps_list = []

    for comp in comps:

        this_comp = []
        
        for t in comp:
            this_comp += list(t)

        comps_list.append(this_comp)

    return comps_list


def modify_minimal(filepath, scratch_dir, dim):

    
    # Load cycle file
    ff = open(filepath, 'r')

    # THIS IS NOT MEMORY EFFICiENT
    scratch_filepath = scratch_dir+'scratch.txt'
    scratch_file = open(scratch_filepath, 'w')

    count = 0
    print('\n')

    for line in ff:

        print(count, end='\r')
        count += 1
    
        line = line.split(',')
        line = line[:-1]

        if dim == 1:
            comps = connected_2simplices(line)
        elif dim == 2:
            comps = connected_3simplices(line)
        else:
            print("ERROR IN DIMENSION IN CHECK CONNECTEDNESS")
            exit()
        
        #if (len(comps)>1):
        #    print('Disconnected boudnaries found')
        
        for boundary in comps:

            for vv in boundary:
                scratch_file.write(str(vv)+',')

            #boundary = ','.join(boundary)
            
            #scratch_file.write(boundary+',\n')
            scratch_file.write('\n')


    scratch_file.close()

    ff.close()

    
    # Remove original cycle file
    subprocess.run(['rm', filepath])

    # Rename scratch file to original cycle filename
    subprocess.run(['mv', scratch_filepath, filepath])

    print('\n')
    
