import numpy as np
import networkx as nx
import cooler



labels = [\
         'fa_dsg_ddel_5kb'\
        ,'fa_dsg_dpn_5kb'\
        ,'fa_dsg_ddel_dpn_5kb'\
        ,'fa_dsg_mnase_5kb'\
        ,'fa_dpn_5kb'\
        ]

for e_idx, exp in enumerate(labels):

    print(exp)

    cfile = '../cooler_files/'+exp + '.cool'

    cf = cooler.Cooler(cfile)

    mat = cf.matrix()

    print(mat[4000:5000, 4000:5000])

    #edge_file = '../cooler_results/edge_files/'+exp+'_edges_ALL.csv'

    #edge_data = np.loadtxt(edge_file, delimiter=',')


    ##print(edge_data)
    ##exit()

    #G = nx.Graph()

    #graph_file = exp+'_edge_graph_5kb.gml'

    #for ee in edge_data:
    #    b1 = int(ee[0])
    #    b2 = int(ee[1])
    #    weight = 1/float(ee[2])

    #    #print(b1, b2, weight)

    #    G.add_edge(b1, b2, weight=weight)

    #nx.write_gml(G, graph_file)

        
    


