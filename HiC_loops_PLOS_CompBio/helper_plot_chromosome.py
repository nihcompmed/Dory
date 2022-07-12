from numba import njit
import math
from matplotlib.patches import Arc
from matplotlib.collections import LineCollection
import matplotlib.pyplot as plt
import random
import numpy as np
import networkx as nx
import itertools as it
from matplotlib.pyplot import cm



@njit
def polarize(n, r, circum):
    theta = n/r
    x = r*math.cos(theta)
    y = r*math.sin(theta)
    return r*math.cos(theta), r*math.sin(theta)

@njit
def rot45(mat):

    n_rows, n_col = mat.shape
    
    for i in range(n_rows):
        x = mat[i, 0]
        y = mat[i, 1]

        x_rot = x*cos45 + y*sin45
        y_rot = -x*sin45 + y*cos45

        mat[i, 0] = x_rot
        mat[i, 1] = y_rot

    return mat

@njit
def rot_theta(x, y, theta):
    return x*math.cos(theta) - y*math.sin(theta)\
            ,x*math.sin(theta) + y*math.cos(theta)

def bez_curve(x1, x2, y1, y2, xb, yb):
    bezier_path = np.arange(0, 1.01, 0.01)

    # Compute and store the Bezier curve points
    x = (1 - bezier_path)** 2 * x1 + 2 * (1 - bezier_path) * bezier_path * xb + bezier_path** 2 * x2
    y = (1 - bezier_path)** 2 * y1 + 2 * (1 - bezier_path) * bezier_path * yb + bezier_path** 2 * y2

    path = []

    for i in range(len(x)):
        path.append((x[i], y[i]))

    #bezier_paths.append(path)

    #bezier_colors.append(color)

    #bezier_colors.append(mpl.colors.to_rgba(color, alpha))

    #bezier_style.append(style)

    #plt.plot(x, y, c = color, linestyle = style, alpha=alpha)

    return path



def bez_curve_polarized(bin1, bin2, bin_coords, color, style, alpha, ax):

        x1 = bin_coords[bin1][0]
        y1 = bin_coords[bin1][1]

        if x1:
            theta1 = math.atan(y1/x1)
        else:
            theta1 = math.pi/2

        x1_flat, y1_flat = rot_theta(x1, y1, 2*math.pi - theta1)

        x2 = bin_coords[bin2][0]
        y2 = bin_coords[bin2][1]
        if x2:
            theta2 = math.atan(y2/x2)
        elif y2 > 0:
            theta2 = math.pi/2
        else:
            theta2 = -math.pi/2

        x2_flat, y2_flat = rot_theta(x2, y2, 2*math.pi - theta2)


        theta_mid = (theta1 + theta2)/2

        x_b = (x1_flat + x2_flat)/2
        y_b = np.abs(x1_flat - x2_flat)/2

        xb_rot, yb_rot = rot_theta(x_b, y_b, theta_mid)

        path = bez_curve(x1, x2, y1, y2, xb_rot, yb_rot)

        return path




def plot_chrom_border(ax, chrom_start, chrom_end):


    circum = chrom_end + 1 - chrom_start

    r_chrom = circum/(2*math.pi)

    width = 0.15

    r_inner = (1-width)*r_chrom
    r_outer = (1+width)*r_chrom

    #r_inner = 0.25
    #r_outer = 0.75

    print(r_inner, r_outer)

    # -80 to 260
    angle_offset = 80

    bin_coords = dict()

    # Get bin coords
    for binn in range(chrom_start, chrom_end + 1):

        bin_polar = (binn - chrom_start)/(chrom_end+1-chrom_start)*(260+80) - 80

        # pick radius randomly
        r_bin = random.uniform(r_inner, r_outer)

        bin_coords[binn] = [r_bin*math.cos(bin_polar*math.pi/180), r_bin*math.sin(bin_polar*math.pi/180)]
        

    # Plot the chromosome border
    ax.add_patch(Arc((0, 0), 2*r_inner, 2*r_inner,
         theta1=-angle_offset, theta2=260, edgecolor='k'\
                 , lw = 1, ls = '--', alpha = 0.75))
    ax.add_patch(Arc((0, 0), 2*r_inner, 2*r_inner,
         theta1=-angle_offset, theta2=260, edgecolor='k'\
                 , lw = 4, ls = '--', alpha = 0.5))

    ax.add_patch(Arc((0, 0), 2*r_outer, 2*r_outer,
         theta1=-angle_offset, theta2=260, edgecolor='k'\
                 , lw = 1, ls = '--', alpha = 0.75))
    ax.add_patch(Arc((0, 0), 2*r_outer, 2*r_outer,
         theta1=-angle_offset, theta2=260, edgecolor='k'\
                 , lw = 4, ls = '--', alpha = 0.5))

    ax.set_xlim([-1.25*r_chrom, 1.25*r_chrom])
    ax.set_ylim([-1.25*r_chrom, 1.25*r_chrom])

    # Plot the chromosome end and start lines

    # Start
    x1 = r_inner*math.cos(-angle_offset*math.pi/180)
    y1 = r_inner*math.sin(-angle_offset*math.pi/180)

    x2 = r_outer*math.cos(-angle_offset*math.pi/180)
    y2 = r_outer*math.sin(-angle_offset*math.pi/180)

    plt.plot([x1, x2], [y1, y2], color='black')

    # Start
    x1 = r_inner*math.cos(260*math.pi/180)
    y1 = r_inner*math.sin(260*math.pi/180)

    x2 = r_outer*math.cos(260*math.pi/180)
    y2 = r_outer*math.sin(260*math.pi/180)

    plt.plot([x1, x2], [y1, y2], color='black')



    # Plot the ticks and the tick labels
    n_ticks = 10
    #binn = begin[chrom_idx] - bin_ref
    binn = 0
    step = int((chrom_end - chrom_start)/n_ticks)
    for i in range(n_ticks):
        theta = binn/r_chrom
        xx = r_inner*math.cos(theta-math.pi*angle_offset/180)
        yy = r_inner*math.sin(theta-math.pi*angle_offset/180)
        plt.scatter(xx, yy, color='k', s=4)
        plt.text(0.88*xx, 0.88*yy, str(binn)\
                , ha = 'center'\
                , va = 'center'\
                )
        binn += step
    

    ax.axis('off')

    #plt.show()

    return ax, bin_coords



def plot_gene_GO(genes, dict_GO, gene_to_GO, genes_to_bins):

        plt.cla()
        plt.clf()
        plt.close()

        all_GO = []
        G = nx.Graph()

        for pair in it.combinations(genes, 2):
            pair_graph = nx.Graph()

            if pair[0] not in gene_to_GO or pair[1] not in gene_to_GO:
                continue

            pair_graph.add_node(pair[0])
            pair_graph.add_node(pair[1])

            for GO in gene_to_GO[pair[0]]['GO']:
                pair_graph.add_edge(pair[0], GO)
            for GO in gene_to_GO[pair[1]]['GO']:
                pair_graph.add_edge(pair[1], GO)


            all_GO = gene_to_GO[pair[0]]['GO']\
                    +gene_to_GO[pair[1]]['GO']
            all_GO = frozenset(all_GO)
            for GO in all_GO:
                for GO_related in dict_GO[GO]['related_GO']:
                    if GO_related in all_GO:
                        pair_graph.add_edge(GO, GO_related)

            paths = nx.all_simple_paths(pair_graph, pair[0], pair[1])
            for path in paths:
                nx.add_path(G, path)


        layout = nx.bipartite_layout(G, genes)
        gene_color = {}
        count_gene = 0
        layout_dict = {}
        cmap = cm.get_cmap('tab10')
        cmap = cmap.colors

        #print(layout)
        for node in layout:
            pos = layout[node]
            xx = pos[0]
            yy = pos[1]
            layout_dict[node] = {'x':xx, 'y':yy}
            if node in genes:
                gene_color[node] = cmap[count_gene]
                count_gene += 1
                color = gene_color[node]
                gene_annotate = node + '\n' + str(genes_to_bins[node]['chrom'])\
                            +'\n'+str(min(genes_to_bins[node]['bins']))\
                            +'\n'+str(max(genes_to_bins[node]['bins']))
                if yy >= 0:
                    plt.text(xx-0.75, yy+0.05, gene_annotate, fontweight='bold', color=color)
                elif yy < 0:
                    plt.text(xx-0.75, yy-0.05, gene_annotate, fontweight='bold', color=color)
                #plt.scatter(xx, yy, c = color, zorder = 100)
            else:
                plt.text(xx+1, yy-0.005, dict_GO[node]['name'], fontsize=14, wrap=True)
                plt.plot([xx, xx+1], [yy, yy], c = 'k', alpha = 0.5, linestyle=':', zorder = 100)
                plt.scatter(xx, yy, c = 'k', alpha = 0.75, s = 6, zorder = 100)


        for edge in G.edges:
            if edge[0] in genes or edge[1] in genes:
                if edge[0] in genes:
                    color = gene_color[edge[0]]
                else:
                    color = gene_color[edge[1]]
                linestyle = '-'
                opacity = 1
                plt.plot([layout_dict[edge[0]]['x'], layout_dict[edge[1]]['x']]\
                        , [layout_dict[edge[0]]['y'], layout_dict[edge[1]]['y']]\
                        , linestyle = linestyle\
                        , alpha = opacity\
                        , color = color)
            else:
                linestyle = '--'
                opacity = 0.5
                color = 'green'
                x1 = layout_dict[edge[0]]['x']
                x2 = layout_dict[edge[1]]['x']
                y1 = layout_dict[edge[0]]['y']
                y2 = layout_dict[edge[1]]['y']

                bez_curve(x1, x2\
                        , y1, y2\
                        , x1 + 1, (y1+y2)/2
                        )
            
            #print('EDGE', edge[0], edge[1])
            #input('w')
        #plt.figure(figsize())
        ax = plt.gca()
        ax.xaxis.set_visible(False)
        ax.yaxis.set_visible(False)
        ax.axis('off')
        plt.xlim([-2, 4])
        #plt.margins(0,0)
        plt.gca().xaxis.set_major_locator(plt.NullLocator())
        plt.gca().yaxis.set_major_locator(plt.NullLocator())
        #plt.title(list(chrom)[0])
        
        plt.show()
        #print(list(chrom)[0], cycle_counter)
        #plt.savefig(pwd+'/gene_graphs/'+list(chrom)[0]+cycle_type+str(cycle_counter)+'.pdf', format='pdf')
        #input('w')
        plt.clf()

        #cycle_counter +=1 
        #plt.show()



        #print(layout)
        #input('w')
        #GO_subgraph = hv.Graph.from_networkx(sub_G, layout).opts(tools=['hover'])
        #GO_subgraph = GO_subgraph.options(width=1000, height=1000)
        #show(hv.render(GO_subgraph))
        #input('w')


