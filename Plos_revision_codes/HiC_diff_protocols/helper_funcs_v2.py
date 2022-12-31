import matplotlib.pyplot as plt
import numpy as np
import pandas
import math
import h5py
from numba import njit, prange
from numba import njit
import pickle
import itertools
import cooler 
import pyBigWig
from sklearn.neighbors import KernelDensity
import networkx as nx

@njit
def rot45(mat):

    n_rows, n_col = mat.shape

    cos45 = 1/math.sqrt(2)
    sin45 = cos45
    
    for i in range(n_rows):
        x = mat[i, 0]
        y = mat[i, 1]

        x_rot = x*cos45 + y*sin45
        y_rot = -x*sin45 + y*cos45

        mat[i, 0] = x_rot
        mat[i, 1] = y_rot

    return mat


@njit
def search_edges(n_elements, bin1, bin2, pixels, weights, thresh, edges\
                , count, chrom_start_bins):
    
    for i in range(n_elements):
        val = pixels[i]
        bin1_e = bin1[i]
        bin2_e = bin2[i]

        if not (bin1_e < bin2_e) :
            continue
        c1_idx = 0
        while (chrom_start_bins[c1_idx] < bin1_e):
            c1_idx += 1

        c2_idx = 0
        while (chrom_start_bins[c2_idx] < bin2_e):
            c2_idx += 1

        if (c1_idx != c2_idx):
            continue

        c_idx = c1_idx - 1


        bal_val = val*weights[bin1_e]*weights[bin2_e]
        bal_val = 1/bal_val

        if (bal_val < thresh):
            edges[c_idx, count[c_idx], 0] = bin1_e
            edges[c_idx, count[c_idx], 1] = bin2_e
            edges[c_idx, count[c_idx], 2] = bal_val
            count[c_idx] += 1

    return edges, count

@njit
def save_all_edges(jj, bin1, bin2, pixels, all_edges, count):

    for i in range(jj):
        val = pixels[i]
        bin1_e = bin1[i]
        bin2_e = bin2[i]

        # +1 for MaxHiC interactions file indexing
        all_edges[count, 0] = bin1_e+1
        all_edges[count, 1] = bin2_e+1
        all_edges[count, 2] = val

        count += 1

    return all_edges, count







def get_lens_birth_cycs(fname):

    exp1_len = open(fname, 'r')
    line = exp1_len.readline()
    line = line.split(',')
    line = line[:-1]
    line = [int(x) for x in line]
    exp1_len.close()

    return line

def get_cycs_as_edge_lists(fname):

    cycs = open(fname, 'r')

    edges_list = []

    for line in cycs:
        line = line.split(',')
        line = line[:-1]
        line = [int(x) for x in line]
        if len(line) < 3:
            continue
        nn = 0
        cyc = []
        while nn < len(line):
            v1 = line[nn]
            v2 = line[(nn+1)%len(line)]

            g1 = min(v1, v2)
            g2 = max(v1, v2)

            cyc.append([g1, g2])
            nn += 1

        edges_list.append(cyc)

    cycs.close()
    
    return edges_list

def show_graph(G, edge_list):

    for cyc in edge_list:

        for (v1, v2) in cyc:
            G.add_edge(v1, v2)

    return G

def make_line_graph_plot(ax, edge_list, color, rev, scale_reso):


    for cyc in edge_list:

        xx = []
        yy = []

        for (v1, v2) in cyc:

            if not rev:
                xx.append(min(v1, v2))
                yy.append(max(v1, v2))
            else:
                xx.append(max(v1, v2))
                yy.append(min(v1, v2))

        (v1, v2) = cyc[0]
        if not rev:
            xx.append(min(v1, v2))
            yy.append(max(v1, v2))
        else:
            xx.append(max(v1, v2))
            yy.append(min(v1, v2))

        xx = np.array(xx)*scale_reso
        yy = np.array(yy)*scale_reso
        ax.plot(xx, yy, lw=5, color = color, alpha=0.4)


    return ax



def get_hiccup_peaks(fname):

    ff = open(fname, 'r')

    # read header
    ff.readline()

    unique_reso = []

    peaks_dict = dict()

    for line in ff:
        line = line.split('\t')

        #print(line)

        chrom = line[0]

        b1 = int(line[1])
        b2 = int(line[2])

        reso = abs(b1 - b2)

        if reso not in unique_reso:
            unique_reso.append(reso)


        cent1 = min(float(line[6]), float(line[7]))
        cent2 = max(float(line[6]), float(line[7]))


        start1 = int(line[1])
        end1 = int(line[2])

        start2 = int(line[4])
        end2 = int(line[5])

        c_size = int(line[9])
        c_count = int(line[10])

        c_donut = float(line[11])
        c_vertical = float(line[12])
        c_horizontal = float(line[13])
        c_lowleft = float(line[14])

        #print(c_donut, c_vertical, c_horizontal, c_lowleft)
        #exit()

        if chrom not in peaks_dict:
            peaks_dict[chrom] = [[cent1, cent2\
                            , c_size, c_count\
                            , c_donut, c_vertical\
                            , c_horizontal, c_lowleft]]
        else:
            peaks_dict[chrom].append([cent1, cent2\
                            , start1, end1\
                            , start2, end2\
                            , c_size, c_count\
                            , c_donut, c_vertical\
                            , c_horizontal, c_lowleft])
                    

    ff.close()

    return peaks_dict

def get_cycs_hiccups_dist(cycs, hiccups, reso):

    n_cycs = len(cycs)
    n_hiccups = hiccups.shape[0]

    for cyc_idx, cyc in enumerate(cycs):

        print(cyc_idx, end='\r')

        all_measures = []
        all_births = []
        all_deaths = []
        all_genomic_dist = []
        
        for cyc in cycs:
        
            cyc = np.array(cyc)*reso
        
            cyc_measure = math.inf
        
            cyc_birth = 0
        
            genomic_dist = 0
        
            all_verts = []
        
            for row in cyc:
        
        
                v1 = row[0]
                v2 = row[1]
        
                genomic_dist = max(genomic_dist, abs(v2 - v1))
        
                #all_verts += [v1, v2]
        
                #this_edge = global_G[v1][v2]['weight']
        
                #cyc_birth = max(cyc_birth, this_edge)
        
                temp_diff = np.abs(hiccups - row)
        
                maxx = np.max(temp_diff, axis=1)
        
                minn = np.min(maxx)
        
                cyc_measure = min(minn, cyc_measure)
        
            all_measures.append(cyc_measure)
            #all_births.append(cyc_birth)
            all_genomic_dist.append(genomic_dist)
        
            #all_verts = list(frozenset(all_verts))
        
            #cyc_death = -1
            #for u, v in itertools.combinations(all_verts, 2):
        
            #    if not global_G.has_edge(u, v):
            #        continue
        
            #    cyc_death = max(cyc_death, global_G[u][v]['weight'])
            #
            #all_deaths.append(cyc_death)

    print(n_cycs)
    print(n_hiccups)


    all_measures = np.array(all_measures)/reso
    all_genomic_dist = np.array(all_genomic_dist)/reso

    plt.scatter(all_genomic_dist, all_measures, alpha=0.5, color='blue')
    plt.xlabel('maximal genomic interaction in PH loop (5kb)')
    plt.ylabel('nearest hiccup peak to PH loop')

    plt.yscale('log', base=2)
    plt.xscale('log', base=2)

    plt.show()



    #fig, axs = plt.subplots(1, 3)
    #
    #all_births = np.array(all_births)
    #all_deaths = np.array(all_deaths)
    #all_measures = np.array(all_measures)
    #all_genomic_dist = np.array(all_genomic_dist)
    #
    #all_measures = all_measures/PH_reso
    #all_genomic_dist = all_genomic_dist/PH_reso
    #
    #all_pers = all_deaths - all_births
    #
    #axs[0].scatter(all_pers, all_measures)
    #axs[1].scatter(all_genomic_dist, all_measures)
    #axs[2].scatter(all_genomic_dist, all_pers)
    #
    #axs[0].set_xlabel('estimated persistence')
    #axs[0].set_ylabel('measure')
    #
    #axs[1].set_xlabel('max genomic dist')
    #axs[1].set_ylabel('measure')
    #
    #axs[2].set_xlabel('max genomic dist')
    #axs[2].set_ylabel('estimated persistence')
    #
    #plt.show()
    

def get_hiccup_peaks_dict(fname, get_chrom):

    ff = open(fname, 'r')

    # read header
    ff.readline()


    peaks = dict()

    for line in ff:
        line = line.split('\t')
        chrom = line[0]
        if chrom != get_chrom:
            continue

        b1 = int(line[1])
        b2 = int(line[2])

        #reso = abs(b1 - b2)

        #if reso not in unique_reso:
        #    unique_reso.append(reso)


        cent1 = min(float(line[6]), float(line[7]))
        cent2 = max(float(line[6]), float(line[7]))

        cent1 = int(cent1)
        cent2 = int(cent2)+1

        c_size = int(line[9])
        c_count = int(line[10])

        c_donut = float(line[11])
        c_vertical = float(line[12])
        c_horizontal = float(line[13])
        c_lowleft = float(line[14])

        #print(c_donut, c_vertical, c_horizontal, c_lowleft)
        #exit()

        #peaks.append([cent1, cent2\
        #            , c_size, c_count\
        #            , c_donut, c_vertical\
        #            , c_horizontal, c_lowleft]\
        #            )

        peaks[(cent1, cent2)] = {\
                            'c_size':c_size\
                            ,'c_count':c_count\
                            ,'c_donut':c_donut\
                            ,'c_vertical':c_vertical\
                            ,'c_horizontal':c_horizontal\
                            ,'c_lowleft':c_lowleft\
                            }



    ff.close()

    #peaks = np.array(peaks)
    #print(peaks.shape)

    return peaks




def find_nearest_biomarker(binn, biomarker_dict):

    minn = math.inf
    nearest_bin = -1
    nearest_bio_peaks = []

    for key in biomarker_dict:

        dist = abs(binn - key)
        if dist < minn:
            minn = dist
            nearest_bin = [key]
            nearest_bio_peaks = [biomarker_dict[key]]
        elif dist == minn:
            nearest_bin.append(key)
            nearest_bio_peaks.append(biomarker_dict[key])

    return minn, nearest_bin, nearest_bio_peaks

def get_narrow_range(biomarker_narrow_file, chrom, start, end):

    # start and end are in bp

    ff = open(biomarker_narrow_file)

    xx = []
    yy = []

    bed_lines = []

    for line in ff:

        line = line.split('\t')
    
        this_chrom = line[0]


        if this_chrom != chrom:
            continue


        this_start = int(line[1])
        this_end = int(line[2])
        this_peak = int(line[-1])

        this_peak = this_start + this_peak

        if this_peak >= start and this_peak <= end:

            bed_lines.append(line)

            #signalval = float(line[-4])
            #xx.append(this_peak)
            #yy.append(signalval)

    #xx = np.array(xx)
    #yy = np.array(yy)

    return bed_lines

def plot_narrow_range(ax, biomarker_narrow_file, chrom, start, end):

    # start and end are in bp

    ff = open(biomarker_narrow_file)

    xx = []
    yy = []

    for line in ff:

        line = line.split('\t')
    
        this_chrom = line[0]

        if this_chrom != chrom:
            continue


        this_start = int(line[1])
        this_end = int(line[2])
        this_peak = int(line[-1])

        this_peak = this_start + this_peak

        if this_peak >= start and this_peak <= end:
            signalval = float(line[-4])
            xx.append(this_peak)
            yy.append(signalval)

    xx = np.array(xx)
    yy = np.array(yy)

    ax.scatter(xx/1000, yy)
    print(xx, yy)

    return

def plot_bigwig_range(ax, biomarker_bigwig_file, chrom, start, end):


    ff = pyBigWig.open(biomarker_bigwig_file)

    print(start, end)
    
    vals = ff.values(chrom, int(start), int(end))

    ax.plot(np.array(list(range(start, end)))/1000, vals)




def standardize_wrt_gdist1(edges):

    gdist_standard = 1

    # remove NaN from balanced
    nan_idxs = np.isnan(edges[:,2])
    not_nan_idxs = np.argwhere(nan_idxs==False).flatten()
    temp_edges = edges[not_nan_idxs]

    bin_dist = np.abs(temp_edges[:, 1] - temp_edges[:, 0])

    this_idxs = np.argwhere(bin_dist == gdist_standard).flatten()

    print('number of edges', len(this_idxs))

    this_edges = temp_edges[this_idxs, 2]

    #this_edges = np.log2(this_edges)
    outlier_ptile = np.percentile(this_edges, 99.9)

    trunc_edges = this_edges[this_edges <= outlier_ptile]

    mu = np.mean(trunc_edges)
    sigma = np.std(trunc_edges)

    edges[:,2] = (edges[:,2] - mu)/sigma

    nan_idxs = np.isnan(edges[:,2])
    not_nan_idxs = np.argwhere(nan_idxs==False).flatten()

    temp_edges = edges[not_nan_idxs, 2]
    positive_shift = np.amin(temp_edges)

    edges[:, 2] = edges[:,2] + abs(positive_shift)

    print('mu', mu, 'sigma', sigma, 'posshift', positive_shift)

    return edges

def standardize_w_mu_sigma_posshift(edges, mu, sigma, positive_shift):

    edges[:,2] = (edges[:,2] - mu)/sigma + abs(positive_shift)

    return edges


def KLdivergence(x, y):
  """Compute the Kullback-Leibler divergence between two multivariate samples.
  Parameters
  ----------
  x : 2D array (n,d)
    Samples from distribution P, which typically represents the true
    distribution.
  y : 2D array (m,d)
    Samples from distribution Q, which typically represents the approximate
    distribution.
  Returns
  -------
  out : float
    The estimated Kullback-Leibler divergence D(P||Q).
  References
  ----------
  Pérez-Cruz, F. Kullback-Leibler divergence estimation of
continuous distributions IEEE International Symposium on Information
Theory, 2008.
  """
  from scipy.spatial import cKDTree as KDTree

  # Check the dimensions are consistent
  x = np.atleast_2d(x)
  y = np.atleast_2d(y)

  n,d = x.shape
  m,dy = y.shape

  assert(d == dy)


  # Build a KD tree representation of the samples and find the nearest neighbour
  # of each point in x.
  xtree = KDTree(x)
  ytree = KDTree(y)

  # Get the first two nearest neighbours for x, since the closest one is the
  # sample itself.

  # p-norm = 1 for Manhattan
  # p-norm = 2 for Euclidean
  # p-norm = Infinity for Euclidean
  r = xtree.query(x, k=2, p=2)[0][:,1]
  s = ytree.query(x, k=1, p=2)[0]

  nonzero_r = np.argwhere(r!=0).flatten()

  #print(r/s)
  #print(np.amax(r/s), np.amin(r/s))

  #exit()

  AA = r/s

  # 0 can happen if there are multiple points at same location
  AA = AA[np.nonzero(r)]

  # q must have support bigger than p
  # q neq zero where p neq zero

  print('calculating A')
  A = -np.log(AA).sum() * d / n
  print('A is', A)

  print('calculating B')
  B = np.log(m / (n - 1.))
  print('B is', B)

  # There is a mistake in the paper. In Eq. 14, the right side misses a negative sign
  # on the first term of the right hand side.
  return A + B

def KL_symm(x,y):

    xy_div = KLdivergence(x,y)

    yx_div = KLdivergence(y,x)

    #return 0.5 * max(A + B,0.0)
    return 0.5 * (xy_div + yx_div)


def scale_hic_estimate_array(bal_val, mu, sigma, posshift, scale, shift):

    bal_val = np.power(((bal_val - mu)/sigma + posshift), scale) + shift

    return bal_val


def kde2D(x, y, bandwidth, xbins=100j, ybins=100j, **kwargs): 
    """Build 2D kernel density estimate (KDE)."""

    # create grid of sample locations (default: 100x100)
    xx, yy = np.mgrid[x.min():x.max():xbins, 
                      y.min():y.max():ybins]

    xy_sample = np.vstack([yy.ravel(), xx.ravel()]).T
    xy_train  = np.vstack([y, x]).T

    kde_skl = KernelDensity(bandwidth=bandwidth, **kwargs)
    kde_skl.fit(xy_train)

    # score_samples() returns the log-likelihood of the samples
    z = np.exp(kde_skl.score_samples(xy_sample))

    return xx, yy, np.reshape(z, xx.shape)

def get_chrom(binn, chrom_start_bins_list, chrom_names):

    for idx, start in enumerate(chrom_start_bins_list):

        if binn < start:
            return chrom_names[idx-1], binn-chrom_start_bins_list[idx-1]



def find_nearest_peak_bed(chrom, binn, reso, bed_file):

    ddist = math.inf

    target = binn*reso

    ff = open(bed_file, 'r')

    for line in ff:

        line = line.split('\t')
    
        this_chrom = line[0]

        if this_chrom != chrom:
            continue

        this_start = int(line[1])
        this_end = int(line[2])
        this_peak = int(line[-1])

        this_peak = this_start + this_peak

        ddist = min(ddist, abs(target-this_peak))

        #if this_peak >= start and this_peak <= end:
        #    signalval = float(line[-4])
        #    xx.append(this_peak)
        #    yy.append(signalval)

    #xx = np.array(xx)
    #yy = np.array(yy)

    #ax.scatter(xx/1000, yy)
    #print(xx, yy)

    ff.close()

    return ddist




def get_binary_tuples_combination(n):

    sub_result=list(itertools.product([0, 1], repeat=n))
    return sub_result


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



def plot_gene_GO(genes, dict_GO, gene_to_GO, genes_to_bins):

        plt.cla()
        plt.clf()
        plt.close()

        all_GO = []
        G = nx.Graph()

        for pair in itertools.combinations(genes, 2):
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
        cmap = plt.get_cmap('tab10')
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

def plot_loop_on_matrix(ax, cycle, cfile):

    print('plotting')

    cf = cooler.Cooler(cfile)

    mat = cf.matrix()
    
    nn = 0
    xx = []
    yy = []

    minn_bin = np.amin(np.array(cycle))
    maxx_bin = np.amax(np.array(cycle))

    this_mat = mat[minn_bin:maxx_bin+1, minn_bin:maxx_bin+1]

    while nn < len(cycle)-1:
        b1 = cycle[nn]
        b2 = cycle[(nn+1)]

        #print(b1, b2)
        xxx = b1
        yyy = b2

        xx.append(xxx)
        yy.append(yyy)

        nn += 1

    xx.append(xx[0])
    yy.append(yy[0])

    plot_mat = np.vstack((xx, yy)).T
    plot_mat = rot45(plot_mat)

    #print(cycle)

    ax.plot(plot_mat[:,0], plot_mat[:,1], color='black', lw=2)

    #plt.show()
    #return ax

    pts = []

    xx = []
    yy = []
    color = []

    for bb in range(minn_bin, maxx_bin+1):

        for cc in range(minn_bin, maxx_bin+1):

            if not math.isnan(this_mat[bb-minn_bin, cc-minn_bin]):

                weight = this_mat[bb-minn_bin, cc-minn_bin]
                xx.append(bb)
                yy.append(cc)
                color.append([1,0,0,weight])

    plot_mat = np.vstack((xx, yy)).T
    plot_mat = rot45(plot_mat)

    color = np.array(color)
    color[:,3] = color[:,3]/np.amax(color[:,3])

    ax.scatter(plot_mat[:,0], plot_mat[:,1], color=color)

    #plt.show()
    return ax

def get_chrom(binn, chrom_start_bins_list, chrom_names):

    for idx, start in enumerate(chrom_start_bins_list):

        if binn < start:
            return chrom_names[idx-1], binn-chrom_start_bins_list[idx-1]
