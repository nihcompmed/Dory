import numpy as np
from numba import njit
import os
import math
import matplotlib.pyplot as plt
import pickle
import networkx as nx
import itertools
import pydot
from networkx.drawing.nx_pydot import graphviz_layout
import copy
import sklearn
from sklearn.metrics import pairwise_distances

from scipy.spatial import cKDTree

from scipy import stats

#from scipy.stats import qmc

import gc
#import pydory as dory

#@njit
def get_cover_of_cycle(cycle, locs):

    cycle_locs = locs[cycle]

    minns = np.min(cycle_locs, axis=0)
    maxxs = np.max(cycle_locs, axis=0)

    idxs_min_check = np.all(locs >= minns, axis=1)
    idxs_max_check = np.all(locs <= maxxs, axis=1)

    idxs = np.argwhere(idxs_min_check*idxs_max_check==True).flatten()

    return idxs




#@njit
def get_cover_of_subset_locs(subset, locs):

    minns = np.min(subset, axis=0)
    maxxs = np.max(subset, axis=0)

    idxs_min_check = np.all(locs >= minns, axis=1)
    idxs_max_check = np.all(locs <= maxxs, axis=1)

    idxs = np.argwhere(idxs_min_check*idxs_max_check==True).flatten()

    return idxs


def get_L0_norm(pers1_diff, pers2_diff):        

    temp1 = pers1_diff.copy()
    temp2 = pers2_diff.copy()
    
    max_len = max(len(temp1), len(temp2))
    temp1.resize(max_len)
    temp2.resize(max_len)

    L0_norm = np.max(np.abs(temp1 - temp2))

    return L0_norm
    
def get_nearest_neigh(n_pts, locs):

    # Using ball tree to get nearest neighbors
    kd_tree = cKDTree(locs)

    # Get the second nearest neighbor (k = 2)
    # because the pt itself is the nearest one
    nearest, dummy = kd_tree.query(locs, k=[2])

    nearest = nearest.flatten()

    return nearest



def get_PH_perturbs(original_locs, max_sig, n_samples, pert_dirr, filetype, thresh, threads, dim):

    print('Getting perturbations...')

    # Divide in ten intervals of halving
    init = max_sig
    n_d_min = 10
    d_mins = []

    birth_thresh = thresh

    for i in range(n_d_min):
        init = init/2
        d_mins.append(init)
    
        d_min_dir = pert_dirr+'/dmin_'+str(i)
        if not os.path.isdir(d_min_dir): 
            os.mkdir(d_min_dir)
    
    
    d_mins_file = pert_dirr+'/d_mins.txt'
    np.savetxt(d_mins_file, np.array(d_mins))
    
    
    # Get the distance of closest neighbor for every point 
    n_pts, embed_dim = original_locs.shape
    
    nearest_file = pert_dirr+'/nearest.npy'
    
    if os.path.isfile(nearest_file):
        # load if nearest is already computed
        print('nearest is already computed')
        nearest = np.load(nearest_file)
    else:
        # compute nearest 
        nearest = np.ones((n_pts,)) * np.inf
        print('computing nearest for num points', n_pts)
        nearest = get_nearest_neigh(n_pts, original_locs, nearest)
        np.save(nearest_file, nearest)
    

    print('min and max in nearest is', np.amin(nearest), np.amax(nearest))

    #print(nearest.shape)
    #plt.cla()
    #plt.clf()
    #plt.hist(np.log10(nearest), bins=100)
    #plt.show()
    #exit()
    
    
    # Get Perturbations
    for i, d_min in enumerate(d_mins):

    
        print(i, 'd_min is', d_min)
    
        d_min_dir = pert_dirr+'/dmin_'+str(i)
    
        if not os.path.isdir(d_min_dir): 
            print('Something went wrong? Quitting')
            exit()
    
        temp_nearest = nearest.copy()

    
        temp_nearest[temp_nearest > d_min] = d_min
    
        for it in range(n_samples):
    
    
            this_perturb_dir = d_min_dir+'/perturb_'+str(it)
            if not os.path.isdir(this_perturb_dir): 
            	os.mkdir(this_perturb_dir)
    
            print('Processing: perturbing number', it, ' out of', n_samples\
                    , 'perturbation idx', i, 'out of', len(d_mins), end='\r')
    
            rand_rad = np.random.random(n_pts)*temp_nearest/3
    
            if embed_dim == 2:
                rand_ang = np.random.random(n_pts)*math.pi*2
                new_locs_x = original_locs[:,0] + rand_rad*np.cos(rand_ang)
                new_locs_y = original_locs[:,1] + rand_rad*np.sin(rand_ang)
                new_locs = np.column_stack((new_locs_x, new_locs_y))
            elif embed_dim == 3:
                rand_theta = np.random.random(n_pts)*math.pi*2
                rand_phi = np.random.random(n_pts)*math.pi
                new_locs_x = original_locs[:,0] + rand_rad*np.cos(rand_theta)*np.sin(rand_phi)
                new_locs_y = original_locs[:,1] + rand_rad*np.sin(rand_theta)*np.sin(rand_phi)
                new_locs_z = original_locs[:,2] + rand_rad*np.cos(rand_phi)
                new_locs = np.column_stack((new_locs_x, new_locs_y, new_locs_z))
    
    
            new_source = this_perturb_dir+'/locs.csv'
    
            np.savetxt(new_source, new_locs, delimiter=',', fmt='%.4f')
    
            temp_target = this_perturb_dir+'/result'
    
            import pydory as dory
            dory.compute_PH(new_source, thresh, filetype, threads, temp_target, dim, 1, 1,
                    thresh, 1)
            del dory
            gc.collect()
    
    print('Done getting perturbations.')
         
    
def get_L0_max_cover_all_perturbs(pert_dirr, analyze_dim, n_dmins, n_samples, original_locs):

    print('Getting L0 and max cover for all perturbs...')

    original_pers = np.loadtxt(pert_dirr+'/originalH'+str(analyze_dim)+'_pers_data.txt', delimiter=',')
    cyc_file = open(pert_dirr+'/originalminimal_V_birth_H'+str(analyze_dim)+'.txt'\
                        , 'r')
    original_max_cover = 0
    for cyc in cyc_file:
        cyc = cyc.split(',')
        cyc = cyc[:-1]
        cyc = [int(x) for x in cyc]

        idxs = get_cover_of_cycle(cyc, original_locs)
        if len(idxs) > original_max_cover:
            original_max_cover = len(idxs)


    pers_diff = original_pers[:,1] - original_pers[:,0]

    pers_diff = -np.sort(-pers_diff)

    d_mins_file = pert_dirr+'/d_mins.txt'
    d_mins = np.loadtxt(d_mins_file, delimiter=',')

    L0_medians = np.zeros((n_dmins,))
    min_max_covers = np.zeros((n_dmins,))
    n_common_features = np.zeros((n_dmins,))

    sig_thresh = np.zeros((n_dmins,))
    labels = []




    # Get Perturbations
    for i in range(n_dmins):

        print('d_min_idx', i)

        d_min_dir = pert_dirr+'/dmin_'+str(i)

        L0_norms = np.zeros((n_samples,))
        max_covers = np.zeros((n_samples,))

        if not os.path.isdir(d_min_dir): 
            print('Something went wrong? Quitting')
            exit()
    
        max_n_features = 0
        all_pers_diff = []

        for it in range(n_samples):

            this_perturb_dir = d_min_dir+'/perturb_'+str(it)
            if not os.path.isdir(this_perturb_dir): 
                print('Something went wrong? Quitting')
                exit()

            temp_target = this_perturb_dir+'/result'

            pers = np.loadtxt(temp_target+'H'+str(analyze_dim)+'_pers_data.txt', delimiter=',')

            this_pers_diff = pers[:,1] - pers[:,0]

            this_pers_diff = -np.sort(-this_pers_diff)

            # This will be used in std dev in number of features
            all_pers_diff.append(this_pers_diff.copy())
            if len(this_pers_diff) > max_n_features:
                max_n_features = len(this_pers_diff)

            # Get L0 norm
            if (len(pers_diff) > len(this_pers_diff)):
                maxlen = len(pers_diff)
                this_pers_diff.resize(maxlen)
            elif (len(pers_diff) < len(this_pers_diff)):
                maxlen = len(this_pers_diff)
                pers_diff.resize(maxlen)

            L0_norm = np.max(np.abs(pers_diff - this_pers_diff))

            L0_norms[it] = L0_norm

            cyc_file = open(temp_target+'minimal_V_birth_H'+str(analyze_dim)+'.txt'\
                                , 'r')

            locs = np.loadtxt(this_perturb_dir+'/locs.csv', delimiter=',')

            max_cover = 0

            for cyc in cyc_file:
                cyc = cyc.split(',')
                cyc = cyc[:-1]
                cyc = [int(x) for x in cyc]

                idxs = get_cover_of_cycle(cyc, locs)
                if len(idxs) > max_cover:
                    max_cover = len(idxs)


            max_covers[it] = max_cover

            
        #print(L0_norms, max_covers)

        #L0_norms = np.array(L0_norms)
        #max_covers = np.array(max_covers)

        #plt.title(str(d_mins[i]))
        #plt.scatter(L0_norms, max_covers, label=str(round(d_mins[i],3)), alpha=0.5)

        L0_medians[i] = np.median(L0_norms)
        #print(L0_norms, L0_medians[i])

        idxss = np.argwhere(L0_norms <= L0_medians[i]).flatten()
        minn_max_cover = np.amin(max_covers[idxss])


        min_max_covers[i] = np.amin(max_covers[idxss])
        labels.append(str(round(d_mins[i],3)))

        #plt.scatter(L0_norm_median, minn_max_cover, label=str(round(d_mins[i],3)), alpha=0.5)



        ####################################################
        #  Significance PD
        ####################################################
        
        pers_array = np.zeros((max_n_features, n_samples))
        
        
        for ii in range(n_samples):
            diff = all_pers_diff[ii].copy()
            diff.resize(max_n_features)
            pers_array[:,ii] = diff
        
        
        all_pts = []
        n_features = []
        
        for ii in range(max_n_features):
            minn = np.amin(pers_array[ii, :])
            #print(minn)
            lens = []
            for j in range(n_samples):
                #print(pers_array[:,j])
                idxs = np.argwhere(pers_array[:, j] >= minn)
                #print(idxs, len(idxs))
                lens.append(len(idxs))
                #input('w')
        
            lens = np.array(lens)
            #print(np.std(lens))
            #plt.scatter(minn, np.std(lens))
            if minn != 0:
                all_pts.append([minn, np.std(lens)])
        
        all_pts = np.array(all_pts)
        minstd = np.amin(all_pts[:,1])
        minstd_idxs = np.where(all_pts[:,1] == minstd)
        temp_min = all_pts[minstd_idxs]
        
        sig_thresh[i] = np.amin(temp_min[:,0])

        if minstd == 0:
            n_common_features[i] = len(np.argwhere(pers_array[:, 0] >= sig_thresh[i]).flatten())
            print('\n d_min, sig thresh, min std, n_features'\
                    , d_mins[i], sig_thresh[i], minstd, n_common_features[i],  '\n')
        else:
            # This is to say that 1/1 will be highest possible
            n_common_features[i] = 0
            print('\n d_min, sig thresh, min std, n_features'\
                    , d_mins[i], sig_thresh[i], minstd, '0  \n')
        
        #dict_all_info[d_min_idx]['sig_pers_thresh'] = sig_pers
        #dict_all_info[d_min_idx]['min_std_sig_thresh'] = minstd
        
        
        #plt.xlabel('persistence threshold')
        #plt.ylabel('std. dev. of num. of features > thresh')
        #
        #plt.scatter(all_pts[:, 0], all_pts[:, 1])


        #plt.savefig(pictures_dirr+'/significant_pers_thresh_dmin_'+str(d_min_idx)+'_H'+str(analyze_dim)+'.png', format='png')
        
        #plt.show()

        #plt.cla()
        #plt.clf()



    #L0_medians = L0_medians/(np.max(L0_medians))

    #min_max_covers = min_max_covers/(original_max_cover)

    #print(max_n_features)

    #n_common_features = np.log2(max_n_features+1)/np.log2(n_common_features+1)

    #dist = L0_medians**2 + min_max_covers**2 + n_common_features**2

    dist = (L0_medians/np.amax(L0_medians))**2 + (min_max_covers/original_max_cover)**2

    min_idx = np.argmin(dist)

    plt.scatter(L0_medians, min_max_covers,alpha=0.5)

    for pt in range(n_dmins):
        plt.annotate(labels[pt]+', '+str(int(n_common_features[pt])), (L0_medians[pt], min_max_covers[pt]))

    plt.xlabel('median_L0_norm (normed)')
    plt.ylabel('max cover size (normed)')
    plt.savefig(pert_dirr+'/robust_significance_analysis.pdf', format='pdf')
    #plt.legend()
    #plt.show()

    print('Done getting L0 and max cover for all perturbs.')

    return min_idx, L0_medians[min_idx], sig_thresh[min_idx]






def get_unique_covers(pert_dirr, d_min_idx, L0_norm, analyze_dim, n_samples, original_locs):

    print('Computing unique covers from the perturbations for d_min_idx', d_min_idx)


    d_min_dir = pert_dirr+'/dmin_'+str(d_min_idx)

    original_pers = np.loadtxt(pert_dirr+'/originalH'+str(analyze_dim)+'_pers_data.txt', delimiter=',')
    original_diff = original_pers[:,1] - original_pers[:,0]
    original_diff = -np.sort(-original_diff)

    unique_covers = []

    #test_covers = frozenset({})

    for it in range(n_samples):

        print(it, n_samples, end='\r')


        this_perturb_dir = d_min_dir+'/perturb_'+str(it)
        if not os.path.isdir(this_perturb_dir): 
            print('Something went wrong? Quitting')
            exit()

        temp_target = this_perturb_dir+'/result'

        pers = np.loadtxt(temp_target+'H'+str(analyze_dim)+'_pers_data.txt', delimiter=',')

        this_pers_diff = pers[:,1] - pers[:,0]

        this_pers_diff = -np.sort(-this_pers_diff)

        pers_diff = original_diff.copy()

        # Get L0 norm
        if (len(pers_diff) > len(this_pers_diff)):
            maxlen = len(pers_diff)
            this_pers_diff.resize(maxlen)
        elif (len(pers_diff) < len(this_pers_diff)):
            maxlen = len(this_pers_diff)
            pers_diff.resize(maxlen)

        this_L0_norm = np.max(np.abs(pers_diff - this_pers_diff))

        # If L0 norm is greater than threshold, then skip this perturbation
        if (this_L0_norm > L0_norm):
            continue


        # Otherwise, check for unique covers

        
        cyc_file = open(temp_target+'minimal_V_birth_H'+str(analyze_dim)+'.txt'\
                            , 'r')

        locs = np.loadtxt(this_perturb_dir+'/locs.csv', delimiter=',')

        for cyc in cyc_file:
            cyc = cyc.split(',')
            cyc = cyc[:-1]
            cyc = [int(x) for x in cyc]

            idxs = get_cover_of_cycle(cyc, locs)

            if frozenset(idxs) not in unique_covers:
                unique_covers.append(frozenset(idxs))

            #test_covers = test_covers.union(frozenset(idxs))


    #plt.cla()
    #plt.clf()
    #plt.scatter(original_locs[:,0], original_locs[:,1], marker='x', alpha=0.5, color='blue')
    #plt.scatter(original_locs[list(test_covers),0], original_locs[list(test_covers),1]\
    #                    , marker='x', alpha=0.5, color='red')
    #
    #plt.show()

    print('Number of unique covers', len(unique_covers))

    return unique_covers


def graph_analysis(unique_covers, sig_thresh, original_locs\
            , filetype, thresh, threads, dim, analyze_dim, pert_dirr):


    original_pers = np.loadtxt(pert_dirr+'/originalH'+str(analyze_dim)+'_pers_data.txt', delimiter=',')
    original_diff = original_pers[:,1] - original_pers[:,0]
    original_diff = original_diff[original_diff >= sig_thresh]
    print('sig in original', original_diff)

    birth_thresh = thresh
    

    # Find PH of covers w.r.t. original locations, if not already found
    cover_info_path = pert_dirr+'/cover_info_dict.p'
    n_sig_dict_path = pert_dirr+'/n_sig_dict.p'

    plt.cla()
    plt.clf()
    plt.scatter(original_locs[:,0], original_locs[:,1]\
                    , alpha=0.5, marker='x', color='blue')
    
    if not os.path.isfile(cover_info_path):
        print('Computing PH of unique covers for the first time')
        # We note the maximum barcode in each cover, and check that it is more than the sig_thresh

        cover_info_dict = {}
        n_sig_dict = {}

        for idx, cover in enumerate(unique_covers):

            cover_locs = original_locs[list(cover)]
            new_source = 'scratch.csv'
            np.savetxt(new_source, cover_locs, delimiter=',', fmt='%.4f')

            temp_target = 'covers'

            import pydory as dory
            dory.compute_PH(new_source, 1000000, filetype, threads, temp_target, dim, 0, 0,
                    birth_thresh, 1)
            del dory
            gc.collect()

            pers = np.loadtxt(temp_target+'H'+str(analyze_dim)+'_pers_data.txt', delimiter=',')

            if len(pers) == 0:
                continue
            if pers.ndim == 1:
                diff_pers = pers[1] - pers[0]
                if diff_pers < sig_thresh-5:
                    continue
                n_sig = 1
            else:
                diff_pers = pers[:,1] - pers[:,0]
                diff_pers = diff_pers[diff_pers >= sig_thresh-5]
                n_sig = diff_pers.shape[0]
                if n_sig == 0:
                    continue


            if n_sig not in n_sig_dict:
                n_sig_dict[n_sig] = []

	    	
            n_sig_dict[n_sig].append(idx)
            
            cover_info_dict[idx] = {}
            cover_info_dict[idx]['cover'] = cover
            cover_info_dict[idx]['num_sig'] = n_sig
            
            pickle.dump( cover_info_dict, open( cover_info_path, "wb" ) )
            pickle.dump( n_sig_dict, open( n_sig_dict_path, "wb" ) )

    else:
        print('Loading pre-computed PH of covers')

    cover_info_dict = pickle.load( open( cover_info_path, "rb" ) )
    n_sig_dict = pickle.load( open( n_sig_dict_path, "rb" ) )


    all_test_cover = frozenset({})
    for c_idx in cover_info_dict:
        cover = cover_info_dict[c_idx]['cover']
        all_test_cover = all_test_cover.union(cover)
    
    plt.scatter(original_locs[list(all_test_cover), 0], original_locs[list(all_test_cover), 1]\
                , alpha=0.5, marker='x', color='red')


    cover_source = pert_dirr+'/scratch.csv'
    np.savetxt(cover_source, original_locs[list(all_test_cover)], delimiter=',', fmt='%.4f')
    cover_target = pert_dirr+'/scratch'
    
    import pydory as dory
    dory.compute_PH(cover_source, 100000, filetype, threads, cover_target, dim, 0, 0, birth_thresh, 1)
    del dory
    gc.collect()


    pers = np.loadtxt(cover_target+'H'+str(analyze_dim)+'_pers_data.txt', delimiter=',')
    diff_pers = pers[:,1] - pers[:,0]
    diff_pers = diff_pers[diff_pers >= sig_thresh-5]
    print(diff_pers)

    plt.show()

    exit()

    # At this point, we have all unique covers with number of sig features.
    # Time to make the subset graph and reject the ones that are subsets of each other and have same
    # number of significant features

    G = nx.DiGraph()

    max_idx = 0
    # Add nodes
    for idx in cover_info_dict:
        G.add_node(idx, num_sig=cover_info_dict[idx]['num_sig'])
        if idx > max_idx:
            max_idx = idx

    H = copy.deepcopy(G)

    print('nodes before', len(H.nodes))

    for n_sig in n_sig_dict:

        # subset check
        for ci_idx, cj_idx in itertools.combinations(list(n_sig_dict[n_sig]), 2):
        
        
        	if ci_idx not in list(H.nodes) or cj_idx not in list(H.nodes):
        		continue
        	
        	ci = cover_info_dict[ci_idx]['cover'] 
        	cj = cover_info_dict[cj_idx]['cover'] 
        	
        	if ci.issuperset(cj):
        		H.remove_node(ci_idx)
        
        	if cj.issuperset(ci):
        		H.remove_node(cj_idx)

			
    print('nodes after', len(H.nodes))


    # subset graph
    for ci_idx, cj_idx in itertools.combinations(list(H.nodes), 2):
    
        ci = cover_info_dict[ci_idx]['cover']
        cj = cover_info_dict[cj_idx]['cover']

        if ci.issuperset(cj):
        	H.add_edge(ci_idx, cj_idx)
        
        if cj.issuperset(ci):
        	H.add_edge(cj_idx, ci_idx)


    # intersection graph
    intersectG = nx.Graph()
    intersectG.add_nodes_from(H.nodes(data=True))


    for (ci_idx, cj_idx) in itertools.combinations(list(intersectG.nodes), 2):

        ci = cover_info_dict[ci_idx]['cover'] 
        cj = cover_info_dict[cj_idx]['cover'] 
        
        if (not ci.issuperset(cj)) and (not cj.issuperset(ci)) and (not ci.isdisjoint(cj)):
                intersectG.add_edge(ci_idx, cj_idx, process=0)
                print('adding edge', ci_idx, cj_idx)
    
    update = 1


    while (update):

        update = 0
        
        print('NEW ITERATION')
        temp_data = []
        for ci_idx, cj_idx, attr in intersectG.edges(data=True):
        	temp_data.append((ci_idx, cj_idx, attr['process']))
        
        
        for ci_idx, cj_idx, process in temp_data:
        
                #print(process)
                #print('checking if edge to be processed')
                if process == 1:
                	continue
                
                ci = cover_info_dict[ci_idx]['cover'] 
                cj = cover_info_dict[cj_idx]['cover'] 

                print('checking', ci_idx, cj_idx)
                
                #print('process edge')
                ck = ci.intersection(cj)

                
                #print('ci len, cj len, ck len', len(ci), len(cj), len(ck))
                
                intersectG[ci_idx][cj_idx]['process'] = 1
                
                if (len(ck) < 3):
                	continue
                
                cover_locs = original_locs[list(ck)]
                
                cover_source = pert_dirr+'/scratch.csv'
                
                np.savetxt(cover_source, cover_locs, delimiter=',', fmt='%.4f')
                
                cover_target = pert_dirr+'/temp_cover'
                
                #print('Computing PH...')
                import pydory as dory
                dory.compute_PH(cover_source, 100000, filetype, threads, cover_target, dim, 0, 0,\
                        birth_thresh, 1)
                del dory
                gc.collect()
                #print('Computed PH.')
                
                pers = np.loadtxt(cover_target+'H'+str(analyze_dim)+'_pers_data.txt', delimiter=',')
                
                if (len(pers) == 0):
                	continue
                
                if (pers.ndim == 1):
                	pers = pers.reshape((1,2))
                
                
                #pers = pers[pers[:,1] != -1]
                
                diff_p = pers[:,1] - pers[:,0]
                
                diff_p = diff_p[diff_p >= sig_thresh]
                
                num_sig = len(diff_p)
                
                # DOUUUUUUUUUUUUUUUBBBBBBBBBBBBBBTTTTTTT
                # If the intersection has no significant features, then remove the edge (ci_idx, cj_idx)
                if num_sig == 0:
                	intersectG.remove_edge(ci_idx, cj_idx)
                	continue
                
                # If ck does not have same number of sig features as both ci and cj
                if num_sig != intersectG.nodes[ci_idx]['num_sig'] and\
                        num_sig != intersectG.nodes[cj_idx]['num_sig']:
                	continue
                
                
                # At this point, ck has same number of sig features as either ci, cj or both

                
                # Add the intersection as a node
                max_idx += 1
                ck_idx = max_idx

                print('Adding node', ck_idx, 'intersection of', ci_idx, cj_idx)

                cover_info_dict[ck_idx] = {}
                cover_info_dict[ck_idx]['cover'] = ck
                cover_info_dict[ck_idx]['num_sig'] = num_sig
                
                intersectG.add_node(ck_idx, num_sig=num_sig)
                
                if num_sig == intersectG.nodes[ci_idx]['num_sig']:
                	print('removing node', ci_idx)
                	intersectG.remove_node(ci_idx)
                
                if num_sig == intersectG.nodes[cj_idx]['num_sig']:
                	print('removing node', cj_idx)
                	intersectG.remove_node(cj_idx)
                
                temp_nodes = list(intersectG.nodes)[:]
                
                # Subset check again because new node has been added
                for node in temp_nodes:
                        if node == ck_idx:
                        	continue
                        node_cover = cover_info_dict[node]['cover']
                        if ck.issubset(node_cover) and num_sig == intersectG.nodes[node]['num_sig']:
                        	print('removing node')
                        	intersectG.remove_node(node)
                
                
                for node in list(intersectG.nodes):
                        if node == ck_idx:
                        	continue
                        node_cover = cover_info_dict[node]['cover']
                        ckk = ck.intersection(node_cover)
                        
                        if (len(ckk) != 0) and (ckk != ck) and (ckk != node):
                        	intersectG.add_edge(ck_idx, node,  process = 0)
                
                update = 1
                break
                
                
                intersectG[ci_idx][cj_idx]['process'] = 1
                        		
            
    exit()
        
    # subset graph
    H = nx.DiGraph()
    H.add_nodes_from(intersectG.nodes(data=True))
    
    unG = nx.Graph()
    unG.add_nodes_from(intersectG.nodes(data=True))
    
    for (ci_idx, cj_idx) in itertools.combinations(list(H.nodes), 2):
        ci = cover_info_dict[ci_idx]['cover']
        cj = cover_info_dict[cj_idx]['cover']
        if ci.issuperset(cj):
            H.add_edge(ci_idx, cj_idx)
            unG.add_edge(ci_idx, cj_idx)
        elif cj.issuperset(ci):
            H.add_edge(cj_idx, ci_idx)
            unG.add_edge(cj_idx, ci_idx)
    
    
    nx.write_gpickle(H, pert_dirr+"/subset_graph.gpickle")
    
    
    components = [c for c in sorted(nx.connected_components(unG), key=len, reverse=False)]
    
    for ii, comp in enumerate(components):
    
        dict_node = {}
    
        print('component', ii)
    
        comp_graph = H.subgraph(list(comp))
    
        new_g = nx.DiGraph()
    
        count = 0
        for node in list(comp_graph.nodes):
            dict_node[node] = count
            new_g.add_node(count, num_sig=comp_graph.nodes[node]['num_sig'])
            count += 1
    
        
    
        for ci_idx, cj_idx in comp_graph.edges:
            new_g.add_edge(dict_node[ci_idx], dict_node[cj_idx])
    
    
        
    
        labels = nx.get_node_attributes(new_g, 'num_sig')
    
        if len(comp) > 1:
            pos = graphviz_layout(new_g, prog="dot")
            nx.draw(new_g,pos=pos,labels=labels,node_size=1000)
        else:
            nx.draw(new_g,labels=labels,node_size=1000)
    
        #plt.show()
        plt.savefig(pert_dirr+'/subset_component_'+str(ii)+'.png', format='png')
        plt.cla()
        plt.clf()
    
    
        #plt.scatter(locs[:, 0], locs[:, 1], alpha=0.1)
    
    
        plt.scatter(original_locs[:, 0], original_locs[:, 1], alpha=0.5, marker='+')

        for c in comp:

            #print(c)

            #print(cover_info_dict[c]['cover'])
    
            #if (H.out_degree(c) != 0):
            #    continue

                #leaf_nodes.append(c)
    
            #plt.scatter(locs[:, 0], locs[:, 1], alpha=0.1)
    
            print(intersectG.nodes[c]['num_sig'], H.out_degree(c))
    
            #new_source = target+'_test_cover.csv'
            #new_locs = locs[list(c)]
    
            #np.savetxt(new_source, new_locs, delimiter=',')
            #new_target = target + '_test_cover'
    
            #pydory.compute_PH(new_source, 10000, filetype, threads, new_target, dim , 1, 1,\
            #birth_thresh, 1)
    
            #pers = np.loadtxt(new_target+'H1_pers_data.txt', delimiter=',')
            #diff = pers[:,1] - pers[:,0]
            #plt.scatter(np.zeros((len(diff),)), diff)
            #plt.show()
            #exit()
    
    
            
    
            plt.scatter(original_locs[list(cover_info_dict[c]['cover']), 0]\
                                    , original_locs[list(cover_info_dict[c]['cover']), 1], alpha=0.5, marker='x')
    
        #plt.show()
        plt.pause(3)
        plt.cla()
        plt.clf()
    
    
            
            
                
            
            
def match_points_2D(source, target):

    pair_dist = sklearn.metrics.pairwise_distances(source, target)

    while pair_dist.shape[0]:
        #idx = np.argmin(pair_dist)
        row, col = np.unravel_index(np.argmin(pair_dist, axis=None), pair_dist.shape)

        print(row, col, pair_dist[row, col])

        pt = source[row]
        nearest_pt = target[col]
        plt.plot([pt[0], nearest_pt[0]], [pt[1], nearest_pt[1]])

        pair_dist = np.delete(pair_dist, row, axis=0)
        pair_dist = np.delete(pair_dist, col, axis=1)

        source = np.delete(source, row, axis=0)
        target = np.delete(target, col, axis=0)



    #idxs = np.argsort(pair_dist)
    #print(idxs.shape)
    #print(idxs)
    #print(idxs[0,0])

    #print(pair_dist.shape)

    #while(len(source)):

    #    d_min = math.inf
    #    n_rows = target.shape[0]
    #    for pt in source:
    #        nearest_row = -1
    #        for row in range(n_rows):
    #            check_pt = target[row]
    #            dist = np.sum((pt-check_pt)**2) 
    #            if dist < d_min:
    #                d_min = dist
    #                nearest_row = row

    #        if nearest_row == -1:
    #            print('error in nearest row')
    #            exit()

    #        nearest_row = int(nearest_row)
    #        nearest_pt = target[nearest_row]
    #        print(nearest_pt)
    #        plt.plot([pt[0], nearest_pt[0]], [pt[1], nearest_pt[1]])
    #        target = np.delete(target, nearest_row, axis=0)

    #plt.show()

    #exit()
                
            
    #minns = np.min(source, axis=0)
    #maxxs = np.max(source, axis=0)

    #idxs_min_check = np.all(target >= minns, axis=1)
    #idxs_max_check = np.all(target <= maxxs, axis=1)

    #idxs = np.argwhere(idxs_min_check*idxs_max_check==True).flatten()

    #target = target[idxs]

    #plt.cla()
    #plt.clf()

    #plt.scatter(source[:,0], source[:,1], alpha=0.5)
    #plt.scatter(target[:,0], target[:,1], alpha=0.75, marker='+', color='red')

    #plt.show()
    #exit()


def recompute_minimal(pert_dirr, d_min_idx, sig_pers, L0_norm\
                                , analyze_dim, n_samples, original_locs\
                                , filetype, dim, threads):


    print('Finding optimal threshold and then recomputing PH of valid perturbations...')

    


    original_target = pert_dirr+'/original'
    # First, get the max birth of the sig pers in the original data set
    original_pers = np.loadtxt(pert_dirr+'/originalH'+str(analyze_dim)+'_pers_data.txt', delimiter=',')

    original_diff = original_pers[:,1] - original_pers[:,0]

    idxs = np.argwhere(original_diff >= sig_pers).flatten()

    sig_pers_pairs = original_pers[idxs]

    optimal_thresh = np.amax(sig_pers_pairs[:,0])

    print('optimal birth threshold is', optimal_thresh)

    # Plot these features in original data set

    diff = original_pers[:,1] - original_pers[:,0]
    idxs = np.argwhere(diff > sig_pers).flatten()
    original_pers = original_pers[idxs]
    original_diff = diff[idxs]
    

    ## Find number of points in each cover
    cycles = open(original_target+'homH'+str(analyze_dim)+'_cycles.txt', "r")
    
    original_cycle_dict = {}
    original_embed = []
    
    original_cycle_cycle = {}

    all_covers = []
    
    # FOR HOMOLOGY CYCLE
    idx = 0
    skip = 0
    for line in cycles:
        if line[0] != 'b' and line[0] != 'h':
            pers = line
            pers = pers.strip('\n')
            pers = pers.split(',')
            pers = [float(x) for x in pers]
            if pers[1] - pers[0] < sig_pers:
                skip = 1
            else:
                skip = 0
    
        if line[0] != 'h':
            continue
        if skip:
            continue
    
        cyc = line
        cyc_locs = []
        cyc = cyc.strip('\n')
        cyc = cyc.split(',')
        cyc = cyc[1:]
        cyc = [int(x) for x in cyc]
    
        original_cycle_cycle[idx] = cyc
    
        cyc_idxs = np.unique(np.array(cyc))
    
        cover_idxs = get_cover_of_cycle(cyc_idxs, original_locs)
        all_covers.append(cover_idxs)
    
        original_cycle_dict[idx] = cover_idxs
        #original_cycle_dict[idx]['cover'] = cover_idxs
        original_embed.append([pers[1], pers[1] - pers[0]])
        #original_cycle_dict[idx]['diff'] = (pers[1], pers[1] - pers[0])
        idx += 1
    
        #print(pers[0], pers[1])
    
        #all_covers.append(cover_idxs)
        skip = 0
        
    
    #plt.cla()
    #plt.clf()

    original_embed = np.array(original_embed)

    sorted_idxs = np.argsort(-original_embed[:,1])
    

    for key in sorted_idxs:

        #print(original_embed[key, 1])
    
        this_ori_cover = frozenset(original_cycle_dict[key])
        
        ori_hom_cyc = original_cycle_cycle[key]
        
        cyc_idx = 0
        while (cyc_idx < len(ori_hom_cyc)):
            pt1 = ori_hom_cyc[cyc_idx]
            pt2 = ori_hom_cyc[cyc_idx+1]
            plt.plot([original_locs[pt1, 0], original_locs[pt2,0]]\
                    , [original_locs[pt1, 1], original_locs[pt2,1]]\
                    , color='black'\
                    , zorder = 1000)
            cyc_idx += 2
        
        
        
        #plt.scatter(original_locs[:,0], original_locs[:,1], alpha=0.5)


        #plt.show()
    
    
    # COMPUTE TILL BIRTH = OPTIMAL THRESH

    source = pert_dirr+'/scratch.csv'
    np.savetxt(source, original_locs, delimiter=',', fmt='%.4f')
    temp_target = pert_dirr+'/scratch'

    birth_thresh = optimal_thresh

    import pydory as dory
    dory.compute_PH(source, optimal_thresh, filetype, threads, temp_target, dim, 1, 1, birth_thresh, 1)
    del dory
    gc.collect()

    #plot_minimal_cycles(original_locs, temp_target, analyze_dim)
    #plt.show()
    #exit()

    # FOR MINIMAL CYCLE
    cycles = open(temp_target+'minimal_V_birth_H'+str(analyze_dim)+'.txt', "r")

    all_minimal_covers = []
    all_minimal_cycles = []
    
    max_cover_size = 0
    
    for cyc in cycles:
        cyc_locs = []
        cyc = cyc.split(',')
        cyc = cyc[:-1]
        cyc = [int(x) for x in cyc]
    
        all_minimal_cycles.append(cyc)
    
        cyc_idxs = np.unique(np.array(cyc))
    
    
        cover_idxs = get_cover_of_cycle(cyc_idxs, original_locs)
        all_minimal_covers.append(cover_idxs)
    
        if len(cover_idxs) > max_cover_size:
            max_cover_size = len(cover_idxs)
    
    
    all_pers = []
    minimal_cyc_idx = {}
    minimal_cyc_cyc_idx = {}
    minimal_embed = np.zeros((len(all_minimal_covers),2))
    max_sig = 0
    
    # If not, then do the following:
    for idx, cover in enumerate(all_minimal_covers):
        print('Computing PH of cover', idx, 'out of', len(all_minimal_cycles), end='\r')
    
        cover_locs = original_locs[cover]
        np.savetxt('scratch.csv', cover_locs, delimiter=',', fmt='%.4f')
    
        temp_source = 'scratch.csv'
        temp_target = pert_dirr+'/scratch'
    
        import pydory as dory
        dory.compute_PH(temp_source, 1000000, filetype, threads, temp_target, dim, 0, 0,\
                1000000, 1)
        del dory
        gc.collect()
    
        pers = np.loadtxt(temp_target+'H'+str(analyze_dim)+'_pers_data.txt', delimiter=',')
    
    
        if len(pers) == 0:
            continue
        if pers.ndim == 1:
            pers = pers.reshape((1,2))
            
        all_pers += list(pers)
    
        diff = pers[:,1] - pers[:,0]
        sorted_idxs = np.argsort(-diff)
    
    
        minimal_cyc_idx[idx] = cover
    
        minimal_cyc_cyc_idx[idx] = all_minimal_cycles[idx]
    
        #minimal_cyc_idx[idx]['cover'] = cover
        #minimal_cyc_idx[idx]['diff'] = (pers[sorted_idxs[0], 1], pers[sorted_idxs[0],1] - pers[sorted_idxs[0],0])
        minimal_embed[idx, 0] = pers[sorted_idxs[0], 1]
        minimal_embed[idx, 1] = pers[sorted_idxs[0],1] - pers[sorted_idxs[0],0]
    
        #plt.scatter(pers[:,0], pers[:,1], color='red', alpha=0.85, marker='+')
        #plt.scatter(pers[:,1], pers[:,1] - pers[:,0], color='red', alpha=0.85, marker='+')
    
        if (np.amax(diff) > max_sig):
            max_sig = np.amax(diff)
    
        
        
    processed = np.zeros((len(all_minimal_covers),))
    
    
    for key in original_cycle_dict:
    
        this_ori_cover = frozenset(original_cycle_dict[key])
    
        ori_hom_cyc = original_cycle_cycle[key]
    
        cyc_idx = 0
        while (cyc_idx < len(ori_hom_cyc)):
            pt1 = ori_hom_cyc[cyc_idx]
            pt2 = ori_hom_cyc[cyc_idx+1]
            plt.plot([original_locs[pt1, 0], original_locs[pt2,0]]\
                    , [original_locs[pt1, 1], original_locs[pt2,1]]\
                    , color='black'\
                    , zorder = 1000)
            cyc_idx += 2
    
    
    
        ori_embed = original_embed[key]
        min_dist = math.inf
        min_key2 = 0
    
        plt.scatter(original_locs[:,0], original_locs[:,1], alpha=0.5)
    
        #plt.scatter(original_locs[list(this_ori_cover), 0], original_locs[list(this_ori_cover), 1]\
        #            , alpha=0.75, marker='x', color='green')
    
    
        for key2 in minimal_cyc_idx:
            if processed[key2]:
                continue
    
            this_cover = frozenset(minimal_cyc_idx[key2])
            if (len(this_cover.intersection(this_ori_cover))/len(this_cover) < 0.5)\
                    and (len(this_ori_cover.intersection(this_cover))/len(this_ori_cover) < 0.5):
                continue
    
    
    
            #plt.scatter(original_locs[list(this_cover), 0], original_locs[list(this_cover), 1]\
            #        , alpha=0.75, marker='o', color='orange')
    
            min_embed = minimal_embed[key2]
            dist = np.sum((ori_embed - min_embed)**2)
    
            if dist < min_dist:
                min_dist = dist
                min_key2 = key2
    
        print(min_dist)
    
        processed[min_key2] = 1
    
    
        this_min_cover = minimal_cyc_idx[min_key2]
    
        this_hom_cyc = minimal_cyc_cyc_idx[min_key2]
    
        cyc_idx = 0
        while (cyc_idx < len(this_hom_cyc)):
            pt1 = this_hom_cyc[cyc_idx]
            pt2 = this_hom_cyc[cyc_idx+1]
            plt.plot([original_locs[pt1, 0], original_locs[pt2,0]]\
                    , [original_locs[pt1, 1], original_locs[pt2,1]]\
                    , color='red'\
                    , zorder = 1000)
            cyc_idx += 2
    
        #plt.scatter(original_locs[list(this_min_cover), 0], original_locs[list(this_min_cover), 1]\
        #            , alpha=0.75, marker='+', color='red')
        
    plt.show()
        
            
        
        
        
def plot_minimal_cycles(locs, target, analyze_dim):        


    # FOR MINIMAL CYCLE
    cycles = open(target+'minimal_V_birth_H'+str(analyze_dim)+'.txt', "r")

    for cyc in cycles:
        cyc_locs = []
        cyc = cyc.split(',')
        cyc = cyc[:-1]
        cyc = [int(x) for x in cyc]
        cyc_idx = 0
        while (cyc_idx < len(cyc)):
            pt1 = cyc[cyc_idx]
            pt2 = cyc[cyc_idx+1]
            plt.plot([locs[pt1, 0], locs[pt2,0]]\
                    , [locs[pt1, 1], locs[pt2,1]]\
                    , color='red'\
                    , zorder = 1000)
            cyc_idx += 2
        
        
def statistical_analysis(pert_dirr, original_locs, d_min_idx, L0_thresh\
                        , filetype, dim, analyze_dim, threads, n_samples):


    print('Statistical analysis of perturbations...')

    d_min_dir = pert_dirr+'/dmin_'+str(d_min_idx)
    
    original_pers = np.loadtxt(pert_dirr+'/originalH'+str(analyze_dim)+'_pers_data.txt', delimiter=',')

    original_diff = original_pers[:,1] - original_pers[:,0]

    original_diff = -np.sort(-original_diff)

    unique_minimal_covers = []

    for pert in range(n_samples):

        print('Checking L0 dist for pert', pert, 'out of', n_samples, end='\r')

        # Step 1: First remove perturbations that are farther than L0 norm

        this_perturb_dir = d_min_dir+'/perturb_'+str(pert)
        if not os.path.isdir(this_perturb_dir): 
            print('Something went wrong? Quitting')
            exit()
        
        temp_target = this_perturb_dir+'/result'

        pers = np.loadtxt(temp_target+'H'+str(analyze_dim)+'_pers_data.txt', delimiter=',')

        this_pers_diff = pers[:,1] - pers[:,0]

        this_pers_diff = -np.sort(-this_pers_diff)

        L0_norm = get_L0_norm(original_diff, this_pers_diff)

        if (L0_norm > L0_thresh):
            continue

        
        # Go over all minimal cycles, compute covers, and check for uniqueness


        # FOR MINIMAL CYCLE
        cycles = open(temp_target+'minimal_V_birth_H'+str(analyze_dim)+'.txt', "r")

        for cyc in cycles:

            cyc_locs = []
            cyc = cyc.split(',')
            cyc = cyc[:-1]
            cyc = [int(x) for x in cyc]
            cyc = list(frozenset(cyc))

            cover_idxs = frozenset(get_cover_of_cycle(cyc, original_locs))
            if len(cover_idxs) > 2450:
                continue
            if cover_idxs not in unique_minimal_covers:
                unique_minimal_covers.append(cover_idxs)


    print('Done with getting unique minimal covers. Computing PH and distr. of significance...')
            
    #point_stats = {}
    #all_barcodes = []
    n_pts = original_locs.shape[0]

    all_freq = np.zeros((n_pts,))
        
    for cover in unique_minimal_covers:

        locs = original_locs[list(cover)]
        source = pert_dirr+'/scratch.csv'
        np.savetxt(source, locs, delimiter=',', fmt='%.4f')
        target = pert_dirr+'/scratch'


        # FOR NOW ASSUMING THAT THIS CAN BE DONE FOR ALL COVERS TILL INIFINITTY
        import pydory as dory
        dory.compute_PH(source, 0, -1, filetype, threads, target, dim, 1, 1, -1, 1)
        del dory
        gc.collect()

        pers = np.loadtxt(target+'H'+str(analyze_dim)+'_pers_data.txt', delimiter=',')

        if len(pers) == 0:
            continue
        if pers.ndim == 1:
            pers = pers.reshape((1,2))

        diff_p = pers[:,1] - pers[:,0]
        #print(diff_p)

        if np.amax(diff_p) > 0.07:
            #plt.scatter(original_locs[:,0], original_locs[:,1], alpha=0.1, color='blue')
            #plt.scatter(original_locs[list(cover),0], original_locs[list(cover),1]\
            #            , alpha=0.4, color='red')
            all_freq[list(cover)] += 1
            #plt.show()

        #diff_p = list(diff_p)

        #all_barcodes += diff_p

        #for pt in cover:
        #    if pt not in point_stats:
        #        point_stats[pt] = diff_p
        #    else:
        #        point_stats[pt] += diff_p


    #return point_stats, all_barcodes
    return all_freq


def get_varying_thresh(pert_dirr, original_source, low_thresh, up_thresh\
        , original_locs, n_grid\
        , filetype, dim, threads, analyze_dim\
        , sig_thresh):
    

    print('Doing varying thresh analysis...')
    all_thresh = np.linspace(low_thresh, up_thresh, n_grid)

    markers = [\
              '$a$'\
            , '$b$'\
            , '$c$'\
            , '$d$'\
            , '$e$'\
            , '$f$'\
            , '$g$'\
            , '$h$'\
            , '$i$'\
            , '$j$'\
            , '$k$'\
            , '$l$'\
            , '$m$'\
            , '$n$'\
            , '$o$'\
            ]
    
    counter = 0

    #plt.scatter(original_locs[:,0], original_locs[:,1], alpha=0.5)

    for idx, thresh in enumerate(all_thresh):

        print('Processing thresh', thresh, ', number', idx, 'out of', n_grid)

        temp_target = pert_dirr+'/scratch'



        import pydory as dory
        dory.compute_PH(original_source, 0, thresh, filetype, threads, temp_target, dim, 1, 1,\
                thresh, 1)
        del dory
        gc.collect()
        
        this_pers = np.loadtxt(temp_target+'H'+str(analyze_dim)+'_pers_data.txt', delimiter=',')

        #undead_idxs = np.argwhere(this_pers[:,1] == -1).flatten()

        #this_pers[undead_idxs, 1] = thresh

        #diff_pers = this_pers[:,1] - this_pers[:,0]

        #cond_idxs = np.argwhere(diff_pers > sig_thresh).flatten()

        #this_pers = this_pers[cond_idxs]

        xx = np.ones(this_pers.shape[0])*thresh
        plt.scatter(xx, this_pers[:,0], alpha=0.5)

        #plt.scatter(this_pers[:,0], this_pers[:,1], alpha=0.5)

        #plt.scatter(this_pers[:,1]-this_pers[:,0], this_pers[:,1], alpha=0.5)
        

        min_cyc_file = open(temp_target+'minimal_V_birth_H'+str(analyze_dim)+'.txt'\
                        , 'r')

        pers_this_thresh = []
        

        for cyc in min_cyc_file:
            cyc = cyc.split(',')
            cyc = cyc[:-1]
            cyc = [int(x) for x in cyc]
            
            idxs = get_cover_of_cycle(cyc, original_locs)
            
            cover_locs = original_locs[idxs]
            
            cover_source = pert_dirr+'/scratch_cover_locs.csv'
            
            np.savetxt(cover_source, cover_locs, delimiter=',', fmt = '%.4f')
            
            cover_target = pert_dirr+'/scratch_cover'
            
            import pydory as dory
            dory.compute_PH(cover_source, 0, thresh, filetype, threads, cover_target, dim, 0, 0,\
                    thresh, 1)
            del dory
            gc.collect()
            
            cover_pers = np.loadtxt(cover_target+'H'+str(analyze_dim)+'_pers_data.txt', delimiter=',')
            
            if (len(cover_pers) == 0):
                continue
            
            if (cover_pers.ndim == 1):
                cover_pers = cover_pers.reshape((1,2))
                    
            #undead_idxs = np.argwhere(cover_pers[:,1] == -1).flatten()

            #cover_pers[undead_idxs, 1] = thresh

            #diff_pers = cover_pers[:,1] - cover_pers[:,0]

            #cond_idxs = np.argwhere(diff_pers > sig_thresh).flatten()

            #cover_pers = cover_pers[cond_idxs]

            if len(cover_pers) == 0:
                continue

            if cover_pers.ndim == 1:
                cover_pers.reshape((1,2))


            for per in cover_pers:
                pers_this_thresh.append(per)

            #plt.scatter(original_locs[idxs, 0], original_locs[idxs, 1]\
            #            , alpha=0.85, marker=markers[counter], label=str(round(thresh, 4)))

        if len(pers_this_thresh) == 0:
            print('empty, continuing')
            continue

        pers_this_thresh = np.array(pers_this_thresh)

        if pers_this_thresh.ndim == 1:
            pers_this_thresh.reshape((1,2))

        print(pers_this_thresh)
                    
        #plt.scatter(pers_this_thresh[:,0], pers_this_thresh[:,1]\
        #            , marker='+', alpha=0.5, color='red')

        #plt.scatter(pers_this_thresh[:,1]-pers_this_thresh[:,0], pers_this_thresh[:,1]\
        #            , marker=markers[counter], alpha=0.95, label=str(round(thresh, 4)))

        xx = np.ones(len(pers_this_thresh))*thresh

        #plt.scatter(pers_this_thresh[:,0], pers_this_thresh[:,0]\
        #            , marker=markers[counter], alpha=0.95, label=str(round(thresh, 4)))
        plt.scatter(xx, pers_this_thresh[:,0]\
                    , alpha=0.7, marker = '+')
        counter += 1

        #plt.legend()

        #plt.xlabel('birth')
        #plt.ylabel('death')
        #plt.title('thresh='+str(round(thresh,4)))
        #plt.pause(0.5)
    #plt.legend()
    plt.show()

    print('Done varying thresh analysis.')







def graph_analysis_new(graph, lower_thresh, thresh, epsil, original_locs, scratch_dirr\
                        , analyze_dim, filetype, threads):

    # At this point, we have all unique covers with number of sig features.
    # Time to make the subset graph and reject the ones that are subsets of each other and have same
    # number of significant features

    all_nodes = list(graph.nodes)

    print('Conducting subset check')
    print('Nodes before', len(all_nodes))
    update = 1
    while (update):

        update = 0

        for i, j in itertools.combinations(all_nodes, 2):

            if graph.nodes[i]['n_sig'] != graph.nodes[j]['n_sig']:
                continue

            assert (i != j), 'They should not be equal'

            if i.issuperset(j):
                all_nodes.remove(i)
                update = 1
                break

            if j.issuperset(i):
                all_nodes.remove(j)
                update = 1
                break

    print('Done with subset check')
    print('Nodes after', len(all_nodes))


    H = nx.DiGraph()
    for node in all_nodes:
        H.add_node(node, n_sig = graph.nodes[node]['n_sig'])

    print('Creating intersection graph')
    # intersection graph
    intersectG = nx.Graph()
    intersectG.add_nodes_from(H.nodes(data=True))

    for (ci, cj) in itertools.combinations(list(intersectG.nodes), 2):

        ck = ci.intersection(cj)

        if len(ck) and (ck!=ci) and (ck!=cj):
            # Note: There can be subsets only if they have different number of features
            intersectG.add_edge(ci, cj, process=0)


    #plt.cla()
    #plt.clf()
    #nx.draw(intersectG)
    #plt.show()

    # We can analyze different components of intersectG independently
    intersection_comps = nx.connected_components(intersectG)


    all_graph = nx.Graph()

    #print("There are", len(list(intersection_comps)), 'components')
    #input('w')

    node_sizes_before = []

    for comp_idx, comp in enumerate(intersection_comps):

        print('Doing component', comp_idx, 'of size', len(comp))
        #input('w')

        comp_graph = nx.Graph(intersectG.subgraph(list(comp)))

        update = 1
            
        comp_nodes = list(comp_graph.nodes)

        if len(comp) > 1:
            node_sizes_before.append(len(comp))


        iteration = 0

        #print('before update', nx.is_frozen(comp_graph))

        while (update):

            print('Iteration', iteration, 'nodes in graph', len(comp_nodes))
            iteration += 1

            update = 0
            
            add_edges = []

            remove_edges = []

             
            for ci, cj in list(comp_graph.edges):

                remove_edges.append([ci, cj])

                assert (ci != cj), 'should not be equal'

                ck = ci.intersection(cj)

                if ck==ci or ck==cj:
                    continue

                if len(ck) < 4:
                    continue

                #assert (ck != cj), 'should not be equal'
                #assert (ck != ci), 'should not be equal'

                assert (len(ck) != 0), 'should not be empty intersection'

                #print(len(ci), len(cj), len(ck))

                # Compute PH
                cover_locs = original_locs[list(ck)]
                
                cover_source = scratch_dirr+'/scratch.csv'
                
                np.savetxt(cover_source, cover_locs, delimiter=',', fmt='%.2f')
                
                cover_target = scratch_dirr+'/temp_cover'
                
                #print('Computing PH...')
                # NOTE: birth_thresh here is thresh

                import pydory as dory

                dory.compute_PH(cover_source, lower_thresh, thresh + epsil\
                                , filetype, threads, cover_target, analyze_dim, 0, 0, thresh, 1)

                del dory
                gc.collect()

                #print('Computed PH.')

                
                pers = np.loadtxt(cover_target+'H'+str(analyze_dim)+'_pers_data.txt', delimiter=',')
                
                if (len(pers) == 0):
                	continue
                
                if (pers.ndim == 1):
                	pers = pers.reshape((1,2))
                
                
                pers[pers[:,1] == -1, 1] = thresh
                
                diff_p = pers[:,1] - pers[:,0]
                
                num_sig = len(diff_p[diff_p >= epsil])

                add_flag = 0

                if num_sig == comp_graph.nodes[ci]['n_sig']:
                    comp_nodes.remove(ci)
                    add_flag = 1

                if num_sig == comp_graph.nodes[cj]['n_sig']:
                    comp_nodes.remove(cj)
                    add_flag = 1

                if add_flag:

                    #print('adding and updating')

                    # subset check
                    for cii in comp_graph.nodes:

                        if cii not in comp_nodes:
                            continue

                        if comp_graph.nodes[cii]['n_sig'] == num_sig:
                            if ck.issubset(cii):
                                comp_nodes.remove(cii)
                            
                    # intersection check
                    for cii in comp_nodes:

                        ckk = cii.intersection(ck)

                        if len(ckk):
                            add_edges.append(cii)


                    update = 1
                    break

            if update:
                #print(num_sig)
                #comp_graph = nx.Graph(comp_graph.subgraph(comp_nodes))
                comp_graph = nx.Graph(comp_graph.subgraph(comp_nodes))
                comp_graph.add_node(ck, n_sig = num_sig)
                comp_nodes.append(ck)
                for cii in add_edges:
                    comp_graph.add_edge(cii, ck)
                for cii, cjj in remove_edges:
                    if cii in comp_graph.nodes and cjj in comp_graph.nodes:
                        comp_graph.remove_edge(cii, cjj)


        all_graph = nx.compose(all_graph, comp_graph)
        #all_nodes += list(comp_graph.nodes)


    ## Copy nodes from all_graph
    #final_g = nx.Graph()
    #final_g.add_nodes_from(all_graph)

    #for ci, cj in itertools.combinations(list(final_g), 2):
    #    if len(ci.intersection(cj)):
    #        final_g.add_edge(ci, cj)

    #return list(nx.connected_components(final_g))




    
    return list(all_graph.nodes)








    ###exit()

    ###node_sizes_after = []
    ###for comp in nx.connected_components(final_g):
    ###    if len(comp) == 1:
    ###        continue
    ###    #for nn in comp:
    ###    node_sizes_after.append(len(comp))

    ###plt.gca()
    ###plt.gcf()
    ###plt.gca().boxplot([node_sizes_before, node_sizes_after])
    ###plt.show()
    ####exit()



    ###plt.cla()
    ###plt.clf()
    ###nx.draw(final_g)
    ###plt.show()
    ###            

    ###exit()
    ###    
    #### subset graph
    ###H = nx.DiGraph()
    ###H.add_nodes_from(intersectG.nodes(data=True))
    ###
    ###unG = nx.Graph()
    ###unG.add_nodes_from(intersectG.nodes(data=True))
    ###
    ###for (ci_idx, cj_idx) in itertools.combinations(list(H.nodes), 2):
    ###    ci = cover_info_dict[ci_idx]['cover']
    ###    cj = cover_info_dict[cj_idx]['cover']
    ###    if ci.issuperset(cj):
    ###        H.add_edge(ci_idx, cj_idx)
    ###        unG.add_edge(ci_idx, cj_idx)
    ###    elif cj.issuperset(ci):
    ###        H.add_edge(cj_idx, ci_idx)
    ###        unG.add_edge(cj_idx, ci_idx)
    ###
    ###
    ###nx.write_gpickle(H, pert_dirr+"/subset_graph.gpickle")
    ###
    ###
    ###components = [c for c in sorted(nx.connected_components(unG), key=len, reverse=False)]
    ###
    ###for ii, comp in enumerate(components):
    ###
    ###    dict_node = {}
    ###
    ###    print('component', ii)
    ###
    ###    comp_graph = H.subgraph(list(comp))
    ###
    ###    new_g = nx.DiGraph()
    ###
    ###    count = 0
    ###    for node in list(comp_graph.nodes):
    ###        dict_node[node] = count
    ###        new_g.add_node(count, num_sig=comp_graph.nodes[node]['num_sig'])
    ###        count += 1
    ###
    ###    
    ###
    ###    for ci_idx, cj_idx in comp_graph.edges:
    ###        new_g.add_edge(dict_node[ci_idx], dict_node[cj_idx])
    ###
    ###
    ###    
    ###
    ###    labels = nx.get_node_attributes(new_g, 'num_sig')
    ###
    ###    if len(comp) > 1:
    ###        pos = graphviz_layout(new_g, prog="dot")
    ###        nx.draw(new_g,pos=pos,labels=labels,node_size=1000)
    ###    else:
    ###        nx.draw(new_g,labels=labels,node_size=1000)
    ###
    ###    #plt.show()
    ###    plt.savefig(pert_dirr+'/subset_component_'+str(ii)+'.png', format='png')
    ###    plt.cla()
    ###    plt.clf()
    ###
    ###
    ###    #plt.scatter(locs[:, 0], locs[:, 1], alpha=0.1)
    ###
    ###
    ###    plt.scatter(original_locs[:, 0], original_locs[:, 1], alpha=0.5, marker='+')

    ###    for c in comp:

    ###        #print(c)

    ###        #print(cover_info_dict[c]['cover'])
    ###
    ###        #if (H.out_degree(c) != 0):
    ###        #    continue

    ###            #leaf_nodes.append(c)
    ###
    ###        #plt.scatter(locs[:, 0], locs[:, 1], alpha=0.1)
    ###
    ###        print(intersectG.nodes[c]['num_sig'], H.out_degree(c))
    ###
    ###        #new_source = target+'_test_cover.csv'
    ###        #new_locs = locs[list(c)]
    ###
    ###        #np.savetxt(new_source, new_locs, delimiter=',')
    ###        #new_target = target + '_test_cover'
    ###
    ###        #pydory.compute_PH(new_source, 10000, filetype, threads, new_target, dim , 1, 1,\
    ###        #birth_thresh, 1)
    ###
    ###        #pers = np.loadtxt(new_target+'H1_pers_data.txt', delimiter=',')
    ###        #diff = pers[:,1] - pers[:,0]
    ###        #plt.scatter(np.zeros((len(diff),)), diff)
    ###        #plt.show()
    ###        #exit()
    ###
    ###
    ###        
    ###
    ###        plt.scatter(original_locs[list(cover_info_dict[c]['cover']), 0]\
    ###                                , original_locs[list(cover_info_dict[c]['cover']), 1], alpha=0.5, marker='x')
    ###
    ###    #plt.show()
    ###    plt.pause(3)
    ###    plt.cla()
    ###    plt.clf()
    




def perturbation_analysis(locs, thresh, epsil, analyze_dim\
                        , n_samples\
                        , filetype, threads\
                        , lower_thresh, comp_dirr):



    if not os.path.isdir(comp_dirr) :
        print('Exiting. Perturbation does not exist')
        exit()



    # The idea is: First find the max. persistence
    # when computing PH for thresh + epsil


    # Step 1: Compute max persistence in this
    # WAIT!!!!!! NOOOO, IT DOES NOT MAKE SENSE TO COMPUTE MAX PERSISTENCE

    source = comp_dirr+'/locs.csv'

    target = comp_dirr+'/'

    print('Computing PH...')
    # NOTE: birth_thresh here is thresh
    import pydory as dory
    dory.compute_PH(source, lower_thresh, thresh + epsil, filetype, threads, target, analyze_dim,\
            1, 1, thresh, 1)
    del dory
    gc.collect()
    print('Computed PH')

    # Get max cover size
    cyc_file = open(target+'minimal_V_birth_H'+str(analyze_dim)+'.txt', 'r')
    max_cover_size = 0
    for cyc in cyc_file:
        cyc = cyc.split(',')
        cyc = cyc[:-1]
        cyc = [int(x) for x in cyc]

        idxs = get_cover_of_cycle(cyc, locs)
        if len(idxs) > max_cover_size:
            max_cover_size = len(idxs)

    print('max cover size is', max_cover_size)
    #exit()
    

    # The computationally expensive part: Get the nearest neighbor of every point
    n_pts = locs.shape[0]
    nearest = np.Inf*np.ones(n_pts,)

    print("Getting nearest for ", n_pts, "points...")
    nearest = get_nearest_neigh(n_pts, locs, nearest)

    nearest = np.sqrt(nearest) 

    max_nearest = np.amax(nearest)

    nearest = nearest/2

    print("Done.")


    # 10 perturbation parameters
    delta_list = [max_nearest/2]
    for ii in range(9):
        delta_list.append(delta_list[-1]/2)

    print(delta_list)

    embed_dim = locs.shape[1]

    pert_dirr = comp_dirr+'/Perturbations'
    if not os.path.isdir(pert_dirr) :
        os.mkdir(pert_dirr)

    np.savetxt(pert_dirr+'/delta_list.csv', np.array(delta_list), delimiter=','\
                    , fmt = '%.4f')


    # Create perturbations for every delta

    # Get Perturbations
    for i, delta in enumerate(delta_list):

        delta_dir = pert_dirr+'/delta_'+str(i)
        if not os.path.isdir(delta_dir): 
            os.mkdir(delta_dir)


        temp_nearest = nearest.copy()

        #print(temp_nearest)

        temp_nearest[temp_nearest > delta] = delta

        #print(temp_nearest)
        #exit()

        for it in range(n_samples):

            this_perturb_dir = delta_dir+'/perturb_'+str(it)

            if not os.path.isdir(this_perturb_dir): 
            	os.mkdir(this_perturb_dir)
    
            print('Processing: perturbing number', it, ' out of', n_samples\
                    , 'perturbation idx', i, 'out of', len(delta_list), end='\r')
    
            rand_rad = np.random.random(n_pts)*temp_nearest
    
            if embed_dim == 2:
                rand_ang = np.random.random(n_pts)*math.pi*2
                new_locs_x = locs[:,0] + rand_rad*np.cos(rand_ang)
                new_locs_y = locs[:,1] + rand_rad*np.sin(rand_ang)
                new_locs = np.column_stack((new_locs_x, new_locs_y))
            elif embed_dim == 3:
                rand_theta = np.random.random(n_pts)*math.pi*2
                rand_phi = np.random.random(n_pts)*math.pi
                new_locs_x = locs[:,0] + rand_rad*np.cos(rand_theta)*np.sin(rand_phi)
                new_locs_y = locs[:,1] + rand_rad*np.sin(rand_theta)*np.sin(rand_phi)
                new_locs_z = locs[:,2] + rand_rad*np.cos(rand_phi)
                new_locs = np.column_stack((new_locs_x, new_locs_y, new_locs_z))


            new_source = this_perturb_dir+'/locs.csv'
    
            np.savetxt(new_source, new_locs, delimiter=',', fmt='%.2f')
    
            temp_target = this_perturb_dir+'/'
    
            import pydory as dory
            # NOTE: birth_thresh here is thresh
            dory.compute_PH(new_source, lower_thresh, thresh + epsil, filetype, threads\
                    , temp_target, analyze_dim, 1, 1, thresh, 1)
            del dory
            gc.collect()





    

def get_optimal_delta(dirr, analyze_dim, thresh, epsil, n_samples):


    # Load original pers
    pers = np.loadtxt(dirr+'/H'+str(analyze_dim)+'_pers_data.txt', delimiter=',')
    pers[pers[:, 1] == -1, 1] = thresh + epsil
    original_pers_diff = pers[:, 1] - pers[:, 0]


    pert_dirr = dirr+'/Perturbations'
    if not os.path.isdir(pert_dirr) :
        print("Perturbation dir for this does not exist. Quitting.")
        exit()

    delta_list = np.loadtxt(pert_dirr+'/delta_list.csv', delimiter=',')
         
    locs = np.loadtxt(dirr+'/locs.csv', delimiter=',')


    all_delta_L0_medians = np.zeros(len(delta_list))
    all_delta_max_covers = np.zeros(len(delta_list))

    labels = []

    # Go over each delta
    for i, delta in enumerate(delta_list):

        labels.append(str(round(delta, 4)))

        delta_dir = pert_dirr+'/delta_'+str(i)

        delta_L0_list = np.zeros(n_samples,)
        delta_max_cover_list = np.zeros(n_samples,)

        # Go over each perturbation
        for it in range(n_samples):

            this_perturb_dir = delta_dir+'/perturb_'+str(it)

            temp_target = this_perturb_dir+'/'

            # Load minimal cycles, find covers, and get max cover size
            cyc_file = open(temp_target+'minimal_V_birth_H'+str(analyze_dim)+'.txt', 'r')

            max_cover_size = 0
            for cyc in cyc_file:
                cyc = cyc.split(',')
                cyc = cyc[:-1]
                cyc = [int(x) for x in cyc]

                idxs = get_cover_of_cycle(cyc, locs)
                if len(idxs) > max_cover_size:
                    max_cover_size = len(idxs)

            # Append to the list of delta_
            delta_max_cover_list[it] = max_cover_size

            # Get L0 norm
            temp_target = this_perturb_dir+'/'

            pers = np.loadtxt(temp_target+'H'+str(analyze_dim)+'_pers_data.txt', delimiter=',')

            if pers.ndim == 1:
                pers = pers.reshape((1,2))

            pers[pers[:, 1] == -1, 1] = thresh + epsil

            this_pers_diff = pers[:,1] - pers[:,0]

            L0_norm = get_L0_norm(original_pers_diff, this_pers_diff)

            delta_L0_list[it] = L0_norm



        # Find median of L0
        all_delta_L0_medians[i] = np.median(delta_L0_list)

        shortlisted_perts = np.where(delta_L0_list <=all_delta_L0_medians[i])[0].flatten()

        all_delta_max_covers[i] = np.amax(delta_max_cover_list[shortlisted_perts])



    print(all_delta_L0_medians, all_delta_max_covers)

    plt.cla()
    plt.clf()
        
    plt.scatter(all_delta_L0_medians, all_delta_max_covers)
    for pt in range(len(delta_list)):
        plt.annotate(labels[pt], (all_delta_L0_medians[pt], all_delta_max_covers[pt]))

    plt.show()

            


            
#def plot_hom_3d(source, target, pers_sig):

def get_edges(locs, thresh):
    """Use KDtree or Balltree implementation 
    """

    kd_tree = cKDTree(locs)

    pairs = kd_tree.query_pairs(r=thresh)

    return pairs



def get_edges_from_locs(source, target_dirr, target_file, thresh, round_to):
    """This function makes edge file from locations, with lengths <= thresh
    
    source is the csv file with the locations
    target_dir is the destination directory to save the edge file
    target_file is name of edge file to be stored
    round_to is the decimal precision of edge lengths
    """

    print('Getting edge file from locs', source, target_dirr, target_file, thresh, round_to)

    locs = np.loadtxt(source, delimiter=',')
    n_pts, dim = locs.shape

    pairs = get_edges(locs, thresh)

    print('Got pairs', len(pairs))

    f = open(target_dirr+'/'+target_file, 'w')

    for i, j in pairs:

        lenn = round(math.sqrt(np.sum((locs[i] - locs[j])**2)), round_to)
    
        f.write(str(int(i)) + ',' + str(int(j)) + ',' + str(lenn))
        f.write("\n")
    
    f.close()

    print('Done')


def get_permutations(source, target, thresh, round_to, n_samples):
    """This function gets the permutations given information about edges
    
    source is the csv file with information of edges in the original data set.
    It should be in the format v1, v2, length (see function get_edges_from_locs)
    Note: The lengths are assumed to already have been rounded up
    
 source   target is the location where permutations are made: target+'/Permutations'
    """

    print('Getting permutations')

    perm_dirr = target+'/Permutations'
    if not os.path.isdir(perm_dirr):
        os.mkdir(perm_dirr)
    
    verts = []
    lens = []
    
    test_count = 0
    
    f = open(source, 'r')
    for line in f:
        line = line.strip('\n')
        line = line.split(',')
        this_lenn = (float)(line[2])
        if this_lenn > thresh:
            continue

        verts.append([(int)(float(line[0])), (int)(float(line[1]))])
        lens.append(round(this_lenn, round_to))
        #lens.append(round(math.sqrt((float)(line[2])), 2))
        #unrounded_lens.append(math.sqrt((float)(line[2])))
    
        #if test_count > 20000:
        #    break
        #test_count += 1
    
    f.close()
    
    verts = np.array(verts, dtype=int)
    lens = np.array(lens)
    
    sort_idxs = np.argsort(lens)
    
    verts = verts[sort_idxs]
    lens = lens[sort_idxs]
    
    #plt.hist(np.log10(unrounded_lens), bins=100)
    #plt.yscale('log')
    #plt.show()
    
    u, indices = np.unique(lens, return_index=True)
    
    n_uni = len(u)
    
    for iters in range(n_samples):

        print('making perm', iters, end='\r')

        this_perm_dirr = perm_dirr+'/perm_'+str(iters)
        if not os.path.isdir(this_perm_dirr):
            os.mkdir(this_perm_dirr)

        perm_verts = np.copy(verts)

        f = open(this_perm_dirr+'/perm_'+str(iters)+'.csv', 'w')
        for idx in range(n_uni-1):
            np.random.shuffle(perm_verts[indices[idx]:indices[idx+1]])
    
            for ii in range(indices[idx], indices[idx+1]):
                f.write(str(perm_verts[ii, 0])+','+str(perm_verts[ii, 1])+','+str(lens[ii]))
                f.write('\n')
    
        f.close()
        #input('w')

    print('\n')


#def reduce_with_trivial_triangles(source, target, thresh):
#"""This function reduces a 1-cycle with triangles with diameter <= thresh
#
#This function assumes that edge file has been created (see function get_edges_from_locs)
#"""

# The idea is that the length will shorten only if two consecutive edges are in a trivial triangle


# First get consecutive edges








def get_max_max_all_minimal(n_sig, source, target\
                            , lower_thresh, thresh, epsil\
                            , filetype, threads, analyze_dim):


    # Compute PH till thresh+epsil (do not compute minimal cycles)

    import pydory as dory
    #NOTE: birth_thresh here is thresh
    dory.compute_PH(source, lower_thresh, thresh+epsil, filetype, threads\
                    , target, analyze_dim, 1, 1, thresh, 1)

    del dory
    gc.collect()

    pers = np.loadtxt(target+'H'+str(analyze_dim)+'_pers_data.txt', delimiter=',')

    if not len(pers):
        return 0

    # Check if maximum barcode length/persistence is less than epsil
    # In which case, return 0 to indicate that no significant feature in this
    # perturbation
    if pers.ndim == 1:
        pers = pers.reshape((1, 2))

    pers[pers[:, 1] == -1, 1] = thresh+epsil

    diff = pers[:, 1] - pers[:, 0]

    sigg = np.argwhere(diff >= epsil).flatten()

    # The number of sig features should be the same
    # as in the unperturbed set
    if len(sigg) != n_sig:
        #print('Number of sig features not equal')
        return 0


    # DO NOT NEED TO DO THIS ANYMORE BECAUSE NOW BIRTH_THRESH IS DEFINED AS THRESH (SEE ABOVE
    # compute_PH command)
    ## Compute PH till thresh -- this is because we will compute minimal cycles here and
    ## we do not want to shorten with cycles born at greater than thresh
    #import pydory as dory
    #dory.compute_PH(source, lower_thresh, thresh, filetype, threads\
    #                , target, analyze_dim, 1, 1, birth_thresh, 1)
    #del dory


    locs = np.loadtxt(source, delimiter=',')
    
    # Go over all minimal cycles, and find max-max of pairwise distance
    cyc_file = open(target+'minimal_V_birth_H'+str(analyze_dim)+'.txt'\
                        , 'r')

    max_len = 0

    for cyc in cyc_file:
        cyc = cyc.split(',')
        cyc = cyc[:-1]
        cyc = [int(x) for x in cyc]


        # Get pairwise distances of pts in cyc
        pdist = sklearn.metrics.pairwise_distances(locs[cyc])

        if np.amax(pdist) > max_len:
            max_len = np.amax(pdist)




    return max_len



def get_tetra_graph(edges_file, thresh):


    print('Getting tetra...')

    f = open(edges_file, 'r')

    neighbors_dict = {}
    all_v = frozenset({})

    for line in f:
        line = line.split(',')

        v1 = int(line[0])
        v2 = int(line[1])
        lenn = float(line[2])

        if lenn <= thresh:
            if v1 not in neighbors_dict:
                neighbors_dict[v1] = frozenset({v2})
            else:
                neighbors_dict[v1] = neighbors_dict[v1].union(frozenset({v2}))

            if v2 not in neighbors_dict:
                neighbors_dict[v2] = frozenset({v1})
            else:
                neighbors_dict[v2] = neighbors_dict[v2].union(frozenset({v1}))

            all_v = all_v.union(frozenset({v1, v2}))

    
    all_v = list(all_v)

    tetra_list = []

    for (v1, v2, v3) in itertools.combinations(all_v, 3):
        common = neighbors_dict[v1].intersection(neighbors_dict[v2])
        if not len(common):
            continue

        common = common.intersection(neighbors_dict[v3])
        if not len(common):
            continue

        for vv in common:
            tetra_list.append(frozenset({v1, v2, v3, vv}))


    return tetra_list

    
    
#def reduce_with_trivial_tetra(boundary, edges_file):
#
#    # Step 1. Create neighborhood info up to thresh
#    # Note that we already ahave edge info, so, do not have to
#    # make that again
#
#    f = open(edges_file, 'r')
#
#    for line in f:
#        line = line.split(',')



def smoothen_trivial_tetra(boundary, G):

    # 1. Make a list of all triangles in boundary, format [{v1, v2, v3}]

    triangles = []

    idx = 0
    lenn = len(boundary)

    plot_triangles = []

    #print(len(boundary))

    while idx < lenn:
        v1 = boundary[idx]
        v2 = boundary[idx+1]
        v3 = boundary[idx+2]

        plot_triangles.append([v1, v2, v3])
        
        triangles.append(frozenset({v1, v2, v3}))


        idx += 3

    #mlab.triangular_mesh(locs[:,0],locs[:,1],locs[:,2]\
    #                , plot_triangles, color=(1,223/255,120/255), opacity=0.5)

    #mlab.triangular_mesh(locs[:,0],locs[:,1],locs[:,2]\
    #                , plot_triangles, color=(0,128/255,255/255), opacity=0.5)

    #mlab.show()

    # Make dual graph
    dual_graph = nx.Graph()

    # Add triangles as nodes
    for t in triangles:

        diam = 0
        for v1, v2 in itertools.combinations(list(t), 2):
            lenn = G[v1][v2]['weight']
            #lenn = np.sum((locs[v1] - locs[v2])**2)
            if lenn > diam:
                diam = lenn

        #diam = math.sqrt(diam)
        dual_graph.add_node(t, diam = diam)

    # Add adjacent triangles as edges of the dual graph
    for t1, t2 in itertools.combinations(triangles, 2):
        t3 = t1.intersection(t2)
        assert len(t3) < 3, 'how is it more than or equal to 3?'
        if len(t3) == 2:
            dual_graph.add_edge(t1, t2)


    update = 1
    while (update):

        #####################
        #### FOR TESTING

        ###plot_triangles = []
        ###for node in dual_graph.nodes:
        ###    plot_triangles.append(list(node))

        ###mlab.figure(figure=None, bgcolor=(0,0,0), fgcolor=None\
        ###                        , engine=None, size=(400, 350))
        ###
        ###mlab.points3d(locs[:,0]\
        ###                        , locs[:,1]\
        ###                        , locs[:,2]\
        ###                        , scale_factor = 0.2)

        ###mlab.triangular_mesh(locs[:,0],locs[:,1],locs[:,2]\
        ###                , plot_triangles, color=(1,223/255,120/255), opacity=0.3)

        ###mlab.show()

        #####################


        update = 0

        # Get all tetraherons, a clique of t1, t2, t3, t4

        # Fist, get all maximal cliques
        #print('Finding cliques')
        cliques = nx.find_cliques(dual_graph)

        # Remove the cliques that have less than 3 nodes
        useful_cliques = [x for x in cliques if len(x) > 2]

        smallest_diam = math.inf
        small_tetra = []
        
        ## IGNORE TETRA REMOVAL

        ## Now, go over cliques of len > 3 to get tetrahedrons
        #for cl in useful_cliques:
        #    if len(cl) < 4:
        #        continue

        #    for t1, t2, t3, t4 in itertools.combinations(list(cl), 4):
        #        diam = min(dual_graph.nodes[t1]['diam']\
        #                , dual_graph.nodes[t2]['diam']\
        #                , dual_graph.nodes[t3]['diam']\
        #                , dual_graph.nodes[t4]['diam']\
        #                )
        #        if diam < smallest_diam:
        #            smallest_diam = diam
        #            small_tetra = [t1, t2, t3,t4]

        ## if there was a smallest tetra, remove it from dual graph
        #if smallest_diam != math.inf:
        #    dual_graph.remove_nodes_from(small_tetra)
        #    update = 1
        #    continue

        # No tetra found in boundary, so, look for triangles now
        for cl in useful_cliques:
            for t1, t2, t3 in itertools.combinations(list(cl), 3):

                # This should not be a tetra
                # So, check if new_face is already in the clique or not
                all_verts = t1.union(t2)
                all_verts = all_verts.union(t3)

                new_face = all_verts - t1
                new_face = new_face.union(all_verts - t2)
                new_face = new_face.union(all_verts - t3)

                if new_face in list(cl):
                    # This means tetrahedron [t1, t2, t3, new_face] is in the boundary. Skip it.
                    continue
                
                # new_face is not in the boundary

                # Find diameter of new_face:
                new_face_diam = 0
                for v1, v2 in itertools.combinations(list(new_face), 2):

                    # If this edge is not in G, then it is more than thresh because of how G is
                    # constructed
                    if not G.has_edge(v1, v2):
                        new_face_diam = math.inf
                        break
                    lenn = G[v1][v2]['weight']
                    #lenn = np.sum((locs[v1] - locs[v2])**2)
                    if lenn > new_face_diam:
                        new_face_diam = lenn
                #new_face_diam = math.sqrt(new_face_diam)

                # Diameter of the tetrahedron is the largest of diameter of its faces
                # (Alternately, could have looked at lengths of all edges in the tetrahedron)
                diam = max(dual_graph.nodes[t1]['diam']\
                        , dual_graph.nodes[t2]['diam']\
                        , dual_graph.nodes[t3]['diam']\
                        , new_face_diam\
                        )

                # Consider only if the diameter of this tetrahedron is <= thresh
                if (diam < smallest_diam):
                    smallest_diam = diam
                    small_tetra = [t1, t2, t3]
                    small_new_face = new_face
                    small_new_face_diam = new_face_diam

        # If there was a smallest triangle with diameter, then remove t1,t2,t3 from dual, BUT...
        # ...have to add the fourth face/triangle
        if smallest_diam != math.inf:

            t1, t2, t3 = small_tetra

            new_face = small_new_face

            # Add the new face
            #print('adding', new_face)
            dual_graph.add_node(new_face, diam = small_new_face_diam)

            #mlab.figure(figure=None, bgcolor=(0,0,0), fgcolor=None\
            #                        , engine=None, size=(400, 350))

            #mlab.triangular_mesh(locs[:,0],locs[:,1],locs[:,2]\
            #        ,[list(new_face)], color=(0,0,1), opacity=0.8)

            # Go over neighbors of t1, t2, t3
            neigh = list(dual_graph.neighbors(t1))
            neigh += list(dual_graph.neighbors(t2))
            neigh += list(dual_graph.neighbors(t3))
            neigh = frozenset(neigh)

            #print('t1', t1)
            #print('t2', t2)
            #print('t3', t3)
            #print('new_face', new_face)

            for n in neigh:
                kk = n.intersection(new_face)
                #print(new_face, n, kk)
                assert len(kk) < 3, 'what??????'
                if len(kk) == 2:
                    dual_graph.add_edge(new_face, n)

            # Remove t1, t2, t3
            #print('removing', t1, t2, t3)
            dual_graph.remove_nodes_from([t1, t2, t3])

            
            #######################
            ##### For testing
            ###triangles = []
            ###triangles.append(list(t1))
            ###triangles.append(list(t2))
            ###triangles.append(list(t3))

            ###mlab.triangular_mesh(locs[:,0],locs[:,1],locs[:,2]\
            ###        ,triangles, color=(1,0,0), opacity=0.3)

            ###
            ###mlab.points3d(locs[:,0]\
            ###                        , locs[:,1]\
            ###                        , locs[:,2]\
            ###                        , scale_factor = 0.2)

            ###mlab.triangular_mesh(locs[:,0],locs[:,1],locs[:,2]\
            ###                , plot_triangles, color=(1,223/255,120/255), opacity=0.3)
            ###mlab.show()

            update = 1
            continue


    # Remove 'disconnected' tetrahedrons

    dual_new = []

    for comp in nx.connected_components(dual_graph):
        # If there are <= 4 triangle -> there is no non-trivial boundary
        if len(comp) <= 4:
            #print('skipping these points')
            continue
        dual_new += comp

    dual_graph = dual_graph.subgraph(dual_new)


    new_boundary_list = []

    for node in dual_graph.nodes:
        new_boundary_list += list(node)

    
    #mlab.figure(figure=None, bgcolor=(0,0,0), fgcolor=None\
    #                        , engine=None, size=(400, 350))
    
    #mlab.points3d(locs[:,0]\
    #                        , locs[:,1]\
    #                        , locs[:,2]\
    #                        , scale_factor = 0.2)
    
    #plot_H2_boundary(locs, new_boundary_list, mlab, (238/255, 75/255, 43/255), 0.5)

    #mlab.show()

    #exit()

    return new_boundary_list




def plot_H2_boundary(locs, boundary, mlab, color, opacity):


    triangles = []
    
    i = 0
    while(i < len(boundary)-1):
    
        triangles.append((boundary[i], boundary[i+1], boundary[i+2]))
    
    
        i += 3
    
    
    mlab.triangular_mesh(locs[:,0],locs[:,1],locs[:,2]\
                    ,triangles, color=color, opacity=opacity)
    



@njit
def cart_to_sphere(locs, ptsnew):


    ptsnew[:, 0] = locs[:,0]**2 + locs[:,1]**2
    ptsnew[:, 1] = np.arctan2(np.sqrt(ptsnew[:, 0]), locs[:,2]**2) # for elevation angle defined from Z-axis down

    ptsnew[:,0] = np.sqrt(ptsnew[:, 0] + locs[:, 2]**2)
    ptsnew[:,2] = np.arctan2(locs[:,1], locs[:,0])

    return ptsnew


def get_uniform_measure_phi_theta(pts):

    #plt.cla()
    #plt.clf()
    #plt.scatter(pts[:,0], pts[:,1])
    #plt.show()
    #exit()

    #kernel = stats.gaussian_kde(pts.T)

    #xmin = np.min(pts[:,0])
    #xmax = np.max(pts[:,0])

    #ymin = np.min(pts[:,1])
    #ymax = np.max(pts[:,1])

    #X, Y = np.mgrid[xmin:xmax:100j, ymin:ymax:100j]

    #positions = np.vstack([X.ravel(), Y.ravel()])
    #Z = np.reshape(kernel(positions).T, X.shape)


    #fig, ax = plt.subplots()

    #ax.imshow(np.rot90(Z), cmap=plt.cm.gist_earth_r,
    #           extent=[xmin, xmax, ymin, ymax])
    #ax.plot(pts[:,0], pts[:,1], 'k.', markersize=3, alpha=0.5)
    #ax.set_xlim([xmin, xmax])
    #ax.set_ylim([ymin, ymax])
    #plt.show()
    #exit()

    l_bounds = [0, -math.pi]
    u_bounds = [math.pi, math.pi]

    space_1 = qmc.scale(pts, l_bounds, u_bounds, reverse=True)

    return qmc.discrepancy(space_1, method='WD', workers = -1)


    #delta = [0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 1, 0.1\
    #        , 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5\
    #        , 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.9, 0.95, 1]
    #measure = []

    #for dd in delta:
    #    n_theta = int(2*math.pi/dd) + 1
    #    n_phi = int(math.pi/dd) + 1

    #    grid = np.zeros((n_theta, n_phi))

    #    for pt in pts:

    #        ### NOTE Theta is [1] and phi is [0]
    #        theta = pt[1]
    #        phi = pt[0]

    #        theta_idx = int((theta+math.pi)/dd)
    #        phi_idx = int(phi/dd)

    #        grid[theta_idx, phi_idx] += 1


    #    n_tot = n_theta*n_phi

    #    n_filled = np.count_nonzero(grid.flatten())

    #    measure.append(n_filled/n_tot)

    #plt.cla()
    #plt.clf()
    #plt.scatter(delta, measure)
    #plt.show()
    #exit()



def plot_spherical_coords_projection(pts):

    plt.cla()
    plt.clf()
    fig, axs = plt.subplots(1, 3)

    axs[0].scatter(pts[:, 1], pts[:, 0])
    axs[0].set_xlabel('phi')
    axs[0].set_ylabel('radius')
    axs[0].set_ylim([0, np.amax(pts[:, 0])])

    axs[1].scatter(pts[:, 2], pts[:, 0])
    axs[1].set_xlabel('theta')
    axs[1].set_ylabel('radius')
    axs[1].set_ylim([0, np.amax(pts[:, 0])])

    axs[2].scatter(pts[:, 2], pts[:, 1])
    axs[2].set_xlabel('theta')
    axs[2].set_ylabel('phi')


    #axs[0].scatter(sphere_pert_locs[smoothened, 1]\
    #                , sphere_pert_locs[smoothened, 0], marker='+', color='red')

    #axs[1].scatter(sphere_pert_locs[smoothened, 2]\
    #                , sphere_pert_locs[smoothened, 0], marker='+', color='red')

    #axs[2].scatter(sphere_pert_locs[smoothened, 2]\
    #                , sphere_pert_locs[smoothened, 1], marker='+', color='red')

    plt.show()



def get_n_sig(source, target, lower_thresh, thresh, epsil\
                , filetype, threads\
                , analyze_dim):


    import pydory as dory
    dory.compute_PH(source, lower_thresh, thresh+epsil, filetype, threads, target, analyze_dim\
                        , 0, 0, thresh, 1)
    del dory
    gc.collect()

    pers = np.loadtxt(target+'H'+str(analyze_dim)+'_pers_data.txt', delimiter=',')

    if len(pers) == 0:
        return 0

    if pers.ndim == 1:
        pers = pers.reshape((1,2))


    pers[pers[:, 1] == -1, 1] = thresh+epsil

    diff = pers[:,1] - pers[:,0]

    sigg = np.argwhere(diff >= epsil).flatten()

    sig_diff = diff[diff >= epsil]

    print('sig barcode lengths', sig_diff, '(recall epsil is', epsil,')')


    return len(sigg)



# Return an edge graph
# data can be list of edges or locs
# thresh is the birth threshold
# filetype: 1 -- locs, 2 -- list of edges


def get_edge_graph(data, filetype, thresh):


    if filetype == 2:
        # list of edges of form v1, v2, dist
        # csv file expected

        G = nx.Graph()

        ff = open(data, 'r')

        for line in ff:
            line = line.strip('\n')
            line = line.split(',')
            v1 = int(line[0])
            v2 = int(line[1])
            dist = float(line[2])
            if dist > thresh:
                continue
            G.add_edge(v1, v2, weight=dist)

        return G
    
    # Otherwise filetype is locs. Use KDTree

    # DO IT LATER

def smoothen_trivial_triangles(boundary, G):

    # Boundary is already in order
    # [v1, v2, ...., vn]


    # update boundary to reduce with local trivial triangles
    while (1):

        nn = 0

        ssize = len(boundary)

        candidate = [0, 0, math.inf, 0]

        while nn < ssize:

            v1 = boundary[nn]
            v2 = boundary[(nn+2)%ssize]
            if G.has_edge(v1, v2):
                dist = G[v1][v2]['weight']
                if dist < candidate[2]:
                    candidate = [v1, v2, dist, (nn + 1)%ssize]

            nn += 1

        if candidate[2] != math.inf:
            #print('deleting', new_boundary[candidate[3]])
            del boundary[candidate[3]]
            #print('After update')
            #print(new_boundary)
        else:
            break


    return boundary
            











