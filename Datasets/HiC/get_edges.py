import numpy as np
import h5py
from numba import njit, prange
import cooler 

@njit
def search_edges(n_elements, bin1, bin2, pixels, weights, thresh, edges, count):
    
    for i in range(n_elements):
        val = pixels[i]
        bin1_e = bin1[i]
        bin2_e = bin2[i]
        if not (bin1_e < bin2_e) :
            continue

        bal_val = val*weights[bin1_e]*weights[bin2_e]
        bal_val = 1/bal_val

        if (bal_val < thresh):
            edges[count, 0] = bin1_e
            edges[count, 1] = bin2_e
            edges[count, 2] = bal_val
            count += 1

    return edges, count


def save_edges(h5, pwd, exp_name, thresh):

    print("\nFinding edges for", exp_name)
    edge_file = pwd+'/'+exp_name+'_'+str(thresh)+'.csv'
    
    n_elements = h5['pixels']['bin1_id'].shape[0]
    weights = h5['bins']['weight'][:]
    
    chunk_size = 50000
    
    edges = np.ones((500000000, 3))
    
    i_start = 0
    i_end = chunk_size
    
    
    count = 0
    while (i_start < n_elements):
        if (i_end > n_elements):
            i_end = n_elements
    
        print('i_start/n_elements', i_start/n_elements, 'count = ', count, end='\r')
    
        bin1 = h5['pixels']['bin1_id'][i_start:i_end]
        bin2 = h5['pixels']['bin2_id'][i_start:i_end]
        pixels = h5['pixels']['count'][i_start:i_end]
    
        edges, count = search_edges(i_end-i_start, bin1, bin2, pixels, weights, thresh, edges, count)
    
        i_start += chunk_size
        i_end += chunk_size
    
    
    print("\nSaving edges", count)
    edges = edges[:count]
    np.savetxt(edge_file, edges, delimiter=',', fmt='%.2f')


thresh = 400
h5 = h5py.File('4DNFIFLDVASC.mcool', 'r')
h5 = h5['resolutions']['1000']
exp_name = 'control'
pwd = '../HiC'
save_edges(h5, pwd, exp_name, thresh)

thresh = 400
h5 = h5py.File('4DNFILP99QJS.mcool', 'r')
h5 = h5['resolutions']['1000']
exp_name = 'auxin'
pwd = '../HiC'
save_edges(h5, pwd, exp_name, thresh)
