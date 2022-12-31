import helper_funcs_v2 as hf2
import numpy as np
import pickle

reso = 5000

targets = [\
         '../fa_dsg_ddel_reso'+str(reso)+'_'\
        ,'../fa_dsg_dpn_reso'+str(reso)+'_'\
        ,'../fa_dsg_mnase_reso'+str(reso)+'_'\
        ,'../fa_dsg_ddel_dpn_reso'+str(reso)+'_'\
        ,'../fa_dpn_reso'+str(reso)+'_'\
        ]

labels = [\
        'FA+DSG, Ddel'\
        ,'FA+DSG, Dpn'\
        ,'FA+DSG, MNase'\
        ,'FA+DSG, Ddel-Dpn'\
        ,'FA, Dpn'\
        ]


b_idx = 1
dim = 1

# Test for Ddel and Dpn

main_target1 = 0
main_target2 = 1

target1 = targets[main_target1] + 'threshidx' + str(b_idx) +'_'
target2 = targets[main_target2] + 'threshidx' + str(b_idx) +'_'

# using smoothed
cyc_file = target1+'smoothen_H'+str(dim)+'.txt'
# Get cycs as list of edges
cycs1 = hf2.get_cycs_as_edge_lists(cyc_file)

# using smoothed
cyc_file = target2+'smoothen_H'+str(dim)+'.txt'
# Get cycs as list of edges
cycs2 = hf2.get_cycs_as_edge_lists(cyc_file)

n_cycs1 = len(cycs1)
n_cycs2 = len(cycs2)

pair_peaks = []

dist_mat = np.ones((n_cycs1, n_cycs2))*(-1)

for c_idx1 in range(n_cycs1):

    cyc1 = np.array(cycs1[c_idx1])
    cyc1_m = np.argmax(np.abs(cyc1[:,0] - cyc1[:,1])).flatten()
    cyc1_maxx = cyc1[cyc1_m]

    for c_idx2 in range(n_cycs2):

        print(c_idx1, c_idx2, end='\r')

        cyc2 = np.array(cycs2[c_idx2])
        cyc2_m = np.argmax(np.abs(cyc2[:,0] - cyc2[:,1])).flatten()
        cyc2_maxx = cyc2[cyc2_m]

        done = 0

        for row1 in cyc1_maxx:

            for row2 in cyc2_maxx:

                vall = max(abs(row1[0] - row2[0]), abs(row1[1] - row2[1]))

                # If at least one of the maximal in each loop anchors are close
                # Then consider these loops as close
                if vall < 10:
                    pair_peaks.append([cyc1, cyc2, vall])
                    dist_mat[c_idx1, c_idx2] = vall
                    done = 1
                    break

            if done:
                break

pickle.dump(pair_peaks, open('ddel_dpn_common_PH_maximal.p', 'wb'))






#for idx, source in enumerate(edge_file_name):
#
#    main_target = targets[idx]
#
#    target = main_target + 'threshidx' + str(b_idx) +'_'
#
#
#    ### using greedy-shortening set
#    #cyc_file = target + 'minimal_V_birth_H'+str(dim)+'.txt'
#    
#    # using smoothed
#    cyc_file = target+'smoothen_H'+str(dim)+'.txt'
#    
#    # Get cycs as list of edges
#    cycs = hf2.get_cycs_as_edge_lists(cyc_file)
#
#
#
#
#
