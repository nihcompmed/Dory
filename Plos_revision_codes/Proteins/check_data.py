import pickle
import persim
import itertools
import networkx as nx
import numpy as np
import matplotlib.pyplot as plt
import time

def get_L0_norm(pdb1, pdb2, pers_dict, thresh, epsil):

    pers1 = pers_dict[pdb1]
    pers2 = pers_dict[pdb2]

    if pers1.ndim == 1:
        pers1 = pers1.reshape((1,2))
    if pers2.ndim == 1:
        pers2 = pers2.reshape((1,2))

    lenn = max(pers1.shape[0], pers2.shape[0])

    pers1[pers1[:,1]==-1, 1] = thresh+epsil
    pers2[pers2[:,1]==-1, 1] = thresh+epsil

    diff1 = pers1[:,1] - pers1[:,0]
    diff2 = pers2[:,1] - pers2[:,0]

    #print(pers1, pers2)
    diff1 = -np.sort(-diff1)
    diff2 = -np.sort(-diff2)

    n_sig1 = len(np.argwhere(diff1 >= epsil).flatten())
    n_sig2 = len(np.argwhere(diff2 >= epsil).flatten())


    diff1 = np.pad(diff1, (0, lenn-len(diff1)), 'constant', constant_values = 0)
    diff2 = np.pad(diff2, (0, lenn-len(diff2)), 'constant', constant_values = 0)


    diff = np.abs(diff1 - diff2)

    return np.amax(diff), n_sig1, n_sig2



#ff = 'all_final_results_v2.p'
#ff = 'all_final_results.p'
ff = 'pers_dict_homologs.p'
#ff = 'shortlisted_sets.p'
#ff = 'sig_pdb.p'

thresh = 10
epsil = 3.5

count = 0

pers_data = pickle.load(open(ff, 'rb'))

homolog_file = 'similar_pdb.csv'

homolog_data = open(homolog_file, 'r')

all_bn = []
all_ws = []
all_sw = []
all_L0 = []

all_times_bn = []
all_times_ws = []
all_times_sw = []
all_times_L0 = []

count = 0

skip_set = 0

for line in homolog_data:
    line = line.split(',')
    line = line[:-1]

    #G = nx.Graph()
    print(count, end='\r')
    count += 1
    #if count == 10:
    #    break

    flag = 1
    for pdb in line:
        if pdb.lower() not in pers_data:
            flag = 0
            break

    if not flag:
        skip_set += 1
        continue


    for p1, p2 in itertools.combinations(line, 2):

        pers1 = pers_data[p1.lower()]
        pers2 = pers_data[p2.lower()]


        pers1[pers1[:,1]==-1, 1] = thresh+epsil
        pers2[pers2[:,1]==-1, 1] = thresh+epsil

        # Bottleneck
        start = time.time()
        bn = persim.bottleneck(pers1, pers2)
        all_times_bn.append(time.time()-start)
        #print(bn)

        # Wasserstein
        start = time.time()
        ws = persim.wasserstein(pers1, pers2)
        all_times_ws.append(time.time()-start)
        #print(ws)

        # Sliced Wasserstein
        start = time.time()
        sw = persim.sliced_wasserstein(pers1, pers2)
        all_times_sw.append(time.time()-start)
        #print(ws)

        start = time.time()
        L0, n_sig1, n_sig2 = get_L0_norm(p1.lower(), p2.lower(), pers_data, thresh, epsil)
        all_times_L0.append(time.time()-start)
        #print(L0_norm)

        all_bn.append(bn)
        all_ws.append(ws)
        all_sw.append(sw)
        all_L0.append(L0)

    #exit()
    #for pdb in line:
    #    print(pdb, pers_data[pdb.lower()])

print('skipped sets', skip_set)

info_dict = dict()
info_dict['vals'] = dict()
info_dict['times'] = dict()

info_dict['vals']['L0'] = all_L0
info_dict['vals']['bn'] = all_bn
info_dict['vals']['ws'] = all_ws
info_dict['vals']['sw'] = all_sw

info_dict['times']['L0'] = all_times_L0
info_dict['times']['bn'] = all_times_bn
info_dict['times']['ws'] = all_times_ws
info_dict['times']['sw'] = all_times_sw

#print(info_dict)


pickle.dump(info_dict, open('PD_distances_pairwise.p', 'wb'))

fig, axs = plt.subplots(1, 3)

axs[0].scatter(all_L0, all_bn, color='blue')
axs[0].set_xlabel('L0 norm')
axs[0].set_ylabel('bottleneck')


axs[1].scatter(all_L0, all_ws, color='blue')
axs[1].set_xlabel('L0 norm')
axs[1].set_ylabel('wasserstein')

axs[2].scatter(all_L0, all_sw, color='blue')
axs[2].set_xlabel('L0 norm')
axs[2].set_ylabel('sliced wasserstein')

plt.tight_layout()

plt.show()
plt.cla()
plt.clf()
plt.close()



#for key in data:
#    count += 1
#    #this = data[key]
#    print(key)
#    #for key2 in this:
#    #    print(key2)
##
##    exit()
#
#print(count)
