import matplotlib.pyplot as plt
import numpy as np
import pandas
import math
import h5py
from numba import njit, prange
from numba import njit
import pickle
import cooler 
import seaborn as sns
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
import helper_funcs_v2 as hf2
import matplotlib.cm as cm


def get_chrom(binn, chrom_start_bins_list, chrom_names):

    for idx, start in enumerate(chrom_start_bins_list):

        if binn < start:
            return chrom_names[idx-1], binn-chrom_start_bins_list[idx-1]

PH_pars = pickle.load(open('../cooler_results/edge_files/PH_pars_info_5kb.p', 'rb'))

hiccup_files = [\
         '../fa_dsg_ddel_hiccups.bedpe'\
        ,'../fa_dsg_dpn_hiccups.bedpe'\
        ,'../fa_dsg_ddel_dpn_hiccups.bedpe'\
        ,'../fa_dsg_mnase_hiccups.bedpe'\
        ,'../fa_dpn_hiccups.bedpe'\
                ]

labels = [\
         'fa_dsg_ddel_5kb'\
        ,'fa_dsg_dpn_5kb'\
        ,'fa_dsg_ddel_dpn_5kb'\
        ,'fa_dsg_mnase_5kb'\
        ,'fa_dpn_5kb'\
        ]

chroms = [
         '1'\
        ,'2'\
        ,'3'\
        ,'4'\
        ,'5'\
        ,'6'\
        ,'7'\
        ,'8'\
        ,'9'\
        ,'10'\
        ,'11'\
        ,'12'\
        ,'13'\
        ,'14'\
        ,'15'\
        ,'16'\
        ,'17'\
        ,'18'\
        ,'19'\
        ,'20'\
        ,'21'\
        ,'22'\
        ,'X'\
        ,'Y'\
        ]

hiccups_dict = dict()

for idx, ff in enumerate(hiccup_files):
    hiccups_dict[labels[idx]] = hf2.get_hiccup_peaks(ff)

reso = 5000

n_chroms = 24

for exp in PH_pars:

    print(exp)
    #exit()

    # cfile
    cfile = '../cooler_files/'+exp + '.cool'

    cf = cooler.Cooler(cfile)

    chrom_start_bins_list = []
    chrom_names = []

    n_chrom = 0
    
    for ch in cf.chromnames:
        
        chrom_names.append(ch)
        beg, end = cf.extent(ch)
        chrom_start_bins_list.append(beg)
        n_chrom += 1

    chrom_start_bins_list.append(cf.extent(cf.chromnames[-1])[1])
    chrom_start_bins = np.array(chrom_start_bins_list)
    


    for idx in range(3):

        print(exp, idx)

        cycle_info = []

        cis_cycles = []

        target = '../cooler_results/PH_results/'+exp+'_gdist'+str(idx+1)+'_ALL_'

        # H1 shortened cycles (and have been checked for connecedness)
        cyc_file = target + 'minimal_V_birth_H1.txt'

        cycs = open(cyc_file, 'r')

        for line in cycs:
            line = line.split(',')
            line = line[:-1]
            line = [int(x) for x in line]
            line = np.array(line)

            # compute max genomic distance along boundary
            max_gdist = max(np.amax(np.abs(line[1:] - line[:-1]))\
                                , abs(line[-1] - line[0]))

            chroms_in_cycle = []
            bin_in_chrom = []
            for binn in line:

                chrr, binn = get_chrom(binn, chrom_start_bins_list, chrom_names)
                chroms_in_cycle.append(chrr)
                bin_in_chrom.append(binn)

            cycle_info.append([line, chroms_in_cycle, bin_in_chrom, max_gdist])

            unique_chroms = list(frozenset(chroms_in_cycle))

            if len(unique_chroms):
                cis_cycles.append([unique_chroms[0], bin_in_chrom, max_gdist])


        save_file =\
                '../cooler_results/PH_results/'+exp+'_gdist'+str(idx+1)+'_ALL_shortened_summary.p'

        pickle.dump(cycle_info, open(save_file, 'wb'))

        # Get distance of nearest hiccup peak

        for cis_cyc in cis_cycles:

            cyc_chrom = cis_cyc[0]
            cyc_bins = cis_cyc[1]
            cyc_gdist = cis_cyc[2]

            if 'chr'+cyc_chrom not in hiccups_dict[exp]:
                continue

            this_hicc_peaks = np.array(hiccups_dict[exp]['chr'+cyc_chrom])

            #print(this_hicc_peaks)

            this_hicc_cents = this_hicc_peaks[:, :2]


            global_minn = math.inf

            # iterate over every edge
            nn = 0
            lenn= len(cyc_bins)
            while (nn < lenn):
                b1 = cyc_bins[nn]
                b2 = cyc_bins[(nn+1)%lenn]

                bb1 = min(b1, b2)
                bb2 = max(b1, b2)

                this_edge = np.array([bb1, bb2])*reso

                ddist = np.abs(this_hicc_cents - this_edge)

                maxx = np.max(ddist, axis=1)

                minn = np.min(maxx)

                global_minn = min(minn, global_minn)

                #print('closest', minn)

                nn += 1

            plt.scatter(cyc_gdist\
                        , global_minn/reso\
                        , color=cm.gnuplot(chroms.index(cyc_chrom)/24)\
                        , alpha=0.8\
                        )
            #print(global_minn/reso, cyc_gdist)
            #input('w')

        plt.xscale('log', base=2)
        plt.yscale('log', base=2)

        plt.xlabel('max genomic dist in boundary')
        plt.ylabel('nearest hiccup peak in 5kb resolution')
        plt.show()
        plt.cla()
        plt.clf()





#        exit()
#
#
#
#            
#
#
#
#        ## Get cycs as list of edges
#        #cycs = hf2.get_cycs_as_edge_lists(cyc_file)
#
#
#        ## 
#
#        #hf2.get_cycs_hiccups_dist(cycs, hiccups_dict[exp], reso)
#
#
#
#
#
#
#
