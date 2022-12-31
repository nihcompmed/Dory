import pickle
import numpy as np
import helper_funcs_v2 as hf2
import math
import seaborn as sns
import matplotlib.pyplot as plt

reso = 5000



edge_file_name = [\
         '../fa_dsg_ddel_chr1_KR_reso'+str(reso)+'.csv'\
        ,'../fa_dsg_dpn_chr1_KR_reso'+str(reso)+'.csv'\
#        ,'../fa_dsn_mnase_chr1_KR_reso'+str(reso)+'.csv'\
        ,'../fa_dsg_ddel_dpn_chr1_KR_reso'+str(reso)+'.csv'\
#        ,'../fa_dpn_chr1_KR_reso'+str(reso)+'.csv'\
        ]


targets = [\
         '../fa_dsg_ddel_reso'+str(reso)+'_'\
        ,'../fa_dsg_dpn_reso'+str(reso)+'_'\
#        ,'../fa_dsg_mnase_reso'+str(reso)+'_'\
        ,'../fa_dsg_ddel_dpn_reso'+str(reso)+'_'\
#        ,'../fa_dpn_reso'+str(reso)+'_'\
        ]

labels = [\
        'FA+DSG, Ddel'\
        ,'FA+DSG, Dpn'\
#        ,'FA+DSG, MNase'\
        ,'FA+DSG, Ddel-Dpn'\
#        ,'FA, Dpn'\
        ]

dim = 1

PH_reso = reso


colors = [\
        '#E69F00'\
        ,'#56B4E9'\
        ,'#009E73'\
        ,'#D55E00'\
        ,'#CC79A7'\
        ]


# Biomarkers
bio_data_dirr = '../Biomarkers/'


bio_flabel = [\
        'ChIPseq_CTCF'\
        ,'Atacseq'\
        ,'CTCF'
        ]




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


for this_label in bio_flabel:

    load_name = bio_data_dirr+this_label+'_reso'+str(reso)+'_parsed.p'
    
    biomarker_dict = pickle.load(open(load_name, 'rb'))
    
    for b_idx in range(3):
    
        fig, axs = plt.subplots(1, 2)
        plt_fname = 'narrowpeaks_PHanchors_'+this_label+'_birth'+str(b_idx)+'.pdf'
    
        for idx, source in enumerate(edge_file_name):
    
            main_target = targets[idx]
    
            target = main_target + 'threshidx' + str(b_idx) +'_'
    
    
            ### using greedy-shortening set
            #cyc_file = target + 'minimal_V_birth_H'+str(dim)+'.txt'
            
            # using smoothed
            cyc_file = target+'smoothen_H'+str(dim)+'.txt'
    
            # Get cycs as list of edges
            # These are on the scale of original resolution
            cycs = hf2.get_cycs_as_edge_lists(cyc_file)
    
            max_anchors = []
    
            # Get max genomic dists along boundary
            cycs_genomic_dists = []
            for cyc in cycs:
                #print(cyc)
                cyc = np.array(cyc)
                ddist = np.abs(cyc[:,0] - cyc[:,1])
    
                maxx = np.argmax(ddist).flatten()
    
                this_max_anchors = cyc[maxx]
    
                max_anchors += list(this_max_anchors)
    
                vval = np.amax(ddist)
                cycs_genomic_dists.append(vval)
    
    
            cycs_genomic_dists = np.array(cycs_genomic_dists)
            #print(len(max_anchors))
            #print(np.max(cycs_genomic_dists))
    
            #exit()
    
            all_near = []
    
            count = 0
    
            for b1, b2 in max_anchors:
    
                print(count, end='\r')
                count += 1
    
                near_bio1 = find_nearest_biomarker(b1, biomarker_dict)
                near_bio2 = find_nearest_biomarker(b2, biomarker_dict)
    
                #all_near += [near_bio1[0], near_bio2[0]]
                all_near.append(min(near_bio1[0], near_bio2[0]))
    
    
            #print(all_near)
    
            all_near = np.array(all_near)
            all_near = np.log2(all_near+1)
    
            sns.distplot(all_near, label = labels[idx], ax=axs[0])
            axs[1].scatter(cycs_genomic_dists, all_near, alpha=0.5)
    
    
        axs[1].set_xscale('log')
        axs[1].set_xlabel('genomic dist')
        axs[1].set_ylabel('nearest ctcf dist (log2)')

        axs[0].set_xlabel('log2(nearest in 5kb resolution + 1)')
        axs[0].set_ylabel('counts')
    
        axs[0].legend()

    
        plt.legend()

        #plt.show()
        plt.savefig(plt_fname, dpi=600)

        plt.cla()
        plt.clf()
        #plt.close()
    
    
    
    
