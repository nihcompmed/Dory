import pydory as dory
import check_connected_helper as cch
import helper_functions as hf 
import helper_funcs_v2 as hf2
import networkx as nx
import matplotlib.pyplot as plt
import os
import numpy as np
import itertools
import seaborn as sns
from sklearn.neighbors import KernelDensity
import scipy
import matplotlib as mpl
import cooler 
import pickle
import math
import networkx as nx
import pandas as pd
from upsetplot import plot as upsetplot


reso = 5000


labels = [\
         'fa_dsg_ddel_5kb'\
        ,'fa_dsg_dpn_5kb'\
        ,'fa_dsg_ddel_dpn_5kb'\
        ,'fa_dsg_mnase_5kb'\
        ,'fa_dpn_5kb'\
        ]

plt_labels = [\
         'FA+DSG, DdeI'\
        ,'FA+DSG, DpnII'\
        ,'FA+DSG, DdeI+DpnII'\
        ,'FA+DSG, Mnase'\
        ,'FA, Dpn'\
        ]


# DO FOR B_IDX = 1 and dim = 1
dim = 1

long_range_thresh = int(math.exp(5)-1)

get_atac = 0
get_ccre = 0

ccre_dict = pickle.load(open('../Biomarkers/ccre_reso'+str(reso)+'.p', 'rb'))

ccre_types = [\
            'Low-DNase'\
            ,'DNase-only'\
            ,'DNase-H3K4me3'\
            ,'dELS'\
            ,'pELS'\
            ,'PLS'\
            ]

data_dirr = '../Biomarkers/'
# BED PEAKS
fnames = [\
        data_dirr+'ENCFF294RSZ.bed'\
        ,data_dirr+'GSE165895_HFFc6_Atacseq.mRp.clN_peaks.narrowPeak'\
        ,data_dirr+'GSE165895_HFFc6_CTCF_CT.mRp.clN_peaks.narrowPeak'\
        ]

flabel = [\
        'ChIPseq_CTCF'\
        ,'Atacseq'\
        ,'CTCF'
        ]

hiccup_files = [\
         '../fa_dsg_ddel_hiccups.bedpe'\
        ,'../fa_dsg_dpn_hiccups.bedpe'\
        ,'../fa_dsg_ddel_dpn_hiccups.bedpe'\
        ,'../fa_dsg_mnase_hiccups.bedpe'\
        ,'../fa_dpn_hiccups.bedpe'\
                ]

hiccups_data = []
for ff in hiccup_files:

    # Result is in base pairs
    # dict[chrom] = [cent1, cent2\
    # , start1, end1\
    # , start2, end2\
    # , c_size, c_count\
    # , c_donut, c_vertical\
    # , c_horizontal, c_lowleft]

    hiccups_data.append(hf2.get_hiccup_peaks(ff))



atac_narrow_file = fnames[1]

# BIGWIG FILES
bigwig_fnames = [\
         data_dirr+'GSE165895_HFFc6_H3K27ac.bigWig'\
        ,data_dirr+'GSE165895_HFFc6_H3K4me3.bigWig'\
        ,data_dirr+'GSE165895_HFFc6_SMC1_CT.mRp.clN.bigWig'\
        ]

bigwig_flabel = [\
        'H3K27ac'\
        ,'H3K4me3'\
        ,'SMC1'\
        ]


# Initialize Pandas multi index
tuples = [\
          (0,0,0,0,0)\
          ,(0,0,0,0,1)\
          ,(0,0,0,1,0)\
          ,(0,0,0,1,1)\
          ,(0,0,1,0,0)\
          ,(0,0,1,0,1)\
          ,(0,0,1,1,0)\
          ,(0,0,1,1,1)\
          ,(0,1,0,0,0)\
          ,(0,1,0,0,1)\
          ,(0,1,0,1,0)\
          ,(0,1,0,1,1)\
          ,(0,1,1,0,0)\
          ,(0,1,1,0,1)\
          ,(0,1,1,1,0)\
          ,(0,1,1,1,1)\
          ,(1,0,0,0,0)\
          ,(1,0,0,0,1)\
          ,(1,0,0,1,0)\
          ,(1,0,0,1,1)\
          ,(1,0,1,0,0)\
          ,(1,0,1,0,1)\
          ,(1,0,1,1,0)\
          ,(1,0,1,1,1)\
          ,(1,1,0,0,0)\
          ,(1,1,0,0,1)\
          ,(1,1,0,1,0)\
          ,(1,1,0,1,1)\
          ,(1,1,1,0,0)\
          ,(1,1,1,0,1)\
          ,(1,1,1,1,0)\
          ,(1,1,1,1,1)\
        ]

index = pd.MultiIndex.from_tuples(tuples, names=plt_labels)


# Initialize set intersection counts to 0
s = pd.Series(np.zeros(len(tuples), dtype=int)\
    , index=index)

graph_dict = dict()

loop_info = dict()
loop_id = 0

full_loop = dict()

for e_idx, exp in enumerate(labels):

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

    # b_idx = 1 here means tau_u^2
    b_idx = 1

    print(exp, b_idx)
    
    # ANCHOR THRESH IS 1 if birth idx is for gene distance 1
    # ANCHOR THRESH IS 2 if birth idx is for gene distance 2
    anchor_thresh = b_idx + 1
    
    
    t1 = '../cooler_results/PH_results/'+exp+'_gdist'+str(b_idx+1)+'_ALL_scaled_'
    
    # Assuming they are checked for connectedness, hence
    # Cycles are stored as v0, v1, v2... where v_i, v_{i+1} is an edge
    cyc_file = t1+'minimal_V_birth_H'+str(dim)+'.txt'
    
    
    cycs = open(cyc_file, 'r')
    
    
    #progress = 0

    
    for line in cycs:

        long_range_anchors = []
    
        #print(progress, end='\r')
        #progress += 1
    
        line = line.split(',')
        line = line[:-1]
        line = [int(x) for x in line]
        line.append(line[0])
    
        line = np.array(line)
    
        cyc_chroms = []

        for binn in line:
            info = hf2.get_chrom(binn, chrom_start_bins_list, chrom_names)
            cyc_chroms.append(info)

        #long_range_anchor_list = []

        nn = 0
        while nn < len(line)-1:
            c1 = cyc_chroms[nn][0]
            bin1 = cyc_chroms[nn][1]
    
            c2 = cyc_chroms[nn+1][0]
            bin2 = cyc_chroms[nn+1][1]

            if c1==c2:
                this_dist = abs(bin2 - bin1)
                if this_dist > long_range_thresh:
                    long_range_anchors.append([bin1, bin2])
                    #if c1 not in graph_dict:
                    #    graph_dict[c1] = nx.Graph()
                    #long_range_anchor_list.append([c1, bin1, bin2])
                    #graph_dict[c1].add_node((e_idx, bin1, bin2))
    
            nn += 1
        
        ## long-range interaction
        #if len(long_range_anchor_list) > 0:
        #    print(long_range_anchor_list)
        
        if len(long_range_anchors) > 1:
            if c1 not in graph_dict:
                graph_dict[c1] = nx.Graph()

            loop_info[loop_id] = long_range_anchors

            graph_dict[c1].add_node((e_idx, loop_id))

            full_loop[loop_id] = line

            loop_id += 1

common_loop_multiple_long_range = []

for key in graph_dict:
    print(key)

    for main_a1, main_a2 in itertools.combinations(list(graph_dict[key].nodes), 2):

        loop1_id = main_a1[1]
        loop2_id = main_a2[1]

        loop1 = loop_info[loop1_id]
        loop2 = loop_info[loop2_id]

        #print(loop1, loop2)

        n_l1 = len(loop1)
        n_l2 = len(loop2)

        if n_l1 !=  n_l2:
            continue

        distmat = np.zeros((n_l1, n_l2), dtype=int)

        for a1_idx, a1 in enumerate(loop1):

            a11 = min(a1[0], a1[1])
            a12 = max(a1[0], a1[1])

            for a2_idx, a2 in enumerate(loop2):

                a21 = min(a2[0], a2[1])
                a22 = max(a2[0], a2[1])

                ddist = max(abs(a11-a21), abs(a12-a22))
                distmat[a1_idx, a2_idx] = ddist

        val1 = np.max(np.min(distmat, axis=0))
        val2 = np.max(np.min(distmat, axis=1))

        vall = max(val1, val2)

        # threshold b_idx+ 1 = 2
        if vall < 3:
            #print('adding edge', main_a1, main_a2)
            graph_dict[key].add_edge(main_a1, main_a2)


# Initialize set intersection counts to 0
s = pd.Series(np.zeros(len(tuples), dtype=int)\
    , index=index)

for chrom in graph_dict:

    cliques = nx.find_cliques(graph_dict[chrom])

    for cl in cliques:

        this_tuple = np.zeros((5,), dtype=int)

        for nn in cl:
            this_tuple[nn[0]] = 1

        this_tuple = tuple(this_tuple)

        s[this_tuple] += 1

for tt in tuples:
    s[tt] = math.log1p(1+s[tt])

upsetplot(s)
plt.ylabel('log(1+Intersection size)')
#plt.show()
plt.savefig('long_range_multiple_loop_valency_bidx'+str(b_idx+1)+'.pdf', dpi=600)

plt.cla()
plt.clf()
plt.close()

loops_to_check = []

idd = 0

for chrom in graph_dict:

    cliques = nx.find_cliques(graph_dict[chrom])

    for cl in cliques:

        this_tuple = np.zeros((5,), dtype=int)

        for nn in cl:
            this_tuple[nn[0]] = 1

        this_tuple = tuple(this_tuple)

        if this_tuple == (1,1,1,1,1):

            idd += 1

            for nn in cl:
                e_idx = nn[0]
                print(e_idx)

                exp = labels[e_idx]

                cfile = '../cooler_files/'+exp + '.cool'

                #graph_file = labels[e_idx]+'_edge_graph_5kb.gml'
                #print('reading graph...')
                #G = nx.read_gml(graph_file)
                #print('read graph')

                loop_id = nn[1]
                this_loop = full_loop[loop_id]
                loops_to_check.append(this_loop)

                #ax = plt.gca()
                #
                #ax = hf2.plot_loop_on_matrix(ax, this_loop, cfile)

                #ax.axis('off')

                ##plt.savefig('multiple_long_range_id'+str(idd)+'_exp'+exp+'.pdf', dpi=600)
                #plt.savefig('multiple_long_range_id'+str(idd)+'_exp'+exp+'.png', dpi=300)
                #plt.cla()
                #plt.clf()
                #plt.close()



# atac-seq enrichment
data_dirr = '../Biomarkers/'


# BED PEAKS
fnames = [\
        data_dirr+'ENCFF294RSZ.bed'\
        ,data_dirr+'GSE165895_HFFc6_Atacseq.mRp.clN_peaks.narrowPeak'\
        ,data_dirr+'GSE165895_HFFc6_CTCF_CT.mRp.clN_peaks.narrowPeak'\
        ]

flabel = [\
        'ChIPseq_CTCF'\
        ,'Atacseq'\
        ,'CTCF'
        ]

atac_narrow_file = fnames[1]

# BIGWIG FILES
bigwig_fnames = [\
         data_dirr+'GSE165895_HFFc6_H3K27ac.bigWig'\
        ,data_dirr+'GSE165895_HFFc6_H3K4me3.bigWig'\
        ,data_dirr+'GSE165895_HFFc6_SMC1_CT.mRp.clN.bigWig'\
        ]

bigwig_flabel = [\
        'H3K27ac'\
        ,'H3K4me3'\
        ,'SMC1'\
        ]

print('finding atac peaks...')
unique_bed_lines = []

for loop in loops_to_check:

    nn = 0
    while nn < len(loop)-1:
        b1 = loop[nn]
        b2 = loop[nn+1]

        if abs(b2 - b1) < 3:
            nn += 1
            continue

        chrr_1, binn_1 = hf2.get_chrom(b1, chrom_start_bins_list, chrom_names)

        start = binn_1*reso
        end = (binn_1+1)*reso
        bed_lines = hf2.get_narrow_range(atac_narrow_file, 'chr'+chrr_1, start, end)

        for ll in bed_lines:
            if ll not in unique_bed_lines:
                unique_bed_lines.append(ll)

        chrr_1, binn_1 = hf2.get_chrom(b2, chrom_start_bins_list, chrom_names)

        start = binn_1*reso
        end = (binn_1+1)*reso
        bed_lines = hf2.get_narrow_range(atac_narrow_file, 'chr'+chrr_1, start, end)

        for ll in bed_lines:
            if ll not in unique_bed_lines:
                unique_bed_lines.append(ll)

        nn += 1


bed_file = 'multiple_long_range_common_to_all_atac.bed'
bed = open(bed_file, 'w')

for ll in unique_bed_lines:
    for entry in ll[:-1]:
        bed.write(entry+'\t')
    bed.write(ll[-1])

bed.close()

exit()

#    print(key, list(graph_dict[key].nodes))

#anchor_dict = dict()
#
#
#for chrom in graph_dict:
#
#    for n1, n2 in itertools.combinations(list(graph_dict[chrom].nodes), 2):
#    
#        b11 = n1[1]
#        b12 = n1[2]
#    
#        b21 = n2[1]
#        b22 = n2[2]
#    
#        ddist = max(abs(b11-b21), abs(b12-b22))
#    
#        if ddist < anchor_thresh + 1:
#            graph_dict[chrom].add_edge(n1, n2)
#    
#
#    cliques = nx.find_cliques(graph_dict[chrom])
#    
#    for cl in cliques:
#    
#        this_tuple = np.zeros((5,), dtype=int)
#        for nn in cl:
#            #exps.append(nn[0])
#            this_tuple[nn[0]] = 1
#        
#        this_tuple = tuple(this_tuple)
#        s[this_tuple] += 1
#
#
#        if this_tuple not in anchor_dict:
#            anchor_dict[this_tuple] = dict()
#
#
#        if chrom not in anchor_dict[this_tuple]:
#            anchor_dict[this_tuple][chrom] = []
#
#        anchor_dict[this_tuple][chrom].append(cl)
#
#
#
#for tt in tuples:
#    s[tt] = math.log1p(s[tt])
#
#
##upsetplot(s)
##plt.ylabel('log(1 + intersection size)')
##plt.show()
#
#
## Get information for ATAC
#if get_atac:
#
#    for this_tuple in anchor_dict:
#    
#        strtt = ''
#        for tt in this_tuple:
#            strtt = strtt+str(tt)
#    
#        unique_bed_lines = []
#    
#        print('Doing tuple', this_tuple)
#    
#        bed_file = 'long_range_combination_'+strtt+'_atac_dim'+str(dim)+'_bidx'+str(b_idx+1)+'.bed'
#    
#        bed = open(bed_file, 'w')
#
#        for chrom in anchor_dict[this_tuple]:
#    
#            for cl in anchor_dict[this_tuple][chrom]:
#    
#                print(cl)
#    
#                for anchor in cl:
#    
#                    bin1 = anchor[1]
#                    bin2 = anchor[2]
#    
#                    start = bin1*reso
#                    end = (bin1+1)*reso
#    
#                    bed_lines = hf2.get_narrow_range(atac_narrow_file, 'chr'+str(chrom), start, end)
#                    for ll in bed_lines:
#                        if ll not in unique_bed_lines:
#                            unique_bed_lines.append(ll)
#    
#                    start = bin2*reso
#                    end = (bin2+1)*reso
#    
#                    bed_lines = hf2.get_narrow_range(atac_narrow_file, 'chr'+str(chrom), start, end)
#                    for ll in bed_lines:
#                        if ll not in unique_bed_lines:
#                            unique_bed_lines.append(ll)
#    
#            for ll in unique_bed_lines:
#                for entry in ll[:-1]:
#                    bed.write(entry+'\t')
#                bed.write(ll[-1])
#    
#        bed.close()
#    
#    
#    
## GET cCRE information
#
#
#if get_ccre:
#
#    for this_tuple in anchor_dict:
#    
#        strtt = ''
#        for tt in this_tuple:
#            strtt = strtt+str(tt)
#    
#        pie_chart = np.array([0,0])
#    
#        for chrom in anchor_dict[this_tuple]:
#    
#            for cl in anchor_dict[this_tuple][chrom]:
#    
#                flag = 0
#    
#                for anchor in cl:
#    
#                    this_ccre_dict = ccre_dict['chr'+chrom]
#    
#                    binn1 = anchor[1]
#                    start1 = binn1
#                    end1 = binn1 + 1
#    
#                    binn2 = anchor[2]
#                    start2 = binn2
#                    end2 = binn2 + 1
#    
#    
#                    #cre_set = []
#                    for key in this_ccre_dict:
#    
#                        ccre_start = key[0]
#                        ccre_end = key[1]
#    
#                        # If either ccre start or end is in the bin of the anchor
#                        if\
#                            (ccre_start >= start1 and ccre_start <= end1)\
#                            or (ccre_end >= start1 and ccre_end <= end1)\
#                            or (ccre_start >= start2 and ccre_start <= end2)\
#                            or (ccre_end >= start2 and ccre_end <= end2)\
#                            :
#    
#                            for cc in this_ccre_dict[key]:
#                                c_idx = ccre_types.index(cc)
#                                #print(c_idx)
#                                #this_tuple[c_idx] = 1
#                                #cre_set.append(c_idx)
#                                if c_idx != 0:
#                                    flag = 1
#                                    pie_chart[1] += 1
#                                    break
#    
#                if not flag:
#                    pie_chart[0] += 1
#    
#        #print(this_tuple, pie_chart) 
#        
#        plt.pie(pie_chart, labels=['no cCRE', 'has cCRE'], colors=['#D41159', '#1A85FF'])
#        plt.title('Number of cliques '+str(np.sum(pie_chart)))
#    
#        out = 'ccre_cis_long_range_pie_chart_combination_'+strtt+'_dim'+str(dim)+'_bidx'+str(b_idx)
#        
#        plt.savefig(out+'.pdf', dpi=600)
#        plt.cla()
#        plt.clf()
#        plt.close()
#
#
## check genes
#
#dict_GO = pickle.load(open('dict_GO.data', 'rb'))
#
##print(dict_GO)
#
#genes_to_bins1kb = pickle.load(open('genes_to_bins1kb.data', 'rb'))
#
#
#gene_to_GO = pickle.load(open('dict_gene_to_GO.data', 'rb'))
#
#
#bins_to_genes = pickle.load(open('all_bins_to_gene_info_hg38_reso'+str(reso)+'.p', 'rb'))
#
#
##for this_tuple in anchor_dict:
#
#this_tuple = (1,1,1,1,1)
#strtt = ''
#for tt in this_tuple:
#    strtt = strtt+str(tt)
#
#print('Doing tuple', this_tuple)
#
#
#for chrom in anchor_dict[this_tuple]:
#
#    this_chrom_bin_gene_dict = bins_to_genes['chr'+str(chrom)]
#
#
#    for cl in anchor_dict[this_tuple][chrom]:
#
#
#        has_GO1 = 0
#        has_GO2 = 0
#
#        found_genes1 = []
#        found_genes2 = []
#
#        for anchor in cl:
#        
#            bin1 = anchor[1]
#            bin2 = anchor[2]
#
#            if bin1 in this_chrom_bin_gene_dict:
#                genes1 = list(frozenset(this_chrom_bin_gene_dict[bin1]))
#                genes1 = [x.strip('\"') for x in genes1]
#                for gg in genes1:
#                    if gg in gene_to_GO:
#                        found_genes1.append(gg)
#                        has_GO1 = 1
#                
#            if not has_GO1:
#                continue
#
#            if bin2 in this_chrom_bin_gene_dict:
#                genes2 = list(frozenset(this_chrom_bin_gene_dict[bin2]))
#                genes2 = [x.strip('\"') for x in genes2]
#                for gg in genes2:
#                    if gg in gene_to_GO:
#                        found_genes2.append(gg)
#                        has_GO2 = 1
#
#            if not has_GO2:
#                continue
#
#        if not has_GO1 or not has_GO2:
#            continue
#
#        print('FOUND A PAIR', genes1, genes2)
#
#        found_genes1 = list(frozenset(found_genes1))
#        found_genes2 = list(frozenset(found_genes2))
#        
#        genes = found_genes1 + found_genes2
#
#        b1 = bin1*reso
#        b2 = bin2*reso
#
#        ddist = dict()
#        for e_idx, flag in enumerate(this_tuple):
#            if flag:
#                ddist[e_idx] = math.inf
#        
#        for e_idx, flag in enumerate(this_tuple):
#            print(e_idx, flag)
#            if flag:
#                for hiccup in hiccups_data[e_idx]['chr'+chrom]:
#
#                    c1 = min(hiccup[0], hiccup[1])
#                    c2 = max(hiccup[0], hiccup[1])
#
#                    val = max(abs(b1-c1), abs(b2-c2))
#
#                    if val < ddist[e_idx]:
#                        ddist[e_idx] = val/reso
#                        this_g_dist = abs(b2-b1)/reso
#
#        print(ddist)
#
#        hf2.plot_gene_GO(genes, dict_GO, gene_to_GO, genes_to_bins1kb)
#
#        exit()
#
#
#
#
#            
#
#
#
#
#
