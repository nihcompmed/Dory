import helper_funcs_v2 as hf2
import pickle
import numpy as np
import math
import matplotlib.pyplot as plt

reso = 5000

hiccup_files = [\
         '../fa_dsg_ddel_hiccups'\
        ,'../fa_dsg_dpn_hiccups'\
        ,'../fa_dsg_mnase_hiccups'\
        ,'../fa_dsg_ddel_dpn_hiccups'\
        ,'../fa_dpn_hiccups'\
                ]

labels = [\
        'FA+DSG, Ddel'\
        ,'FA+DSG, Dpn'\
        ,'FA+DSG, MNase'\
        ,'FA+DSG, Ddel-Dpn'\
        ,'FA, Dpn'\
        ]

# Initialize dictionary and save it
#chrom = 'chr1'
#
#all_peaks = dict()
#
#for idx, ff in enumerate(hiccup_files):
#
#    all_peaks[labels[idx]] = hf2.get_hiccup_peaks_dict(ff+'.bedpe', chrom)
#
#
#pickle.dump(all_peaks, open('all_hiccup_peaks.p', 'wb'))

all_peaks = pickle.load(open('all_hiccup_peaks.p', 'rb'))

# Get common between Ddel and Dpn
label1 = 0
label2 = 1

peaks1 = all_peaks[labels[label1]]
peaks2 = all_peaks[labels[label2]]

print(len(peaks1))
print(len(peaks2))


# 10 kb margin
margin = 10000

common_pairs = []
peaks_to_check = []

for (c11, c12) in peaks1:

    for (c21, c22) in peaks2:

        minn1 = abs(c11 - c21)
        minn2 = abs(c12 - c22)

        vall = max(minn1, minn2)

        if vall > margin:
            continue

        common_pairs.append([peaks1[(c11, c12)], peaks2[(c21, c22)]])


        peaks_to_check.append((c11, c12))
        peaks_to_check.append((c21, c22))


# Load PH loops
PH_nearby_loops = pickle.load(open('ddel_dpn_common_PH_maximal.p', 'rb'))
print(len(PH_nearby_loops))

# Unique PH loops between two experiments
unique = []

for pair in PH_nearby_loops:
    cyc1, cyc2 = pair

    cyc1_m = np.argmax(np.abs(cyc1[:,0] - cyc1[:,1])).flatten()
    cyc1_maxx = cyc1[cyc1_m]

    cyc2_m = np.argmax(np.abs(cyc2[:,0] - cyc2[:,1])).flatten()
    cyc2_maxx = cyc2[cyc2_m]

    #print(cyc1_maxx, cyc2_maxx)

    for row in cyc1_maxx:
        if [row[0], row[1]] not in unique:
            unique.append([row[0], row[1]])

    for row in cyc2_maxx:
        if [row[0], row[1]] not in unique:
            unique.append([row[0], row[1]])

    
print(len(unique))

vals = []
genomic_dist = []

long_range = []

for idx, PH_anchor in enumerate(unique):

    print(idx, end='\r')

    g_dist = abs(PH_anchor[0] - PH_anchor[1])

    genomic_dist.append(g_dist)

    if g_dist > 1000:
        long_range.append(PH_anchor)

    vall = math.inf

    for key in peaks_to_check:
        dist = max(abs(PH_anchor[0]*reso - key[0]), abs(PH_anchor[1]*reso - key[1])) 
        vall = min(dist, vall)

    #print(vall)
    vals.append(vall/reso)



plt.scatter(genomic_dist, vals, color='blue', alpha=0.2)
plt.xscale('log', base=10)
plt.yscale('log', base=10)
plt.xlabel('maximal genomic distance in PH loop boundary (units 5kb)')
plt.ylabel('nearest distance (units 5kb)')


plt.show()

plt.cla()
plt.clf()
plt.close()

fig, axs = plt.subplots(2, 2)

long_range = np.array(long_range)

axs[0,1].scatter(long_range[:,0],long_range[:,1], color='blue', alpha=0.2)
axs[0,1].set_xscale('log', base=10)
axs[0,1].set_yscale('log', base=10)

xlims = axs[0,1].get_xlim()
ylims = axs[0,1].get_ylim()

print(xlims, ylims)


#plt.show()
#plt.cla()
#plt.clf()
#plt.close()

# Biomarkers
bio_data_dirr = '../Biomarkers/'
flabel = [\
        'ChIPseq_CTCF'\
        ,'Atacseq'\
        ,'CTCF'
        ]

load_name = bio_data_dirr+flabel[2]+'_reso'+str(reso)+'_parsed.p'

biomarker_dict = pickle.load(open(load_name, 'rb'))


xx = []
yy = []

for binn in range(int(xlims[0]), int(xlims[1])+1):

    if binn in biomarker_dict:
        #print(biomarker_dict[binn])
        #exit()
        for vall in biomarker_dict[binn]:
            xx.append(binn)
            yy.append(float(vall))

axs[1,1].scatter(xx, yy, alpha=0.1, color='blue')

xx = []
yy = []

for binn in range(int(ylims[0]), int(ylims[1])+1):

    if binn in biomarker_dict:
        #print(biomarker_dict[binn])
        #exit()
        for vall in biomarker_dict[binn]:
            xx.append(binn)
            yy.append(float(vall))

axs[0,0].scatter(yy, xx, alpha=0.5, marker='x', color='red')


plt.xscale('log', base=10)
plt.yscale('log', base=10)
plt.show()






    

















