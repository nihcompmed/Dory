import numpy as np
import pickle

ff = open('../Biomarkers/ENCFF435PVY.bed')

ccre_dict = dict()

reso = 5000

ccre_types = []

for line in ff:
    line = line.split('\t')

    chrom = line[0]

    start = int(line[1])
    end = int(line[2])

    ccre = line[-2]

    #if ccre not in ccre_types:
    #    ccre_types.append(ccre)

    #print(chrom, start, end, ccre)
    #exit()

    start_scale = int(start/reso)
    end_scale = int(end/reso)

    if chrom not in ccre_dict:
        ccre_dict[chrom] = dict()

    if (start_scale, end_scale) not in ccre_dict[chrom]:
        ccre_dict[chrom][(start_scale, end_scale)] = [ccre]
    else:
        if ccre not in ccre_dict[chrom][(start_scale, end_scale)]:
            ccre_dict[chrom][(start_scale, end_scale)].append(ccre)

    #print(ccre_dict)

    #input('w')

ff.close()

pickle.dump(ccre_dict, open('../Biomarkers/ccre_reso'+str(reso)+'.p', 'wb'))





