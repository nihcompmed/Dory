import numpy as np
import pickle

ff = open('ENCFF159KBI.gtf', 'r')

reso = 5000

all_info = dict()

all_info_bins = dict()
# [chromosome][typ][name][bins:[], desc:'']

count = 0

for line in ff:
    print(count, end='\r')
    count += 1
    #line = line.split('\t')
    #print(line)
    line = line.split(';')
    #print(line)

    info = line[0].split('\t')


    chrom = info[0]
    typ = info[2]
    bin1 = int(int(info[3])/reso)
    bin2 = int(int(info[4])/reso)
    binns = frozenset(list(range(bin1, bin2+1)))

    remain = line[1:]

    for field in remain:
        field = field.strip(' ')
        field = field.split(' ')
        if field[0] == 'gene_name':
            gene_name = field[1]
        elif field[0] == 'gene_type':
            gene_typ_desc  = field[1]
    
    #print('typ desc', gene_typ_desc, 'name', gene_name)

    if chrom not in all_info:
        all_info[chrom] = dict()

    if typ not in all_info[chrom]:
        all_info[chrom][typ] = dict()

    if gene_name not in all_info[chrom][typ]:
        all_info[chrom][typ][gene_name] = {'binns':binns, 'desc':gene_typ_desc}
    else:
        all_info[chrom][typ][gene_name]['binns'] =\
                        all_info[chrom][typ][gene_name]['binns'].union(binns)

    if chrom not in all_info_bins:
        all_info_bins[chrom] = dict()

    for bb in range(bin1, bin2+1):
        if bb not in all_info_bins[chrom]:
            all_info_bins[chrom][bb] = []

        all_info_bins[chrom][bb].append(gene_name)


# types
#gene
#transcript
#exon
#CDS
#start_codon
#stop_codon
#UTR
#Selenocysteine

pickle.dump(all_info, open('all_gene_to_bins_info_hg38_reso'+str(reso)+'.p', 'wb'))

pickle.dump(all_info_bins, open('all_bins_to_gene_info_hg38_reso'+str(reso)+'.p', 'wb'))

#binn_gene_count = dict()
#
#for gene in all_info['chr1']['gene']:
#    #print(gene, all_info['chr1']['gene'][gene]['binns'])
#    #input('w')
#
#    for binn in all_info['chr1']['gene'][gene]['binns']:
#        if binn not in binn_gene_count:
#            binn_gene_count[binn] = 1
#        else:
#            binn_gene_count[binn] += 1
#
#print(binn_gene_count)






