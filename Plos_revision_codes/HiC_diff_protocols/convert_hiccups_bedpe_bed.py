import numpy as np


hiccup_files = [\
         '../fa_dsg_ddel_hiccups'\
        ,'../fa_dsg_dpn_hiccups'\
        ,'../fa_dsg_mnase_hiccups'\
        ,'../fa_dsg_ddel_dpn_hiccups'\
        ,'../fa_dpn_hiccups'\
                ]


for ff in hiccup_files:

    bedpe_fname = ff + '.bedpe'

    bed1_fname = ff+'_left.bed'
    bed2_fname = ff+'_right.bed'

    bedpe = open(bedpe_fname, 'r')

    bed1 = open(bed1_fname, 'w')
    bed2 = open(bed2_fname, 'w')

    for line in bedpe:

        line = line.split('\t')

        chrom = line[0]
        if chrom != 'chr1':
            continue

        start1 = int(line[1])
        end1 = int(line[2])

        start2 = int(line[4])
        end2 = int(line[5])

        pos1 = (start1 + end1)/2
        pos2 = (start2 + end2)/2

        start1 = pos1 - 5000
        end1 = pos1 + 5000
        
        start2 = pos2 - 5000
        end2 = pos2 + 5000

        bed1.write('chr1' + '\t' + str(int(start1)) + '\t' + str(int(end1)) + '\n')
        bed2.write('chr1' + '\t' + str(int(start2)) + '\t' + str(int(end2)) + '\n')

    bed1.close()
    bed2.close()









