import pyBigWig


fnames = [\
        '../Biomarkers/GSE165895_HFFc6_H3K4me3.bigWig'\
        ]

flabel = [\
        'H3k4me3'\
        ]

reso = 5000

ff = pyBigWig.open(fnames[0])

print(ff.header())

#ff = open(fnames[0], 'r')
#
#for line in ff:
#    line = line.split('\t')
#    print(line)
#    exit()
#
#ff.close()




