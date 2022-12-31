import pickle

data_dirr = '../Biomarkers/'


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

reso = 5000

for label_idx in range(len(fnames)):

    ff = open(fnames[label_idx], 'r')
    
    unique_bin_dist = []
    
    bin_signal_dict = dict()
    
    count = 0
    
    for line in ff:
        line = line.split('\t')
    
        chrom = line[0]
    
        if chrom != 'chr1':
            continue
    
        count += 1
    
        start = int(int(line[1])/reso)
        end = int(int(line[2])/reso)
        signalval = line[-4]

        peak = line[-1]
    
        if start not in bin_signal_dict:
            bin_signal_dict[start] = [signalval]
        else:
            bin_signal_dict[start].append(signalval)
    
        if end not in bin_signal_dict:
            bin_signal_dict[end] = [signalval]
        else:
            bin_signal_dict[end].append(signalval)
    
        bin_dist = end-start
    
        if bin_dist not in unique_bin_dist:
            unique_bin_dist.append(bin_dist)
    
    
    ff.close()
    
    uni_lens = []
    for key in bin_signal_dict:
        if len(bin_signal_dict[key]) not in uni_lens:
            uni_lens.append(len(bin_signal_dict[key]))
    
    save_name = data_dirr+flabel[label_idx]+'_reso'+str(reso)+'_parsed.p'
    
    pickle.dump(bin_signal_dict, open(save_name, 'wb'))




