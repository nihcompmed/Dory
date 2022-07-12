import numpy as np
import cooler
import h5py
import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn as sns
from numba import njit
import pickle
import networkx as nx
from matplotlib import cm

cmap = mpl.cm.get_cmap('gnuplot_r')


def get_chrom(binn, chrom_start_bins_list):

    for chrom, start in chrom_start_bins_list:

        if start >= binn:
            break
        else:
            this_chrom = chrom

    return this_chrom

def get_info(exp_name, dim):

    cfile = exp_name+'.cool'
    
    cf = cooler.Cooler(cfile)
    
    chrom_start_bins = {}
    chrom_start_bins_list = []
    
    for ch in cf.chromnames:
        beg, end = cf.extent(ch)
        chrom_start_bins[ch] = beg
        chrom_start_bins_list.append([ch, beg])
    
    chrom_start_bins_list.append([ 'last_bin', cf.extent(cf.chromnames[-1])[1] ])
    
    # Keep it simple
    # Initiate a chorm dict (except last bin)
    chrom_dict = dict()
    for [chrom, binn] in chrom_start_bins_list[:-1]:
        chrom_dict[chrom] = {'cis':dict(), 'trans':dict()}
    
    # Go over birth cycles
    
    smoothened = exp_name + '/smooth_H1_estimated_pers.txt' 
    
    ff = open(smoothened, 'r')
    
    all_info = []
    
    cis_lengths = []
    trans_lengths = []
    
    
    for line in ff:
        if line[0] != 's':
            continue

        line = line.split(',')
        line = line[1:-1]
        line = [int(x) for x in line]
    
        lenn = int(len(line))

        # Ignore degenerate cycles
        # REDUNDANT NOW (since code that estimates persistence takes care of it)
        if lenn < 3:
            continue
    
        ## We consider cycles of length at least 5
        #if lenn < 5:
        #    continue
    
        nn = 0
    
        chroms_in_cycle = frozenset()
    
        while nn < len(line):
            bin1 = line[nn] 
            bin2 = line[(nn+1)%lenn] 
    
            bin1_chrom = get_chrom(bin1, chrom_start_bins_list)
            bin2_chrom = get_chrom(bin2, chrom_start_bins_list)
    
            chroms_in_cycle = chroms_in_cycle.union(frozenset({bin1_chrom, bin2_chrom}))
    
            nn += 1
    
    
        all_info.append([chroms_in_cycle, lenn])
    
        if len(chroms_in_cycle) == 1:
            cis_lengths.append(lenn)
            chrom = list(chroms_in_cycle)[0]
            if lenn not in chrom_dict[chrom]['cis']:
                chrom_dict[chrom]['cis'][lenn] = 1
            else:
                chrom_dict[chrom]['cis'][lenn] += 1
        else:
            trans_lengths.append(lenn)
            for chrom in list(chroms_in_cycle):
                if lenn not in chrom_dict[chrom]['trans']:
                    chrom_dict[chrom]['trans'][lenn] = 1
                else:
                    chrom_dict[chrom]['trans'][lenn] += 1


    return all_info, chrom_dict, cis_lengths, trans_lengths

chrom_list = [\
'chr1',\
'chr2',\
'chr3',\
'chr4',\
'chr5',\
'chr6',\
'chr7',\
'chr8',\
'chr9',\
'chr10',\
'chr11',\
'chr12',\
'chr13',\
'chr14',\
'chr15',\
'chr16',\
'chr17',\
'chr18',\
'chr19',\
'chr20',\
'chr21',\
'chr22',\
'chrX',\
'chrY',\
]

chrom_markers = {\
         'chr1':'o'\
        ,'chr2':'o'\
        ,'chr3':'o'\
        ,'chr4':'o'\
        ,'chr5':'o'\
        ,'chr6':'o'\
        ,'chr7':'o'\
        ,'chr8':'o'\
        ,'chr9':'o'\
        ,'chr10':'o'\
        ,'chr11':'o'\
        ,'chr12':'o'\
        ,'chr13':'o'\
        ,'chr14':'o'\
        ,'chr15':'o'\
        ,'chr16':'o'\
        ,'chr17':'o'\
        ,'chr18':'o'\
        ,'chr19':'o'\
        ,'chr20':'o'\
        ,'chr21':'o'\
        ,'chr22':'o'\
        ,'chrX':'$X$'\
        ,'chrY':'$Y$'\
        }
chrom_alpha = {\
         'chr1':0.5\
        ,'chr2':0.5\
        ,'chr3':0.5\
        ,'chr4':0.5\
        ,'chr5':0.5\
        ,'chr6':0.5\
        ,'chr7':0.5\
        ,'chr8':0.5\
        ,'chr9':0.5\
        ,'chr10':0.5\
        ,'chr11':0.5\
        ,'chr12':0.5\
        ,'chr13':0.5\
        ,'chr14':0.5\
        ,'chr15':0.5\
        ,'chr16':0.5\
        ,'chr17':0.5\
        ,'chr18':0.5\
        ,'chr19':0.5\
        ,'chr20':0.5\
        ,'chr21':0.5\
        ,'chr22':0.5\
        ,'chrX':0.9\
        ,'chrY':0.9\
        }


dim = 1

# with auxin
exp_name = '4DNFI2X45Z5L'

all_info, chrom_dict_auxin, cis_lengths, trans_lengths = get_info(exp_name, dim)

# control
exp_name = '4DNFITH9KJF1'

all_info, chrom_dict_control, cis_lengths, trans_lengths = get_info(exp_name, dim)

print(len(cis_lengths), len(trans_lengths))


def get_cis_trans_lens(data):

    cis_data = data['cis']
    trans_data = data['trans']

    cis_lenns = []
    cis_counts = []
    for lenn in cis_data:
        cis_lenns.append(lenn)
        cis_counts.append(cis_data[lenn])

    trans_lenns = []
    trans_counts = []
    for lenn in trans_data:
        trans_lenns.append(lenn)
        trans_counts.append(trans_data[lenn])

    return cis_lenns, cis_counts, trans_lenns, trans_counts


global_maxx = 0

#fig, axs = plt.subplots(1, 2, sharey=True)

#plt.figure(1)

for chrom_idx, chrom in enumerate(chrom_list):

    print('Doing', chrom)

    data = chrom_dict_control[chrom]

    control_cis_lenns, control_cis_counts, control_trans_lenns, control_trans_counts = get_cis_trans_lens(data)


    maxx_lenn = max(max(control_cis_lenns), max(control_trans_lenns))
    
    data = chrom_dict_auxin[chrom]

    auxin_cis_lenns, auxin_cis_counts, auxin_trans_lenns, auxin_trans_counts = get_cis_trans_lens(data)

    maxx_lenn = max(maxx_lenn, max(auxin_cis_lenns), max(auxin_trans_lenns))

    global_maxx = max(maxx_lenn, global_maxx)


    control_cis_array = np.zeros((maxx_lenn+1,), dtype=int)
    auxin_cis_array = np.zeros((maxx_lenn+1,), dtype=int)

    control_trans_array = np.zeros((maxx_lenn+1,), dtype=int)
    auxin_trans_array = np.zeros((maxx_lenn+1,), dtype=int)

    for idx, lenn in enumerate(control_cis_lenns):

        val = control_cis_counts[idx]

        control_cis_array[lenn] = val

    for idx, lenn in enumerate(auxin_cis_lenns):

        val = auxin_cis_counts[idx]

        auxin_cis_array[lenn] = val

    for idx, lenn in enumerate(control_trans_lenns):

        val = control_trans_counts[idx]

        control_trans_array[lenn] = val

    for idx, lenn in enumerate(auxin_trans_lenns):

        val = auxin_trans_counts[idx]

        auxin_trans_array[lenn] = val


    zero_control_cis = frozenset(np.argwhere(control_cis_array==0).flatten())
    zero_auxin_cis = frozenset(np.argwhere(auxin_cis_array==0).flatten())

    common = zero_control_cis.intersection(zero_auxin_cis)

    plot_idxs = list(frozenset(list(range(4, maxx_lenn+1))) - common)

    cis_diff = control_cis_array - auxin_cis_array

    plt.figure(1)
    plt.scatter(plot_idxs, cis_diff[plot_idxs], marker=chrom_markers[chrom]\
            , color = cmap(chrom_idx/len(chrom_list)), alpha=chrom_alpha[chrom])


    zero_control_trans = frozenset(np.argwhere(control_trans_array==0).flatten())
    zero_auxin_trans = frozenset(np.argwhere(auxin_trans_array==0).flatten())

    common = zero_control_trans.intersection(zero_auxin_trans)

    plot_idxs = list(frozenset(list(range(4, maxx_lenn+1))) - common)

    trans_diff = control_trans_array - auxin_trans_array

    plt.figure(2)
    plt.scatter(plot_idxs, trans_diff[plot_idxs], marker=chrom_markers[chrom]\
            , color = cmap(chrom_idx/len(chrom_list)), alpha=chrom_alpha[chrom])

    #plt.show()

    #exit()


#print(global_maxx)
#exit()



plt.figure(1)
plt.hlines(0, 4, global_maxx+1, colors='black')
plt.xlabel('lengths', fontsize=14)
plt.ylabel('diff in counts (control - auxin)', fontsize=14)
plt.title('cis')
plt.ylim([-2000, 2000])
plt.xscale('symlog', base=10)
plt.yscale('symlog', base=10)
#fig = plt.gcf()
#ticks = list(np.array(list(range(len(chrom_list))))/len(chrom_list))
#cbar = fig.colorbar(cm.ScalarMappable(cmap=cmap), ax=plt.gca(), ticks=ticks)
#cbar.ax.set_yticklabels(chrom_list)
plt.savefig('figures/HiC_cis_diff_lens.pdf', dpi=600)
#plt.show()

plt.figure(2)
plt.hlines(0, 4, global_maxx+1, colors='black')
plt.xlabel('lengths', fontsize=14)
plt.ylim([-2000, 2000])
plt.xscale('symlog', base=10)
plt.yscale('symlog', base=10)
plt.title('trans')
plt.savefig('figures/HiC_trans_diff_lens.pdf', dpi=600)


