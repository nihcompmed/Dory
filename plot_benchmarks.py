import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.rcParams['text.usetex'] = True


mem = [[159 ,	178 ,	102 ,	186 ,	110]]
mem.append([775 ,	605 ,	597 ,	606 ,	598  ])
mem.append([121 ,	148 ,	141 ,	275 ,	269  ])
mem.append([198 ,	208 ,	200 ,	4860 ,	4850 ])
mem.append([1670 ,	976 ,	954 ,	5600 ,	5600  ])


mem = np.array(mem)


times = [[2.51	,1.54	,1.85	,1.57	,1.92]]

times.append([7.9	,22.75	,28.9	,16.2	,18.15])
times.append([6.89	,4.32	,5.3	,4.4	,5.44])
times.append([14.4	,6.39	,8.93	,8.95	,10.92])
times.append([70.95	,30.54	,49.16	,36.5	,50.87])



dataset = ['dragon(1)', 'fract(2)', 'o3(2)', 'torus4(1)', 'torus4(2)']

times = np.array(times)

labels = [   
'Ripser'  
,'Dory, 4'  
,'Dory, 1'  
,'DoryNS, 4' 
,'DoryNS, 1' ]

markers = [
'P'
,'>'
,'<'
,'^'
,'v']

colors = [
         'tab:blue'
        ,'tab:green'
        ,'tab:orange'
        ,'tab:purple'
        ,'tab:red'
        ]

first = 1

for row in range(5):
    flag = 0
    title = ''
    
    for col in range(5):

        print('dataset', dataset[row], 'code ', labels[col], 'mem ', mem[row, col], 'time ', times[row, col])
    
        if (mem[row,col] != -1):
            plt.scatter(times[row,col], mem[row, col], marker = markers[col], s = 750,
                    label=labels[col], alpha=0.75, color=colors[col])
        else:
            flag = 1
            title += labels[col] + ', '
    
    #if flag:
    #    if title != '':
    #        title = title[:-2]
    #    title += ' cannot compute'
    #else:
    #    title = 'All can compute'

    #plt.xscale('log', base=2)
    #plt.yscale('log', base=2)
    plt.xlabel('Time taken (sec)', fontsize=20)
    plt.ylabel('Memory taken (MB)', fontsize=20)

    if first:
        #title = 'Gudhi cannot compute'
        plt.legend(prop={'size': 18}, loc='lower right', framealpha=1)
        first = 0

    #plt.title(title, fontsize=24)
    #plt.rc('xtick',labelsize=16)
    #plt.rc('ytick',labelsize=16)
    ax = plt.gca()
    ax.tick_params(axis = 'both', which = 'major', labelsize = 24)
    ax.tick_params(axis = 'both', which = 'minor', labelsize = 24)

    #if (row > 2):
    #    plt.yscale('log')
    #    plt.ylabel('Memory taken (MB), log scale', fontsize=20)

    
    plt.tight_layout()
    
    plt.savefig('figures/'+dataset[row]+'_benchmarks.pdf', format='pdf')
    #plt.show()
    plt.cla()
    plt.clf()







