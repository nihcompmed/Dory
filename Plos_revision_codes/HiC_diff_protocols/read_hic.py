import hicstraw
import matplotlib.pyplot as plt
import numpy as np

hic_files = [\
        "../4DNFII84FBKM.hic"\
        ,"../4DNFIMROE6N4.hic"\
        ,"../4DNFIPC7P27B.hic"\
        ,"../4DNFI5SUHCGZ.hic"\
        ,"../4DNFIFLJLIS5.hic"\
        ]

edge_file_name = [\
         '../fa_dsg_ddel_chr1_KR'\
        ,'../fa_dsg_dpn_chr1_KR'\
        ,'../fa_dsn_mnase_chr1_KR'\
        ,'../fa_dsg_ddel_dpn_chr1_KR'\
        ,'../fa_dpn_chr1_KR'\
        ]

## For ebserved/expected
#edge_file_name = [\
#         '../fa_dsg_ddel_chr1_KR_oe'\
#        ,'../fa_dsg_dpn_chr1_KR_oe'\
#        ,'../fa_dsn_mnase_chr1_KR_oe'\
#        ,'../fa_dsg_ddel_dpn_chr1_KR_oe'\
#        ,'../fa_dpn_chr1_KR_oe'\
#        ]



#hic = hicstraw.HiCFile("../4DNFII84FBKM.hic")

#print(hic)


#print(hic.getChromosomes())
#print(hic.getGenomeID())
#print(hic.getResolutions())

#normalizations = ['NONE'\
#                , 'KR'\
#                , 'VC'\
#                ]

#meta = hicstraw.read_metadata("../4DNFII84FBKM.hic")

#print(meta)

#for item in hic:
#    print(item)

reso = 5000

#for norm in normalizations:
#    print('loading', norm)
#    try:
#        result_KR = hicstraw.straw('observed', norm, hic_files[0], '1', '1', 'BP', reso)
#    except:
#        print('not found')
#
#
#
#exit()

#mzd = hic.getMatrixZoomData('1', '1', "observed", "VC", "BP", 5000)

#reso = 10000

for idx, hic_file in enumerate(hic_files):

    result_KR = hicstraw.straw('observed', 'KR', hic_file, '1', '1', 'BP', reso)

    ## For observed/expected
    #result_KR = hicstraw.straw('oe', 'KR', hic_file, '1', '1', 'BP', reso)
    
    lenn = len(result_KR)

    print(hic_file)
    
    
    #bin_dist = np.zeros((lenn,))
    #arr_KR = np.zeros((lenn,))

    ff = open(edge_file_name[idx]+'_reso'+str(reso)+'.csv', 'w')
    
    for i in range(lenn):
    
        if i%100000 == 0:
            print(i, end='\r')

        bin1 = result_KR[i].binX/reso
        bin2 = result_KR[i].binY/reso
        val = 1/result_KR[i].counts
    
        if bin1 < bin2:
            ff.write(str(int(bin1))+','+str(int(bin2))+','+str(val)+'\n')
        elif bin1 > bin2:
            ff.write(str(int(bin2))+','+str(int(bin1))+','+str(val)+'\n')
    

    ff.close()


    
    #g_dist = list(range(1,11))
    #
    #ptiles = [10, 50, 90]
    #
    #results = []
    #
    #for gg in g_dist:
    #
    #    idxs = np.argwhere(bin_dist==gg).flatten()
    #
    #    vals = 1/arr_KR[idxs]
    #
    #    ptile_vals = np.percentile(vals, ptiles)
    #    results.append(ptile_vals)
    #
    #
    #results = np.array(results)
    #
    #
    #for p_idx in range(len(ptiles)):
    #
    #    plt.plot(g_dist, results[:, p_idx], label=str(ptiles[p_idx])+'%-tile')
    
    
    
    
    
    
