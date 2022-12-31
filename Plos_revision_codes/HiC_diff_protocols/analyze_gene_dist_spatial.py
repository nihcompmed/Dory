import matplotlib.pyplot as plt
import numpy as np


reso = 5000
#reso = 10000

edge_file_name = [\
         '../fa_dsg_ddel_chr1_KR_reso'+str(reso)+'.csv'\
        ,'../fa_dsg_dpn_chr1_KR_reso'+str(reso)+'.csv'\
        ,'../fa_dsn_mnase_chr1_KR_reso'+str(reso)+'.csv'\
        ,'../fa_dsg_ddel_dpn_chr1_KR_reso'+str(reso)+'.csv'\
        ,'../fa_dpn_chr1_KR_reso'+str(reso)+'.csv'\
        ]
plt_fname = 'gene_dist_vs_spatial_ptiles.pdf'

# For observed/expected
edge_file_name = [\
         '../fa_dsg_ddel_chr1_KR_oe_reso'+str(reso)+'.csv'\
        ,'../fa_dsg_dpn_chr1_KR_oe_reso'+str(reso)+'.csv'\
        ,'../fa_dsn_mnase_chr1_KR_oe_reso'+str(reso)+'.csv'\
        ,'../fa_dsg_ddel_dpn_chr1_KR_oe_reso'+str(reso)+'.csv'\
        ,'../fa_dpn_chr1_KR_oe_reso'+str(reso)+'.csv'\
        ]
plt_fname = 'gene_dist_vs_spatial_ptiles_oe.pdf'

labels = [\
        'FA+DSG, Ddel'\
        ,'FA+DSG, Dpn'\
        ,'FA+DSG, MNase'\
        ,'FA+DSG, Ddel-Dpn'\
        ,'FA, Dpn'\
        ]


g_dists = list(range(1, 6))

plt_g_dists = np.array(g_dists)*reso/1000

ptiles = [10, 50, 90]

fig, axs = plt.subplots(1, 3, figsize=(12, 8))


for e_idx, edge_file in enumerate(edge_file_name):

    print(edge_file)

    data = np.loadtxt(edge_file, delimiter=',')


    bin_dist = np.abs(data[:,0] - data[:,1])

    vals = data[:,2]

    this_data = []


    for g in g_dists:

        idxs = np.argwhere(bin_dist==g).flatten()

        vv = vals[idxs]

        this_ptiles = np.percentile(vv, ptiles)

        this_data.append(this_ptiles)

    this_data = np.array(this_data)

    for p_idx, pt in enumerate(ptiles):

        axs[p_idx].plot(plt_g_dists, this_data[:, p_idx], label=labels[e_idx])

    print(this_data[:6])

axs[0].legend()
axs[1].legend()
axs[2].legend()

axs[0].set_xlabel('genomic distance (kb)')
axs[0].set_ylabel('spatial estimate')

axs[0].set_title('10%-tile')
axs[1].set_title('50%-tile')
axs[2].set_title('90%-tile')



#plt.show()

plt.savefig(plt_fname, dpi=600)


##########################
# OBSERVED
##########################

##########################
# Reso 5000
##########################
#../fa_dsg_ddel_chr1_KR.csv
#[[0.00127402 0.00165758 0.00240757]
# [0.00275244 0.00375676 0.00534985]
# [0.00440118 0.00630977 0.00879996]
# [0.00610019 0.00888489 0.01237322]
# [0.00793804 0.01140681 0.01622359]]
#../fa_dsg_dpn_chr1_KR.csv
#[[0.00144694 0.00198671 0.00276559]
# [0.00278329 0.00388867 0.00545233]
# [0.00436364 0.00616116 0.00868856]
# [0.00597983 0.0084514  0.01214002]
# [0.00759132 0.0107644  0.01561109]]
#../fa_dsn_mnase_chr1_KR.csv
#[[0.00047721 0.00057735 0.00073935]
# [0.00254124 0.00323818 0.00449164]
# [0.00442511 0.00582752 0.00836386]
# [0.00638778 0.00852553 0.01257781]
# [0.00832092 0.0113443  0.01703869]]
#../fa_dsg_ddel_dpn_chr1_KR.csv
#[[0.00101757 0.00134216 0.00185776]
# [0.00255027 0.00346039 0.00480261]
# [0.00432124 0.00593789 0.00835398]
# [0.00613818 0.00845516 0.01211395]
# [0.0079926  0.01099161 0.01615373]]
#../fa_dpn_chr1_KR.csv
#[[0.00340766 0.00488648 0.00710244]
# [0.00538282 0.00775345 0.01101183]
# [0.00737594 0.01058858 0.01509577]
# [0.00943415 0.01349069 0.01939283]
# [0.01146682 0.01637872 0.02395455]]
#

##########################
# Reso 10000
##########################
#../fa_dsg_ddel_chr1_KR_reso10000.csv
#[[0.00065091 0.00088101 0.00119952]
# [0.00159729 0.00223261 0.00303552]
# [0.00252335 0.00353491 0.00497419]
# [0.00336814 0.00480362 0.00700289]
# [0.00413446 0.00600674 0.00914071]]
#../fa_dsg_dpn_chr1_KR_reso10000.csv
#[[0.00062385 0.000855   0.00117106]
# [0.00147386 0.00206366 0.00283449]
# [0.00232408 0.00325093 0.00457273]
# [0.0031051  0.0044004  0.00640587]
# [0.00381453 0.00550143 0.0083232 ]]
#../fa_dsn_mnase_chr1_KR_reso10000.csv
#[[0.00035919 0.0004427  0.00059415]
# [0.00173559 0.00227874 0.00314951]
# [0.0028835  0.00386612 0.00558935]
# [0.00390811 0.00540546 0.00817287]
# [0.00477677 0.00690492 0.01087386]]
#../fa_dsg_ddel_dpn_chr1_KR_reso10000.csv
#[[0.00050855 0.00067343 0.00091857]
# [0.0015022  0.00203807 0.00278725]
# [0.0024731  0.00335243 0.00474511]
# [0.00335806 0.00462213 0.00686027]
# [0.00415534 0.00585644 0.00906259]]
#../fa_dpn_chr1_KR_reso10000.csv
#[[0.00127328 0.00180698 0.00249192]
# [0.00236419 0.0033213  0.00455864]
# [0.00343167 0.00480044 0.00669773]
# [0.00441146 0.00622183 0.00891011]
# [0.00529317 0.00755307 0.01120925]]



###########################
# Observed/expected
#oe
###########################

Resolution 5000
../fa_dsg_ddel_chr1_KR_oe_reso5000.csv
[[0.64966384 0.88465305 1.21655353]
 [0.64357308 0.89264337 1.2398689 ]
 [0.63693453 0.89746528 1.25687819]
 [0.63362493 0.89735147 1.27764433]
 [0.6287245  0.8937422  1.30856971]]
../fa_dsg_dpn_chr1_KR_oe_reso5000.csv
[[0.64210054 0.88163416 1.22727297]
 [0.63770099 0.8909626  1.24922551]
 [0.63351628 0.89448111 1.26140975]
 [0.63128125 0.89220118 1.28160309]
 [0.62855256 0.89127982 1.29257944]]















