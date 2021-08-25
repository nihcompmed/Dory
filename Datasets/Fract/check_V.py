import numpy as np

redV = np.loadtxt('redV1coH2_V_data.txt')
noredV = np.loadtxt('noredV1coH2_V_data.txt')

print(np.amax(noredV - redV))
