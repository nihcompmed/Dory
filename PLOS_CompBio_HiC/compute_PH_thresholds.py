import pydory as dory
import check_connected_helper as cch
import helper_functions as hf 
import networkx as nx
import matplotlib.pyplot as plt
import time


# Select experiment
exp_name = 'control'
#exp_name = 'auxin'

source = exp_name+'/edges.csv'
target = exp_name+'/'

filetype = 2

dim = 1
threads = 4
thresh = 150

print('Doing thresh', thresh)

thresh_target = target+'diff_thresholds/thresh_'+str(thresh)+'_'
start_time = time.time()
dory.compute_PH(source, 0, thresh, filetype, threads, thresh_target, dim, 1, 1, thresh, 1)
end_time = time.time()

print('time taken', round(end_time - start_time,2))


