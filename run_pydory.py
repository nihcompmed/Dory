import pydory

###########
# Source
# string
# Accepts comma separated files
###########
source = '.csv'

###########
# filetype
# integer
# 0 : distance matrix
# 1 : point-cloud
# 2 : sparse-format as: point1, point2, edge length
###########
filetype = 

###########
# Threshold for parameter
# float
# Note: For thresh infinite, give larger than maximum pairwise distance
###########
thresh = 

###########
# Number of threads
###########
threads = 

###########
# Dimension
# integer
# 1 : up to and including H1
# 2 : up to and including H2
###########
dim = 


###########
# Target
# string
###########
# The resulting PD will be stored in target folder, appended with the target filename
# For example, if target = 'folder/file_' then following files will be stored in 'folder'
# file_H0_pers_data.txt, file_H1_pers_data.txt, ....
###########
target = ''


###########
# Final command
# The order of parameters has to be maintained as shown
###########
pydory.compute_PH(source, thresh, filetype, threads, target, dim)

