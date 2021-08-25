Folder for HiC datasets and sparse edge files.

Download mcool file for Hi-C control from https://data.4dnucleome.org/files-processed/4DNFIFLDVASC/ (~12 GB) and for Hi-C auxin from https://data.4dnucleome.org/files-processed/4DNFILP99QJS/ (~13 GB) to the HiC folder.

Run python3 get_edges.py to generate the filtration in sparse format: vertex, vertex, edge-length. This extraction is optimized using Python packages hdf5 and numba.
