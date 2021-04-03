# Dory
Overcoming barriers to computing Persistent Homology

Repository for https://arxiv.org/abs/2103.05608.


### Python installation for Dory
Use `pip3 install pydory` to install the Python (v. 3.5+) wrapper for Dory. It will require gcc-10.2 with openMP support. The resulting persistence-pairs are stored as files H0_pers_data.txt, H1_pers_data.txt, and H2_pers_data.txt. The features that do not die have death -1.

## Reproducing benchmarks

Dory: Run the python file pydory_benchmarks.py for benchmarks. This wrapper has little overhead, and peak memory usage and computation time is similar to the C source code.

Ripser: The benchmarks in https://arxiv.org/abs/2103.05608 were produced using Ripser as on March 22, 2021. It has since been updated. We provide a copy of the legacy version (copied from `https://github.com/Ripser/ripser` on March 22) in this repository. We used `c++ -std=c++11 ripser.cpp -o ripser -O3` to compile Ripser.

Gudhi (v. 3.4.0) can be installed using `pip3 install Gudhi==3.4.0`. The python file to reproduce benchmarks is Gudhi_benchmarks.py in Benchmarking_codes.

### Data sets
All data sets, except for Hi-C, are provided in the folder Datasets.

### Hi-C data sets
The steps to obtain and process Hi-C data sets are as follows:
1. Download mcool file for Hi-C control from https://data.4dnucleome.org/files-processed/4DNFIFLDVASC/ (~12 GB) and for Hi-C auxin from https://data.4dnucleome.org/files-processed/4DNFILP99QJS/ (~13 GB) to the HiC folder.
2. Run `python3 get_edges.py` to generate the filtration in sparse format: vertex, vertex, edge-length. This extraction is optimized using Python packages hdf5 and numba.
3. Run `python3 pydory_benchmarks.py` with relevant commands uncommented to get PH for HiC.
