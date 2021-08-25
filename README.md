# Dory
Overcoming barriers to computing Persistent Homology

Repository for https://arxiv.org/abs/2103.05608.


### PyDory
Pydory is a lightweight Python wrapper for Dory. See pydory_getPD.py for example.
Use `pip3 install pydory` to install the Python (v. 3.5+) wrapper for Dory. It will require gcc-10.2 with openMP support. The resulting persistence-pairs are stored as files H0_pers_data.txt, H1_pers_data.txt, and H2_pers_data.txt. The death of features that do not die is stored as -1.


## Reproducing benchmarks

Dory: Run the script Dory_script.sh for timing benchmarks.
Ripser: The benchmarks in https://arxiv.org/abs/2103.05608 were produced using Ripser v 1.2.1 (March, 2021). Use the script ripser_script.sh for timing benchmarks.

Instruments in macOS or Valgrind can be used to determine peak memory usage.

Gudhi (v. 3.4.0) can be installed using `pip3 install Gudhi==3.4.0`. The python file to reproduce benchmarks is Gudhi_benchmarks.py in Benchmarking_codes.

### Data sets
All data sets, except for Hi-C, are provided in the folder Datasets.

### Hi-C data sets
The steps to obtain and process Hi-C data sets are as follows:
1. Download mcool file for Hi-C control from https://data.4dnucleome.org/files-processed/4DNFIFLDVASC/ (~12 GB) and for Hi-C auxin from https://data.4dnucleome.org/files-processed/4DNFILP99QJS/ (~13 GB) to the HiC folder.
2. Run `python3 get_HiC_edges.py` to generate the filtration in sparse format: vertex, vertex, edge-length. This extraction is optimized using Python packages hdf5 and numba.
3. Run `python3 pydory_benchmarks.py` with relevant commands uncommented to get PH for HiC.
