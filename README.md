# Dory
Computing Persistent Homology (VR-filtration) and tight cycle representatives (loops and voids) for large data sets

**Please make sure that the latest version 1.1.109 of pyDory is installed.**

Repository for https://www.sciencedirect.com/science/article/pii/S1877750324000838 and https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1010341.


### PyDory
PyDory is a lightweight Python wrapper for Dory.

Use `pip3 install pydory` to install the Python (v. 3.5+) wrapper for Dory. It will require gcc-10.2 with openMP support. If gcc-10.2 is not the default gcc version, then `CC=gcc-10.2 python3 -m pip install pydory` should work.

**See pyDory_getPD.py for an example of a python script with details on usage.**
It is currently not compatible with Jupyter notebook, and should be run via command prompt or terminal.

The resulting persistence-pairs are stored as files H0_pers_data.txt, H1_pers_data.txt, and H2_pers_data.txt. The death of features that do not die is stored as -1.

CAUTION: Please make sure floating-point numbers are  stored in decimal notation and not scientific notation in the source files. 


## Reproducing benchmarks

Dory: Run the script Dory_script.sh for timing benchmarks.
Ripser: Use the script ripser_script.sh for timing benchmarks.

Instruments in macOS or Valgrind can be used to determine peak memory usage.

Gudhi (v. 3.4.0) can be installed using `pip3 install Gudhi==3.4.0`. The python file to reproduce benchmarks is Gudhi_benchmarks.py in Benchmarking_codes.

### Data sets
All data sets, except for Hi-C, are provided in the folder Datasets.

### Hi-C data sets
The steps to obtain and process Hi-C data sets are as follows:
1. Download mcool file for Hi-C control from https://data.4dnucleome.org/files-processed/4DNFIFLDVASC/ (~12 GB) and for Hi-C auxin from https://data.4dnucleome.org/files-processed/4DNFILP99QJS/ (~13 GB) to the Datasets/HiC folder.
2. Run `python3 get_HiC_edges.py` to generate the filtration in sparse format: vertex, vertex, edge-length. This extraction is optimized using Python packages hdf5 and numba.
3. Run `python3 pydory_getPD.py` with relevant commands uncommented to get PH for HiC.
