#!/bin/sh
echo "Compiling Dory"
gcc-10.2 Dory/Dory.c -O3 -fopenmp -pthread -o Dory/Dory.o 
echo "Compiling DoryNS (Comb idx)"
gcc-10.2 Dory/Dory.c -O3 -fopenmp -pthread -D COMBIDX -o Dory/DoryNS.o 
echo "Dragon, Dory 4 threads"
time Dory/Dory.o Datasets/Dragon/dragon2000_locs.csv 1000 1 4 Datasets/Dragon/Dory_th4_ 1 0 0 1
echo "Dragon, Dory 1 thread"
time Dory/Dory.o Datasets/Dragon/dragon2000_locs.csv 1000 1 1 Datasets/Dragon/Dory_th1_ 1 0 0 1
echo "Fract, Dory 4 threads"
time Dory/Dory.o Datasets/Fract/fractal_r_distmat.csv 1000 0 4 Datasets/Fract/Dory_th4_ 2 0 0 1
echo "Fract, Dory 1 thread"
time Dory/Dory.o Datasets/Fract/fractal_r_distmat.csv 1000 0 1 Datasets/Fract/Dory_th1_ 2 0 0 1
echo "o3, Dory 4 threads"
time Dory/Dory.o Datasets/o3/o3_8192.csv 1 1 4 Datasets/o3/Dory_th4 2 0 0 1
echo "o3, Dory 1 thread"
time Dory/Dory.o Datasets/o3/o3_8192.csv 1 1 1 Datasets/o3/Dory_th1 2 0 0 1
echo "torus4 (H1), Dory 4 threads"
time Dory/Dory.o Datasets/torus4/torus4_locs.csv 0.15 1 4 Datasets/torus4/Dory_H1_th4_ 1 0 0 1
echo "torus4 (H1), Dory 1 thread"
time Dory/Dory.o Datasets/torus4/torus4_locs.csv 0.15 1 1 Datasets/torus4/Dory_H1_th1_ 1 0 0 1
echo "torus4 (H2), Dory 4 threads"
time Dory/Dory.o Datasets/torus4/torus4_locs.csv 0.15 1 4 Datasets/torus4/Dory_H2_th4_ 2 0 0 1
echo "torus4 (H2), Dory 1 thread"
time Dory/Dory.o Datasets/torus4/torus4_locs.csv 0.15 1 1 Datasets/torus4/Dory_H2_th1_ 2 0 0 1
echo "HiC (control), Dory 4 threads"
time Dory/Dory.o Datasets/HiC/control_500.csv 400 2 4 Datasets/HiC/Dory_control_th4_ 2 0 0 1
echo "HiC (control), Dory 1 thread"
time Dory/Dory.o Datasets/HiC/control_500.csv 400 2 1 Datasets/HiC/Dory_control_th1_ 2 0 0 1
echo "HiC (auxin), Dory 4 threads"
time Dory/Dory.o Datasets/HiC/auxin_500.csv 400 2 4 Datasets/HiC/Dory_auxin_th4_ 2 0 0 1
echo "HiC (auxin), Dory 1 thread"
time Dory/Dory.o Datasets/HiC/auxin_500.csv 400 2 1 Datasets/HiC/Dory_auxin_th1_ 2 0 0 1
echo "Dragon, DoryNS 4 threads"
time Dory/DoryNS.o Datasets/Dragon/dragon2000_locs.csv 1000 1 4 Datasets/Dragon/DoryNS_th4_ 1 0 0 1
echo "Dragon, DoryNS 1 thread"
time Dory/DoryNS.o Datasets/Dragon/dragon2000_locs.csv 1000 1 1 Datasets/Dragon/DoryNS_th1_ 1 0 0 1
echo "Fract, DoryNS 4 threads"
time Dory/DoryNS.o Datasets/Fract/fractal_r_distmat.csv 1000 0 4 Datasets/Fract/DoryNS_th4_ 2 0 0 1
echo "Fract, DoryNS 1 thread"
time Dory/DoryNS.o Datasets/Fract/fractal_r_distmat.csv 1000 0 1 /Datasets/Fract/DoryNS_th1_ 2 0 0 1
echo "o3, DoryNS 4 threads"
time Dory/DoryNS.o Datasets/o3/o3_8192.csv 1 1 4 Datasets/o3/DoryNS_th4 2 0 0 1
echo "o3, DoryNS 1 thread"
time Dory/DoryNS.o Datasets/o3/o3_8192.csv 1 1 1 Datasets/o3/DoryNS_th1 2 0 0 1
echo "torus4 (H1), DoryNS 4 threads"
time Dory/DoryNS.o Datasets/torus4/torus4_locs.csv 0.15 1 4 Datasets/torus4/DoryNS_H1_th4_ 1 0 0 1
echo "torus4 (H1), DoryNS 1 thread"
time Dory/DoryNS.o Datasets/torus4/torus4_locs.csv 0.15 1 1 Datasets/torus4/DoryNS_H1_th1_ 1 0 0 1
echo "torus4 (H2), DoryNS 4 threads"
time Dory/DoryNS.o Datasets/torus4/torus4_locs.csv 0.15 1 4 Datasets/torus4/DoryNS_H2_th4_ 2 0 0 1
echo "torus4 (H2), DoryNS 1 thread"
time Dory/DoryNS.o Datasets/torus4/torus4_locs.csv 0.15 1 1 Datasets/torus4/DoryNS_H2_th1_ 2 0 0 1
make ripser
time ripser/ripser Datasets/Dragon/dragon2000_locs.csv --dim 2 --threshold 1000
