#!/bin/sh
echo "Compiling Ripser"
make ripser
echo "Dragon, Ripser"
time ripser/ripser Datasets/Dragon/dragon2000_locs.csv --dim 1 --threshold 1000 --format point-cloud
echo "Fractal, Ripser"
time ripser/ripser Datasets/Fract/fractal_r_distmat.csv --dim 2 --threshold 1000 --format distance
echo "o3, Ripser"
time ripser/ripser Datasets/o3/o3_8192.csv --dim 2 --threshold 1 --format point-cloud
echo "torus (H1), Ripser"
time ripser/ripser Datasets/torus4/torus4_locs.csv --dim 1 --threshold 0.15 --format point-cloud
echo "torus (H2), Ripser"
time ripser/ripser Datasets/torus4/torus4_locs.csv --dim 2 --threshold 0.15 --format point-cloud
