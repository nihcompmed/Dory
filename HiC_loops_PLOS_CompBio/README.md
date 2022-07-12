#Hi-C data sets

##Process Hi-C data sets

1. Download pairs file for Hi-C control from https://data.4dnucleome.org/files-processed/4DNFITH9KJF1/ (~33.97 GB) and for Hi-C auxin from https://data.4dnucleome.org/files-processed/4DNFI2X45Z5L/ (~33.05 GB) to this folder.
2. Download cooler CLI (https://cooler.readthedocs.io/en/latest/quickstart.html). Run `pip3 install cooler` in command terminal.
3. Generate Hi-C contact matrix.
```
CHROMSIZES_FILE='4DNFI823LSII.chrom.sizes'
BINSIZE=1000
OUTPUT_FILE='4DNFITH9KJF1.cool'
PAIRS_FILE='4DNFITH9KJF1.txt.gz' 
zcat < $PAIRS_FILE \
    | tr ' ' '\t' \
    | cooler cload pairs -c1 2 -p1 3 -c2 5 -p2 6 --chunksize=2000000 $CHROMSIZES_FILE:$BINSIZE - $OUTPUT_FILE
```
and
```
CHROMSIZES_FILE='4DNFI823LSII.chrom.sizes'
BINSIZE=1000
OUTPUT_FILE='4DNFI2X45Z5L.cool'
PAIRS_FILE='4DNFI2X45Z5L.txt.gz' 
zcat < $PAIRS_FILE \
    | tr ' ' '\t' \
    | cooler cload pairs -c1 2 -p1 3 -c2 5 -p2 6 --chunksize=2000000 $CHROMSIZES_FILE:$BINSIZE - $OUTPUT_FILE
```
4. Balance the HiC-contact matrices.
`cooler balance -p 6 -c 1000000 4DNFITH9KJF1.cool` and `cooler balance -p 6 -c 1000000 4DNFI2X45Z5L.cool`.
5. Run `python3 get_HiC_edges.py` to generate the filtration in sparse format: vertex, vertex, edge-length. This extraction is optimized using Python packages hdf5 and numba.

##Computing PH and H1 loops

1. Run `python3 compute_PH.py` to compute PH, birth-cycles, shorten them, and smooth them.
2. Run `python3 estimate_persistence.py` to determine possibly significant H1 loops.
3. Plotting figures are 'plot_' files.

##Cycle and gene analysis

This is done for birth-cycles computed up to threshold of 150.

1. Run `python3 compute_PH_thresholds.py` to compute PH and birth-cycles.
2. Run `python3 plot_cyc_on_matrix.py` to plot loops (cis-chromosome1) on contact matrix.
3. Run `python3 plot_chromosome_arc.py` to plot loops with long-range interaction on chromosome arc.

