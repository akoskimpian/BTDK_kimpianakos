#!/bin/bash

# Define the path to the chunks directory
CHUNKS_DIR="/mnt/c/TDK/consensus_chunks"
OUTPUT_DIR="/mnt/c/TDK/fseq2_results"
mkdir -p $OUTPUT_DIR

# Loop through each chunk and run fseq2
for chunk in $CHUNKS_DIR/*; do
    fseq2 callpeak $chunk -chrom_size_file hg38.chrom.sizes -name $(basename $chunk) -standard_narrowpeak -l 600 -f 0 -p_thr 0.01 -tp 4.0
    mv $(basename $chunk)_peaks.narrowPeak $OUTPUT_DIR/
done
