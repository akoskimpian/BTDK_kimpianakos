#!/bin/bash

# Define the directory with the fseq2 results and output file
RESULTS_DIR="/mnt/c/TDK/fseq2_results"
MERGED_FILE="/mnt/c/TDK/fseq2_merged_peaks.narrowPeak"

# Merge all narrowPeak files in the results directory into one file
cat $RESULTS_DIR/*.narrowPeak | bedtools sort -i - | bedtools merge -i - -c 4,5,6 -o distinct > $MERGED_FILE

echo "Merged narrowPeak files into: $MERGED_FILE"
