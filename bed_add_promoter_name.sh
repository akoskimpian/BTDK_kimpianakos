#!/bin/bash

# Input files
PEAKS_BED="ENCODE_peaks_in_promoters.bed"  # Original peaks BED file
PROMOTERS_BED="promoters_only_important_columns.bed"  # Promoters BED file with gene names
OUTPUT_BED="ENCODE_peaks_in_promoters_2.bed"  # Output file

# Step 1: Use bedtools intersect to map peaks to promoters and replace the peak names with promoter gene names
bedtools intersect -a "$PEAKS_BED" -b "$PROMOTERS_BED" -wa -wb | \
awk 'BEGIN{OFS="\t"} {print $1, $2, $3, $10, $5, $6}' > "$OUTPUT_BED"

# Step 2: Check if the output file is created and contains data
if [ -s "$OUTPUT_BED" ]; then
    echo "Peaks with promoter names saved to $OUTPUT_BED"
else
    echo "No peaks matched the promoters. The output file is empty."
fi
