#!/bin/bash

# Input files
PROMOTERS_FILE="promoters_only_important_columns.bed"
OPEN_CHROMATIN_FILE="colon_intersections_final.bed"

# Output file
OUTPUT_FILE="colon_open_promoters_test.bed"

# Use bedtools intersect to find overlaps and keep only the overlapping regions within the promoters
bedtools intersect -a $OPEN_CHROMATIN_FILE -b $PROMOTERS_FILE -wa -wb | \
awk '{gsub(/^,/, "", $7); print $1"\t"$2"\t"$3"\t"$7}' > $OUTPUT_FILE

echo "Filtered open regions within promoters saved to $OUTPUT_FILE"

