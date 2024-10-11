#!/bin/bash

# Input and output WIG file names
input_wig="consensus_median.wig"
filtered_wig="filtered_consensus_median.wig"

# Step 1: Filter out random scaffolds and associated fixedStep entries
echo "Filtering random scaffolds and associated fixedStep entries from $input_wig..."
awk '
    # If the line contains a random scaffold and no "step", delete the line
    /random|Un|KI|GL/ && !/step/ { next }

    # If the line contains "step" and a random scaffold, set a flag to skip subsequent lines
    /random|Un|KI|GL/ && /step/ { skip = 1; next }

    # If we are in a skipping state, continue skipping until a line with 4 fields (4 words)
    skip {
        if (NF == 4) {
            skip = 0
        }
        next
    }

    # Print the line if we are not skipping
    { print }
' "$input_wig" > "$filtered_wig"

# Step 2: Convert filtered WIG to BigWig
echo "Converting filtered WIG to BigWig..."
wigToBigWig "$filtered_wig" hg38.chrom.sizes consensus_median_filtered.bw

if [ $? -eq 0 ]; then
    echo "Conversion successful. Output saved to consensus_median_filtered.bw."
else
    echo "Error during BigWig conversion."
fi
