#!/bin/bash

# Add UCSC tools to PATH
export PATH=/mnt/c/TDK/UCSC_tools:$PATH

# Define directories
BIGWIG_DIR="/mnt/c/TDK"  
BEDGRAPH_DIR="$BIGWIG_DIR/bedgraph"
BED_DIR="$BIGWIG_DIR/bed"

# Create directories if they do not exist
mkdir -p "$BEDGRAPH_DIR"
mkdir -p "$BED_DIR"

# BigWig file path
bigwig_file="/mnt/c/TDK/consensus_median_filtered.bw"

# Convert the BigWig file to BedGraph and then to BED
echo "Processing file: $bigwig_file"
if [ ! -f "$bigwig_file" ]; then
    echo "File does not exist: $bigwig_file"
    exit 1
fi

base=$(basename "$bigwig_file" .bw)
bedgraph="$BEDGRAPH_DIR/$base.bedgraph"

echo "Converting $bigwig_file to BedGraph..."
bigWigToBedGraph "$bigwig_file" "$bedgraph"

if [ $? -ne 0 ]; then
    echo "Failed to convert $bigwig_file to BedGraph"
    exit 1
else
    echo "Successfully converted $bigwig_file to $bedgraph"
fi

echo "Converting $bedgraph to BED..."
# Include the fourth column (score) in the BED file
awk '{OFS="\t"; print $1, $2, $3, $4}' "$bedgraph" > "$BED_DIR/$base.bed"

if [ $? -ne 0 ]; then
    echo "Failed to convert $bedgraph to BED"
    exit 1
else
    echo "Successfully converted $bedgraph to $BED_DIR/$base.bed"
fi

# Remove the BedGraph file after conversion
rm "$bedgraph"
echo "Removed $bedgraph"

echo "Conversion completed. BED file is in $BED_DIR/$base.bed."
