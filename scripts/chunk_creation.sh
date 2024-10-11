#!/bin/bash

# Define directories
BIGWIG_FILE="/mnt/c/TDK/consensus_median_filtered.bw"
CHUNK_DIR="/mnt/c/TDK/consensus_chunks_test"
mkdir -p $CHUNK_DIR

# Convert BigWig to BED format (the script assumes you already have 'bigWigToBedGraph' installed)
TEMP_BED="$CHUNK_DIR/temp.bed"
bigWigToBedGraph $BIGWIG_FILE stdout | awk 'BEGIN{OFS="\t"} {print $1, $2, $3, "NAME", $4, "+"}' > $TEMP_BED

# Function to generate two-letter chunk names (aa, ab, ac, ..., zz)
generate_chunk_name() {
    local index=$1
    local letters=({a..z})
    local first_letter_index=$(( index / 26 ))
    local second_letter_index=$(( index % 26 ))
    echo "${letters[first_letter_index]}${letters[second_letter_index]}"
}

# Split BED file into chunks with 10 million lines per file
split -l 1000000 -d $TEMP_BED $CHUNK_DIR/chunk_

# Rename chunks with two-letter suffixes (aa, ab, ..., zz) and .bed extension
counter=0
for chunk in $CHUNK_DIR/chunk_*; do
    chunk_name=$(generate_chunk_name $counter)
    mv "$chunk" "${CHUNK_DIR}/chunk_${chunk_name}.bed"
    counter=$((counter + 1))
done

# Clean up temporary files
rm $TEMP_BED

echo "Chunks created in $CHUNK_DIR with names chunk_aa.bed, chunk_ab.bed, etc."
