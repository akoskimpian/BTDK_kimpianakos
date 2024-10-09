#!/bin/bash

# Input files
BED_FILE="ENCODE_peaks_in_promoters_with_names.bed"  # Your input BED file with promoters and gene names
BLACKLIST_FILE="hg38-blacklist.v2.bed"  # The blacklist BED file to filter out unwanted regions
FASTA_FILE="TFLink_Homo_sapiens_bindingSites_LS_v1.0.fasta"  # Your input FASTA file with TF binding sites
OUTPUT_FILE="tf_unique_gene_counts_blacklist.txt"  # Output file for unique gene counts per TF
PROMOTER_BINDING_OUTPUT="promoters_with_tf_binding_blacklist.txt"  # Output file for promoters with at least one TF binding
TEMP_DIR="./temp"  # Temporary directory for intermediate files

# Create temporary directory
mkdir -p "$TEMP_DIR"

# Step 1: Subtract regions overlapping with blacklist from the BED file
echo "Filtering blacklist regions from the BED file..."
bedtools subtract -a "$BED_FILE" -b "$BLACKLIST_FILE" > "$TEMP_DIR/ENCODE_peaks_in_promoters_with_names_blacklisted.bed"

# Step 2: Extract valid TF positions from FASTA file and print for debugging
echo "Extracting genomic positions from FASTA file..."
grep '^>' "$FASTA_FILE" | grep 'chr' | \
awk -F'[;:-]' '{if($5 != "") print $5"\t"$6"\t"$7"\t"$3}' | \
awk 'BEGIN{OFS="\t"} {if ($2 ~ /^[0-9]+$/ && $3 ~ /^[0-9]+$/) print $1, $2, $3, $4}' > "$TEMP_DIR/extracted_positions.bed"

# Step 3: Run bedtools intersect to find overlapping regions and capture gene names
echo "Running bedtools intersect to find overlapping regions..."
bedtools intersect -a "$TEMP_DIR/extracted_positions.bed" -b "$TEMP_DIR/ENCODE_peaks_in_promoters_with_names_blacklisted.bed" -f 1.0 -wa -wb | \
awk 'BEGIN{OFS="\t"} {if ($8 != "") print $4, $8}' > "$TEMP_DIR/intersected_tf_gene.txt"

# Step 4: Extract unique gene names for each TF
echo "Counting unique gene names per TF..."
awk '{
    tf = $1
    gene = $2
    if (tf != "" && gene != "") {
        key = tf "-" gene
        if (!seen[key]++) {
            print tf "\t" gene
        }
    }
}' "$TEMP_DIR/intersected_tf_gene.txt" | sort -u > "$TEMP_DIR/tf_gene_pairs.txt"

# Step 5: Count the number of unique genes each TF binds to
awk '{
    count[$1]++
} END {
    for (tf in count) {
        print tf "\t" count[tf]
    }
}' "$TEMP_DIR/tf_gene_pairs.txt" > "$OUTPUT_FILE"

# ---- New functionality to count genes with at least one TF binding ----

# Step 6: Extract unique gene names with at least one TF binding (skip empty lines)
echo "Counting promoters with at least one TF binding..."
awk 'NF > 0 {print $2}' "$TEMP_DIR/tf_gene_pairs.txt" | sort -u > "$TEMP_DIR/unique_genes_with_tf_binding.txt"

# Step 7: Count the number of unique genes with at least one TF binding
UNIQUE_GENE_COUNT=$(wc -l < "$TEMP_DIR/unique_genes_with_tf_binding.txt")

# Step 8: Save the unique gene count to a separate file
echo "Number of genes with at least one TF binding: $UNIQUE_GENE_COUNT"
echo -e "Total genes with at least one TF binding:\t$UNIQUE_GENE_COUNT" > "$PROMOTER_BINDING_OUTPUT"

# Step 9: Final output
if [ -s "$OUTPUT_FILE" ]; then
    echo "TF unique gene counts saved to $OUTPUT_FILE"
else
    echo "No overlaps found. The output file is empty."
fi

if [ -s "$PROMOTER_BINDING_OUTPUT" ]; then
    echo "Promoters with TF binding saved to $PROMOTER_BINDING_OUTPUT"
else
    echo "No promoters with TF binding found. The output file is empty."
fi

# Clean up temporary files
rm -r "$TEMP_DIR"
