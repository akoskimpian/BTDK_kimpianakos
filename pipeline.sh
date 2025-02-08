#!/bin/bash

# Display usage information
usage() {
    echo "Usage: $0 --input <file|dir|list.txt> --output <output_dir> --fasta <binding_sites.fasta>"
    echo
    echo "Required arguments:"
    echo "  --input       Input BED file(s): a file, a directory, or a text file with paths."
    echo "  --output      Directory to store output files."
    echo "  --fasta       FASTA file containing TF binding sites."
    echo
    echo "Optional arguments:"
    echo "  -h, --help    Show this help message and exit."
    exit 1
}

# Parse command-line arguments
while [[ "$#" -gt 0 ]]; do
    case $1 in
        --input) INPUT="$2"; shift ;;
        --output) OUTPUT_DIR="$2"; shift ;;
        --fasta) FASTA_FILE="$2"; shift ;;
        -h|--help) usage ;;
        *) echo "Unknown parameter passed: $1"; usage ;;
    esac
    shift
done

# Validate inputs
if [[ -z "$INPUT" || -z "$OUTPUT_DIR" || -z "$FASTA_FILE" ]]; then
    echo "Error: Missing required arguments."
    usage
fi

# Create output directories
mkdir -p "$OUTPUT_DIR"
RESULTS_DIR="$OUTPUT_DIR"
SORTED_DIR="$OUTPUT_DIR/sorted"
mkdir -p "$RESULTS_DIR" "$SORTED_DIR"

# Define URLs for required files
PROMOTERS_URL="https://raw.githubusercontent.com/akoskimpian/BTDK_kimpianakos/main/promoters.bed"
BLACKLIST_URL="https://raw.githubusercontent.com/akoskimpian/BTDK_kimpianakos/main/blacklist.bed"

# Define local file paths
PROMOTERS_FILE="$RESULTS_DIR/promoters.bed"
BLACKLIST_FILE="$RESULTS_DIR/blacklist.bed"

# Function to download if the file does not exist
download_if_missing() {
    local file=$1
    local url=$2
    if [[ ! -f "$file" ]]; then
        curl -sS -o "$file" "$url" || { echo "Failed to download $file"; exit 1; }
    fi
}

# Download required files
download_if_missing "$PROMOTERS_FILE" "$PROMOTERS_URL"
download_if_missing "$BLACKLIST_FILE" "$BLACKLIST_URL"

# Standardize chromosome names in promoters file
awk '{ if($1 !~ /^chr/) $1="chr"$1; print }' "$PROMOTERS_FILE" > "$PROMOTERS_FILE.tmp" && mv "$PROMOTERS_FILE.tmp" "$PROMOTERS_FILE"

# Step 1: Gather input files
INPUT_FILES=()
if [[ -f "$INPUT" && "$INPUT" == *.txt ]]; then
    while IFS= read -r line; do INPUT_FILES+=("$line"); done < "$INPUT"
elif [[ -d "$INPUT" ]]; then
    INPUT_FILES=("$INPUT"/*.bed)
elif [[ -f "$INPUT" ]]; then
    INPUT_FILES+=("$INPUT")
else
    echo "Error: Invalid input. Provide a file, directory, or file list."
    exit 1
fi

# Step 2: Sort input BED files
SORTED_PATHS_FILE="$RESULTS_DIR/sorted_paths.txt"
rm -f "$SORTED_PATHS_FILE"
for bedfile in "${INPUT_FILES[@]}"; do
    if [[ -f "$bedfile" ]]; then
        filename=$(basename "$bedfile" .bed)
        sorted_bedfile="$SORTED_DIR/${filename}_sorted.bed"
        awk 'BEGIN { OFS="\t" } { print $1, $2, $3 }' "$bedfile" | sort -k1,1 -k2,2n > "$sorted_bedfile"
        echo "$sorted_bedfile" >> "$SORTED_PATHS_FILE"
    fi
done

# Step 3: Create consensus open chromatin regions
CONSENSUS_FILE="$RESULTS_DIR/consensus.bed"
bedtools multiinter -i $(cat "$SORTED_PATHS_FILE") > "$RESULTS_DIR/temp_intersections.bed"
awk -v total_files=$(wc -l < "$SORTED_PATHS_FILE") 'BEGIN { threshold = int((total_files + 1) / 2) } {
    if ($4 >= threshold && $3 - $2 >= 5) {
        print $1"\t"$2"\t"$3
    }
}' "$RESULTS_DIR/temp_intersections.bed" | bedtools merge -i - > "$CONSENSUS_FILE"

# Standardize chromosome names in consensus file
awk '{ if($1 !~ /^chr/) $1="chr"$1; print }' "$CONSENSUS_FILE" > "$CONSENSUS_FILE.tmp" && mv "$CONSENSUS_FILE.tmp" "$CONSENSUS_FILE"

# Step 4: Intersect consensus with promoters
PROMOTERS_SORTED="$RESULTS_DIR/promoters.sorted.bed"
sort -k1,1 -k2,2n "$PROMOTERS_FILE" > "$PROMOTERS_SORTED"
OPEN_PROMOTERS="$RESULTS_DIR/open_promoters.bed"
bedtools intersect -a "$CONSENSUS_FILE" -b "$PROMOTERS_SORTED" -wa -wb | awk 'NF>=7 { print $1 "\t" $2 "\t" $3 "\t" $7 }' > "$OPEN_PROMOTERS"

# Step 5: Extract TF binding site positions from FASTA
TF_POSITIONS="$RESULTS_DIR/tf_positions.bed"
grep '^>' "$FASTA_FILE" | grep 'chr' | awk -F'[;:-]' '{if($5 != "") print $5"\t"$6"\t"$7"\t"$3}' | awk 'BEGIN{OFS="\t"} {if ($2 ~ /^[0-9]+$/ && $3 ~ /^[0-9]+$/) print $1, $2, $3, $4}' > "$TF_POSITIONS"

# Step 6: Find TF-gene overlaps
TF_GENE_OUTPUT="$RESULTS_DIR/tf_gene_pairs.txt"
bedtools intersect -a "$TF_POSITIONS" -b "$OPEN_PROMOTERS" -wa -wb | awk 'BEGIN{OFS="\t"} {print $4, $8}' | sort -u > "$TF_GENE_OUTPUT"

# Step 7: Count unique genes per TF
TF_UNIQUE_GENE_COUNTS="$RESULTS_DIR/tf_unique_gene_counts.txt"
awk '{count[$1]++} END {for (tf in count) print tf "\t" count[tf]}' "$TF_GENE_OUTPUT" > "$TF_UNIQUE_GENE_COUNTS"

# Cleanup
rm -r "$RESULTS_DIR/temp_intersections.bed"
echo "All tasks completed successfully."
