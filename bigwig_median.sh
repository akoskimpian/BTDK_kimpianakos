#!/bin/bash

# Define the input BigWig files
input_files=(
    "/mnt/c/TDK/colon_bigwigs/ENCFF163IDS.bigwig"
    "/mnt/c/TDK/colon_bigwigs/ENCFF243ZNH.bigwig"
    "/mnt/c/TDK/colon_bigwigs/ENCFF275VQX.bigwig"
    "/mnt/c/TDK/colon_bigwigs/ENCFF299MFR.bigwig"
    "/mnt/c/TDK/colon_bigwigs/ENCFF355IHN.bigwig"
    "/mnt/c/TDK/colon_bigwigs/ENCFF449DFP.bigwig"
    "/mnt/c/TDK/colon_bigwigs/ENCFF469WCM.bigwig"
    "/mnt/c/TDK/colon_bigwigs/ENCFF575JBQ.bigwig"
    "/mnt/c/TDK/colon_bigwigs/ENCFF628ARL.bigwig"
    "/mnt/c/TDK/colon_bigwigs/ENCFF686EAS.bigwig"
    "/mnt/c/TDK/colon_bigwigs/ENCFF807CHN.bigwig"
    "/mnt/c/TDK/colon_bigwigs/ENCFF969UDD.bigwig"
)

# Define the output file
output_file="consensus_median_test.wig"

# Create the command to compute the median
cmd="wiggletools median"

# Append each input file to the command
for file in "${input_files[@]}"; do
    cmd+=" $file"
done

# Redirect the output to the specified file
eval "$cmd > $output_file"

# Print a message indicating completion
echo "Consensus median computed and saved to $output_file"
