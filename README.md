# TFLink 2.0 ATAC-seq Pipeline



This repository contains scripts developed for processing ATAC-seq data to identify transcription factor - target gene interactions in a tissue-specific manner. The pipeline uses ATAC-seq data to predict the interactions likely to occur in a given organ/tissue/cell type.

<br />

**Prerequisites**

The following tools are required to run the scripts in this pipeline:

-wiggletools

-wigToBigWig

-bigWigToBedGraph

-fseq2

-bedtools

<br />

**Input Files**

ATAC-seq BED files: Represent open chromatin regions identified in tissue-specific samples

Promoter regions BED file: Define promoter locations (e.g. from Ensembl; GRCh38.p14)

Blacklist regions: Known problematic regions to exclude from analysis (e.g. from ENCODE)

FASTA containing TF binding sites: Contains transcription factor binding site sequences and coordinates (e.g. downloadable from TFLink, originating from ChIP-seq data)

<br />

**Output Files**

Consensus BED file: Open chromatin regions present in at least 50% of the input samples

TF-gene interaction pairs: File containing TF names and their corresponding target genes

Unique genes per TF: File listing each TF and the number of genes it interacts with
<br />


**Example usage**

```
./pipeline.sh --input <input_file|directory|list.txt> --output <output_dir>  
```
<br />

**Notes**

Developed by Ákos Kimpián as part of Eszter Ari's research group


