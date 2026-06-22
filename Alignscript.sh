#!/bin/bash
set -euo pipefail

# Usage check
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <reference.fasta> <reads.fasta>"
    exit 1
fi

# Input files
REFERENCE=$1
READS=$2

# Output directory
OUTDIR="results"
mkdir -p "$OUTDIR"

# Base names for output files
BASE_NAME=$(basename "$READS" .fasta)
BWA_SAI="$OUTDIR/$BASE_NAME.bwa.sai"
BWA_SAM="$OUTDIR/$BASE_NAME.bwa.sam"
SORTED_BAM="$OUTDIR/$BASE_NAME.bwa.sort.bam"
DEPTH_TXT="$OUTDIR/$BASE_NAME.depth.txt"
BED_FILE="$OUTDIR/$BASE_NAME.bed"
PMD_TXT="$OUTDIR/$BASE_NAME.pmd.txt"
ALLELE_FREQ_OUT="$OUTDIR/$BASE_NAME"

# Index the reference
bwa index "$REFERENCE"

# Align reads to the reference
bwa aln -n 0.01 -o 2 -l 1024 -t 12 "$REFERENCE" "$READS" > "$BWA_SAI"

# Generate SAM file
bwa samse "$REFERENCE" "$BWA_SAI" "$READS" > "$BWA_SAM"

# Convert SAM to sorted BAM and index
samtools view -h "$BWA_SAM" | samtools sort -o "$SORTED_BAM"
samtools index "$SORTED_BAM"

# Generate depth information
samtools depth -a "$SORTED_BAM" > "$DEPTH_TXT"

# Convert BAM to BED
bedtools bamtobed -i "$SORTED_BAM" > "$BED_FILE"

# Run metaDMG-cpp PMD analysis
# NOTE: verify that stdout redirection is the correct output method for your
# version of metaDMG-cpp — some versions expect an -out flag instead
metaDMG-cpp pmd "$SORTED_BAM" > "$PMD_TXT"

# Index the reference for ANGSD
samtools faidx "$REFERENCE"

# Run ANGSD for allele frequency estimation
angsd -out "$ALLELE_FREQ_OUT" -i "$SORTED_BAM" -GL 1 -doMaf 8 -doMajorMinor 4 -ref "$REFERENCE" -doCounts 1

echo "Workflow completed successfully. Output files are in $OUTDIR/"
