#!/bin/bash
# Script: calculate_methylated_c_positions.sh
# Author: Ruining Dong
# Date: 2024-06-26
# Summary: This script parse MM tags from Biomodal bam file and label CpG index per read.

source /g/data/pq08/software/mambaforge/etc/profile.d/conda.sh
conda activate /g/data/pq08/software/mambaforge/envs/methyl_env
module load samtools

# input_bam=/g/data/pq08/projects/biomodal/data_bucket/test-runxyz/nf-results/duet-1.2.1_rd6078_2024-06-06_1408_5bp/sample_outputs/bams/CEG9330132-19-01.genome.GRCh38Decoy_primary_assembly.dedup.bam
# output_bed=/g/data/pq08/projects/biomodal/data_bucket/test-runxyz.MM.bed
# output_anno_bed=/g/data/pq08/projects/biomodal/data_bucket/test-runxyz.MM.annotated.bed
# threads=4

# Function to display usage information
usage() {
    echo "Usage: $0 [--threads <threads>] [-o <output_file>] <input_bam>"
    exit 1
}

THREADS=1
OUTPUT_FILE=""

# Parse the command-line arguments
while [[ "$#" -gt 0 ]]; do
    case $1 in
        --threads)
            THREADS=$2
            shift
            ;;
        -o)
            OUTPUT_FILE=$2
            shift
            ;;
        *)
            if [[ -z "$INPUT_BAM" ]]; then
                INPUT_BAM=$1
            else
                usage
            fi
            ;;
    esac
    shift
done

# Check if the input BED file is provided
if [[ -z "$INPUT_BAM" ]]; then
    usage
fi

# mamba activate methyl_env
parse_mm_tags() {
    local input_bam=$1
    local threads=$2 
    samtools view -@ $threads $input_bam | awk '{
        OFS="\t";
        mm_tag = ($15 ~ /^MM:Z:/ ? $15 : ($16 ~ /^MM:Z:/ ? $16 : "MM tag not found"));
        print $3, $4, $4 + length($10) - 1, mm_tag
    }' > tmp_mm_tags &&
    wgbstools convert -L tmp_mm_tags --parsable --drop_empty -@ $threads &&
    rm tmp_mm_tags
}
export -f parse_mm_tags

if [[ -z "$OUTPUT_FILE" ]]; then
    parse_mm_tags $INPUT_BAM $THREAD
else
    parse_mm_tags $INPUT_BAM $THREADS > $OUTPUT_FILE
fi

