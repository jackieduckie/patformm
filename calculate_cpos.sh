#!/bin/bash
# Script: calculate_methylated_c_positions.sh
# Author: Ruining Dong
# Date: 2024-06-26
# Summary: This script processes an input BED file to calculate the positions of methylated 'C's 
# based on the MM tag, converts them to CpG index values using wgbstools, and handles invalid genomic regions.
# input bed is generated by another script: patformm


# mamba activate methyl_env
source /g/data/pq08/software/mambaforge/etc/profile.d/conda.sh
conda activate /g/data/pq08/software/mambaforge/envs/methyl_env
module load samtools
patformm_path=/g/data/pq08/projects/biomodal/patformm

# Function to display usage information
usage() {
    echo "Usage: $0 [--threads <threads>] [-o <output_file>] <input_bed_file>"
    exit 1
}

# Default number of threads
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
            if [[ -z "$INPUT_BED" ]]; then
                INPUT_BED=$1
            else
                usage
            fi
            ;;
    esac
    shift
done

# Check if the input BED file is provided
if [[ -z "$INPUT_BED" ]]; then
    usage
fi


# Create a temporary directory for split files
export TMPDIR=/scratch/pq08/rd6078/tmp
tmp_dir=$(mktemp -d)
# tmp_dir=/scratch/pq08/rd6078/tmp/tmp.smPmXC0Ndp
# echo "WARNING: resuming previously failed run for testing purposes: $tmp_dir"

# Split the input BED file into smaller chunks
# echo "WARNING: temporarily commenting below commands to resume failed run"
echo "split BED to $tmp_dir"
split -l $(( $(wc -l < "$INPUT_BED") / THREADS + 1 )) "$INPUT_BED" "$tmp_dir/bed_chunk_"

# Process each chunk in parallel and Combine the results
echo "process temp bed chunks"
find "$tmp_dir" -maxdepth 1 -name 'bed_chunk_*' | xargs -n 1 -P "$THREADS" -I {} $patformm_path/process_chunk.sh {} &&
	cat "$tmp_dir"/*.out > $tmp_dir/tmp_final_output.bed


# find "$tmp_dir" -name 'bed_chunk_*' | xargs -n 1 -P "$THREADS" bash -c 'process_chunk "$0"' _
# find "$tmp_dir" -name 'bed_chunk_*' | xargs -n 1 -P "$THREADS" -I {} bash -c 'process_chunk "$@"' _ {}



if [[ -z "$OUTPUT_FILE" ]]; then
    awk -F'\t' '$8 != "NA" {print $1"\t"$5"\t"$8}' $tmp_dir/tmp_final_output.bed | sort -k2,2n -k3,3 | uniq -c | awk '{print $2"\t"$3"\t"$4"\t"$1}'
else
    awk -F'\t' '$8 != "NA" {print $1"\t"$5"\t"$8}' $tmp_dir/tmp_final_output.bed | sort -k2,2n -k3,3 | uniq -c | awk '{print $2"\t"$3"\t"$4"\t"$1}' > $OUTPUT_FILE && 
	    gzip $OUTPUT_FILE
fi

# rm -r "$tmp_dir"
# rm tmp_output.bed
