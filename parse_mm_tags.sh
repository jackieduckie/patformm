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
# tmp_mm_tags="tmp_mm_tags_1_10kbin"
mm_dir="/g/data/pq08/projects/biomodal/patformm"

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
    # echo "Running: samtools view -@ $threads $input_bam"
    samtools view -@ $threads $input_bam | awk '{
        OFS="\t";
        mm_tag = ($15 ~ /^MM:Z:/ ? $15 : ($16 ~ /^MM:Z:/ ? $16 : "MM_tag_not_found"));
        print $3, $4, $4 + length($10) - 1, mm_tag
    }'  
}
export -f parse_mm_tags

# echo "WARNING: temporarily skipping generating tmp_mm_tags for testing purposes"

tmp_mm_tags="$(basename $INPUT_BAM).tmp_mm_tags"
echo "Running: parse_mm_tags $INPUT_BAM $THREADS > ${mm_dir}/${tmp_mm_tags}"
parse_mm_tags $INPUT_BAM $THREADS > ${mm_dir}/${tmp_mm_tags}

# Count the number of rows in the input BED file
num_rows=$(wc -l < "${mm_dir}/$tmp_mm_tags")

export TMPDIR=/scratch/pq08/rd6078/tmp
tmp_dir=$(mktemp -d )

# Check if the number of rows exceeds 2 million
if [[ "$num_rows" -le 2000000 ]]; then
    echo "The input BED file has $num_rows rows, which is less than or equal to 2 million rows. No need to split."
    wgbstools convert -L "${mm_dir}/$tmp_mm_tags" --parsable --drop_empty --out_path "${tmp_dir}/${tmp_mm_tags}.concat.bed"
else
    # Split tmp_mm_tags in to 1M row small beds
    echo "Spliting mm_tags bed file to temp dir $tmp_dir."
    split -l 1000000 "$tmp_mm_tags" "$tmp_dir/tmp_mm_tags_chunk_"
    # Process each chunk in parallel
    echo "Annotate CpG id with wgbstools"
    find "$tmp_dir" -name 'tmp_mm_tags_chunk_*' | xargs -n 1 -P "$THREADS" -I {} wgbstools convert -L "{}" --parsable --drop_empty --out_path "{}.out" &&
    cat "$tmp_dir"/tmp_mm_tags_chunk_*.out > "${tmp_dir}/${tmp_mm_tags}.concat.bed"
fi


if [[ -z "$OUTPUT_FILE" ]]; then
    echo "tmp files are located at: ${tmp_mm_tags} and ${tmp_dir}:" &&
	    echo ${tmp_dir}/${tmp_mm_tags}.concat.bed
else
    cp ${tmp_dir}/${tmp_mm_tags}.concat.bed $OUTPUT_FILE &&
	    rm ${tmp_mm_tags} &&
	    rm -r ${tmp_dir} &&
	    echo "tmp files ${tmp_mm_tags} and ${tmp_dir} deleted"
fi

