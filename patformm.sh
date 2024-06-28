#!/bin/bash

input_bam=/g/data/pq08/projects/biomodal/data_bucket/test-runxyz/nf-results/duet-1.2.1_rd6078_2024-06-06_1408_5bp/sample_outputs/bams/CEG9330132-19-01.genome.GRCh38Decoy_primary_assembly.dedup.bam
output_bed=/g/data/pq08/projects/biomodal/data_bucket/test-runxyz.MM.bed
output_anno_bed=/g/data/pq08/projects/biomodal/data_bucket/test-runxyz.MM.annotated.bed
threads=4

# mamba activate methyl_env
samtools view -@ $threads $input_bam | awk '{
    OFS="\t";
    mm_tag = ($15 ~ /^MM:Z:/ ? $15 : ($16 ~ /^MM:Z:/ ? $16 : "MM tag not found"));
    print $3, $4, $4 + length($10) - 1, mm_tag
}' > $output_bed 

wgbstools convert -L $output_bed --parsable --drop_empty -@ $threads > $output_anno_bed


