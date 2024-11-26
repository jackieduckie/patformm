#PBS -l ncpus=8
#PBS -l mem=80GB
#PBS -l jobfs=20GB
#PBS -q normal
#PBS -l walltime=24:00:00
#PBS -l storage=scratch/pq08+gdata/pq08
#PBS -P pq08
#PBS -m bea
#PBS -M ruining.dong@unimelb.edu.au

patformm_path=/g/data/pq08/projects/biomodal/patformm
patformm=$patformm_path/patformm
# input_path=/g/data/pq08/projects/biomodal/data_bucket/wgs_cups_1/nf-results/duet-1.2.1_rd6078_2024-05-27_2131_5bp/sample_outputs/bams
input_path=$patformm_path
bam=BM-BHE008-mini.bam
cd $patformm_path
threads=1

# while IFS=$'\t' read -r input_path bam; do
    input_bam=${input_path}/${bam}
    output_bed=$patformm_path/${bam}.bed
    output_pat=$patformm_path/${bam}.fast.pat
    # echo "running patformm parse_mm_tags on $input_bam" &&
    # $patformm parse_mm_tags --threads $threads -o $output_bed $input_bam &&    
    echo "running patformm calculate_cpos on $bam" &&
    $patformm calculate_cpos --threads $threads -o ${output_pat} ${output_bed} &&
    echo "all done."
# done < $patformm_path/bam_list
