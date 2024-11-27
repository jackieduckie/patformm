import os

# Define the path to the patformm executable and the base path
patformm_path = "/g/data/pq08/projects/biomodal/patformm"
patformm = os.path.join(patformm_path, "patformm")
wgbstools = "/g/data/pq08/software/mambaforge/envs/methyl_env/bin/wgbstools"

configfile: "config_271124.yaml"
# configfile: "config.yaml"
# configfile: "test_config.yaml"

# Define the rule to run all processes
rule all:
    input:
        expand("output/{sample}.done", sample=config["bam_files"].keys()),
        expand("output_homog/{sample}/uxm.bed.gz", sample=config["bam_files"].keys()),
        expand("output/{sample}.beta.gz", sample=config["bam_files"].keys())

# Rule to parse_mm_tags
rule parse_mm_tags:
    input:
        lambda wildcards: "{input_path}/{bam}".format(input_path=config["bam_files"][wildcards.sample]["input_path"], bam=config["bam_files"][wildcards.sample]["bam"])
        # bam=lambda wildcards: os.path.join(wildcards.input_path, wildcards.bam)
    output:
        bed=temp("/scratch/pq08/rd6078/patformm_tmp/{sample}.bed")
        # bed=lambda wildcards: os.path.join(patformm_path, f"{wildcards.bam}.bed")
    params:
        patformm=patformm
    threads: 8
    resources:
        mem_mb=64000,
        walltime=36000,
        ncpus=8,
        jobfs="5G"
    shell:
        """
        echo "running patformm parse_mm_tags on {input}" &&
        {params.patformm} parse_mm_tags --threads {threads} -o {output.bed} {input}
        """

# Rule to calculate_cpos
rule calculate_cpos:
    input:
        bed=rules.parse_mm_tags.output.bed 
        # bed="{bam}.bed"
    output:
        pat_gz="output/{sample}.pat.gz",
        done="output/{sample}.done"
    params:
        patformm=patformm,
        chunk_size=10000000,
        pat="output/{sample}.pat"
    threads: 8
    log: 
        stderr="/g/data/pq08/projects/biomodal/patformm/snakemake_logs/{sample}.stderr",
        stdout="/g/data/pq08/projects/biomodal/patformm/snakemake_logs/{sample}.stdout"
    # cluster: 
    #     queue="hugemem"
    resources:
        mem_mb=100000,
        walltime=144000,
        jobfs="5G"
    shell:
        """
        echo "running patformm calculate_cpos on {input.bed}" &&
        {params.patformm} calculate_cpos --threads {threads} --chunk-size {params.chunk_size} -o {params.pat} {input.bed} &&
        touch {output.done} &&
        echo "{wildcards.sample} done."
        """

# Rule to perform homogeneity analysis
rule homogeneity_analysis:
    input:
        pat_gz=rules.calculate_cpos.output.pat_gz
    output:
        uxm="output_homog/{sample}/uxm.bed.gz"
    params:
        output_dir="output_homog",
        wgbstools=wgbstools,
        bed_file="data/wgbs_probes_methylonco_segments.chunk.bed.gz"
    threads: 8
    resources:
        mem_mb=64000,
        walltime=36000,
        jobfs="5G"
    log:
        stderr="/g/data/pq08/projects/biomodal/patformm/snakemake_logs/{sample}_homogeneity.stderr",
        stdout="/g/data/pq08/projects/biomodal/patformm/snakemake_logs/{sample}_homogeneity.stdout"
    shell:
        """
        echo "Running homogeneity analysis for {wildcards.sample}" &&
        {params.wgbstools} homog \
        {input.pat_gz} \
        -b {params.bed_file} \
        -o {params.output_dir} \
        --thresholds 0.25,0.75
        """

# Rule to convert pat to beta values
rule pat2beta:
    input:
        pat_gz=rules.calculate_cpos.output.pat_gz
    output:
        beta="output/{sample}.beta"
    params:
        wgbstools=wgbstools
    threads: 8
    resources:
        mem_mb=64000,
        walltime=36000,
        jobfs="5G"
    log:
        stderr="/g/data/pq08/projects/biomodal/patformm/snakemake_logs/{sample}_pat2beta.stderr",
        stdout="/g/data/pq08/projects/biomodal/patformm/snakemake_logs/{sample}_pat2beta.stdout"
    shell:
        """
        echo "Converting pat to beta values for {wildcards.sample}" &&
        {params.wgbstools} pat2beta \
        -o output/ \
        {input.pat_gz}
        """

