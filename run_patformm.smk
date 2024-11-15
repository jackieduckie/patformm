import os

# Define the path to the patformm executable and the base path
patformm_path = "/g/data/pq08/projects/biomodal/patformm"
patformm = os.path.join(patformm_path, "patformm")

configfile: "test_config.yaml"

# Define the rule to run all processes
rule all:
    input:
        expand("output/{sample}.done", sample=config["bam_files"].keys())
        # expand("{bam}.pat.gz", bam=[bam for _, bam in bam_list])

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
        {params.patformm} calculate_cpos --threads {threads} -o {params.pat} {input.bed} &&
        touch {output.done} &&
        echo "{wildcards.sample} done."
        """

