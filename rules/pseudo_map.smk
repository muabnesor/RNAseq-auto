def get_fastqs(wildcards):
    fastq_first = f"{trim_galore_dir}/{wildcards.sample}/{fastq_dict[wildcards.sample]['fastq_first'].split('/')[-1].replace(config['fastq1_suffix'], config['fastq1_suffix_trimmed'])}"
    fastq_second = f"{trim_galore_dir}/{wildcards.sample}/{fastq_dict[wildcards.sample]['fastq_second'].split('/')[-1].replace(config['fastq2_suffix'], config['fastq2_suffix_trimmed'])}"
    return [fastq_first, fastq_second]


rule pseudo_map:
    input:
        fastq_files = f"{trim_galore_dir}/{{sample}}/",
        index_dir = expand_path(reference_dir, "salmon_index"),
        transcript_file = expand_path(reference_dir, config["references"]["transcripts_file"])

    output:
        quant_dir = directory(f"{pseudo_map_dir}/{{sample}}")

    singularity:
        f"{container_dir}/{config['containers']['salmon_image']}"

    threads: config["quant"]["salmon_quant_threads"]

    params:
        slurm_log_dir = f"{str(slurm_logdir_pseudo_map)}",
        fastq_files = get_fastqs

    shell:
        "mkdir -p {output.quant_dir} && "
        "salmon quant -p {threads} "
        "-l A "
        "-i {input.index_dir} "
        "--seqBias --gcBias "
        "--validateMappings "
        "-o {output.quant_dir} "
        "-1 {params.fastq_files[0]} "
        "-2 {params.fastq_files[1]} "


rule merge_quant_tpm:
    input:
        expand(f"{pseudo_map_dir}/{{sample}}", sample=sample_names)

    output:
        f"{pseudo_map_dir}/merged_tpm.sf"

    singularity:
        f"{container_dir}/{config['containers']['salmon_image']}"

    params:
        slurm_log_dir = f"{str(slurm_logdir_pseudo_map)}",

    threads: 1

    shell:
        "salmon quantmerge "
        "--quants {input} "
        "-c tpm "
        "-o {output}"


rule merge_quant_numreads:
    input:
        expand(f"{pseudo_map_dir}/{{sample}}", sample=sample_names)

    output:
        f"{pseudo_map_dir}/merged_numreads.sf"

    singularity:
        f"{container_dir}/{config['containers']['salmon_image']}"

    params:
        slurm_log_dir = f"{str(slurm_logdir_pseudo_map)}",

    threads: 1

    shell:
        "salmon quantmerge "
        "--quants {input} "
        "-c numreads "
        "-o {output}"


rule salmon_deseq2:
    input:
        gene_db = expand_path(reference_dir, "genes.sqllite"),
        count_files = expand(f"{pseudo_map_dir}/{{sample}}", sample=sample_names),
        sample_data = f"{sample_data}",
        transcripts_gtf = expand_path(reference_dir, config["references"]["transcripts_file"])

    output:
        f"{results_dir}/salmon-deseq2/salmon-deseq2.html"

    singularity: f"{container_dir}/{config['containers']['R_image']}"

    params:
        slurm_log_dir = f"{str(slurm_logdir_count)}",
        count_dir = f"{pseudo_map_dir}"

    threads: 1

    script: "../scripts/salmon-deseq2.Rmd"


rule salmon_drimseq:
    input:
        gene_db = expand_path(reference_dir, "genes.sqllite"),
        count_files = expand(f"{pseudo_map_dir}/{{sample}}", sample=sample_names),
        sample_data = f"{sample_data}",
        transcripts_gtf = expand_path(reference_dir, config["references"]["transcripts_file"])

    output:
        f"{results_dir}/salmon-drimseq/salmon-drimseq.html"

    singularity: f"{container_dir}/{config['containers']['R_image']}"

    params:
        slurm_log_dir = f"{str(slurm_logdir_count)}",
        count_dir = f"{pseudo_map_dir}"

    threads: 1

    script: "../scripts/salmon-drimseq.Rmd"
