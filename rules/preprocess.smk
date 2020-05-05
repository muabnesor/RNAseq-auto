
def get_fastqs(wildcards):
    fastq_first = f"{fastq_dict[wildcards.sample]['fastq_first']}"
    fastq_second = f"{fastq_dict[wildcards.sample]['fastq_second']}"
    return [fastq_first, fastq_second]

rule trim_all:
    input:
        expand(f"{trim_galore_dir}/{{sample}}", sample=sample_names)

rule trim:
    input:
        get_fastqs

    output:
        fastq_first = f"{trim_galore_dir}/{{sample}}/{{sample}}_{{lane}}{config['fastq1_suffix_trimmed']}",
        fastq_second = f"{trim_galore_dir}/{{sample}}/{{sample}}_{{lane}}{config['fastq2_suffix_trimmed']}",

    singularity:
        f"{container_dir}/{config['containers']['preprocess_image']}"
    params:
        slurm_log_dir = f"{str(slurm_logdir_preprocess)}",
        trim_dir = f"{trim_galore_dir}/{{sample}}/"

    shell:
        "mkdir -p {params.trim_dir} && trim_galore --illumina --paired --fastqc -o {params.trim_dir} "
        "{input}"

rule trimmed_multiqc:
    input:
        expand(f"{trim_galore_dir}/{{sample}}/", sample=sample_names)

    output:
        directory(f"{trim_galore_dir}/multiqc")

    singularity:
        f"{container_dir}/{config['containers']['preprocess_image']}"

    params:
        slurm_log_dir = f"{str(slurm_logdir_preprocess)}"

    shell:
        "mkdir -p {output} && "
        "multiqc {input} -o {output}/"
