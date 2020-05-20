
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
        trim_dir = directory(f"{trim_galore_dir}/{{sample}}/"),

    singularity:
        f"{container_dir}/{config['containers']['preprocess_image']}"
    params:
        slurm_log_dir = f"{str(slurm_logdir_preprocess)}",

    shell:
        "mkdir -p {output.trim_dir} && trim_galore --stringency 5 --illumina --paired --fastqc -o {output.trim_dir} "
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
