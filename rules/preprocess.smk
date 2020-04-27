
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
        f"{trim_galore_dir}/{{sample}}"

    singularity:
        f"{container_dir}/{config['preprocess_image']}"

    params:
        slurm_log_dir = f"{analysis_dir}/logs/trim/slurm/"

    shell:
        "trim_galore --illumina --paired --fastqc -o {output} "
        "{input}"

rule trimmed_multiqc:
    input:
        expand(f"{trim_galore_dir}/{{sample}}/fastqc.txt", sample=sample_names)

    output:
        f"{trim_galore_dir}/multiqc.html"

    singularity:
        f"{container_dir}/{config['preprocess_image']}"

    params:
        slurm_log_dir = f"{analysis_dir}/logs/trim/slurm/"

    shell:
        "multiqc {input} -o {output}"
