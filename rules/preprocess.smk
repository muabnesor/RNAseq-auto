def get_fastqs(wildcards):
    fastq_first = f"{fastq_dict[wildcards.sample]['fastq_first']}"
    fastq_second = f"{fastq_dict[wildcards.sample]['fastq_second']}"
    return [fastq_first, fastq_second]

rule trim:

    singularity:
        f"{container_dir}/{config['preprocess_image']}"

    input:
        get_fastqs

    output:
        f"{trim_galore_dir}/{{sample}}"

    shell:
        "trim_galore --illumina --paired --fastqc -o {output} "
        "{input}"
