
def get_fastqs(wildcards):
    if not config["pre_trim_fastq"]:
        fastq_first = f"{fastq_dict[wildcards.sample]['fastq_first']}"
        fastq_second = f"{fastq_dict[wildcards.sample]['fastq_second']}"
        return [fastq_first, fastq_second]
    fastq_first = f"{trim_galore_dir}/{wildcards.sample}/{fastq_dict[wildcards.sample]['fastq_first'].split('/')[-1].replace(config['fastq1_suffix'], config['fastq1_suffix_trimmed'])}"
    fastq_second = f"{trim_galore_dir}/{wildcards.sample}/{fastq_dict[wildcards.sample]['fastq_second'].split('/')[-1].replace(config['fastq2_suffix'], config['fastq2_suffix_trimmed'])}"
    return [fastq_first, fastq_second]

rule star_align_all:
    input:
        expand(f"{align_dir}/{{sample}}.bam.bai", sample=sample_names)

rule star_align:
    input:
        fastqs = get_fastqs,
        reference_dir = expand_path(reference_dir, "star_genome")

    output:
        star_dir = directory(f"{align_dir}/{{sample}}"),
        bam_file = f"{align_dir}/{{sample}}.bam"
    params:
        slurm_log_dir = f"{str(slurm_logdir_align)}"

    singularity:
        f"{container_dir}/{config['containers']['star_image']}"

    threads: config['align']["star_align_threads"]

    shell:
        "mkdir -p {output.star_dir} &&"
        "STAR --genomeDir {input.reference_dir} "
        "--runThreadN {threads} "
        "--readFilesIn {input.fastqs} "
        "--outFileNamePrefix {output.star_dir}/ "
        "--outSAMtype BAM SortedByCoordinate "
        "--outSAMunmapped Within "
        "--outSAMattributes Standard "
        "--readFilesCommand zcat "
        "--twopassMode Basic "
        "--chimSegmentMin 20 && "
        "mv {output.star_dir}/Aligned.sortedByCoord.out.bam {output.bam_file}"

rule bam_index:
    input:
        bam_file = f"{align_dir}/{{sample}}.bam"
    output:
        bai_file = f"{align_dir}/{{sample}}.bam.bai"
    params:
        slurm_log_dir = f"{str(slurm_logdir_align)}",
        align_dir = f"{align_dir}"

    threads: config["bam_index"]["bam_index_threads"]
    singularity:
        f"{container_dir}/{config['containers']['samtools_image']}"
    shell:
        "samtools index -@ {threads} {input.bam_file}"
