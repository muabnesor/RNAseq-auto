
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
        expand(f"{align_dir}/{{sample}}", sample=sample_names)

rule star_align:
    input:
        fastqs = get_fastqs,
        reference_dir = expand_path(reference_dir, "star_genome")

    output:
        align_dir = directory(f"{align_dir}/{{sample}}")
    params:
        slurm_log_dir = f"{str(slurm_logdir_align)}"

    singularity:
        f"{container_dir}/{config['star_image']}"

    threads: config["star_align_threads"]

    shell:
        "mkdir -p {output.align_dir} "
        "STAR --genomeDir {input.reference_dir} "
        "--runThreadN {threads} "
        "--readFilesIn {input.fastqs} "
        "--outFileNamePrefix {output.align_dir}/ "
        "--outSamType BAM SortedByCoordinate "
        "--outSAMunmapped Within "
        "--outSAMattributes Standard "
        "--readFilesCommand zcat "
        "--twopassMode Basic "
        "--chimSegmentMin 20"
