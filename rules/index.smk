rule star_index:
    input:
        genome_fasta = expand_path(reference_dir, config["references"]["genome_file"]),
        annotation_gtf = expand_path(reference_dir, config["references"]["transcripts_file"])

    output:
        genome_dir = directory(expand_path(reference_dir, "star_genome"))

    params:
        index_read_length = str(config["index"]["read_length"]),
        slurm_log_dir = f"{str(slurm_logdir_index)}"

    singularity:
        f"{container_dir}/{config['containers']['star_image']}"

    threads: config["index"]["star_index_threads"]

    shell:
        "mkdir {output.genome_dir} && "
        "STAR --runThreadN {threads} "
        "--runMode genomeGenerate "
        "--genomeDir {output.genome_dir} "
        "--genomeFastaFiles {input.genome_fasta} "
        "--sjdbGTFfile {input.annotation_gtf} "
        "--sjdbOverhang {params.index_read_length}"

rule salmon_index:
    input:
        cdna_fasta = expand_path(reference_dir, config["references"]["cdna_file"])

    output:
        index_dir = directory(expand_path(reference_dir, "salmon_index"))

    params:
        slurm_log_dir = f"{str(slurm_logdir_index)}"

    singularity:
        f"{container_dir}/{config['containers']['salmon_image']}"

    threads: config["index"]["salmon_index_threads"]

    shell:
        "salmon index -p {threads} -t {input.cdna_fasta} -i {output.index_dir}"
