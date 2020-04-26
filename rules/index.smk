rule star_index:
    input:
        genome_fasta = expand_path(reference_dir, config["genome_file"]),
        annotation_gtf = expand_path(reference_dir, config['transcripts_file'])

    output:
        genome_dir = expand_path(reference_dir, "star_genome")

    params:
        index_read_length = str(config["read_length"]),
        n_cores = cluster_config["star_index"]["n"],
        slurm_log_dir = f"{analysis_dir}/logs/star_index/slurm/"

    singularity:
        f"{container_dir}/{config['star_image']}"

    shell:
        "mkdir {output.genome_dir} && "
        "STAR --runThreadN {params.n_cores} "
        "--runMode genomeGenerate "
        "--genomeDir {output.genome_dir} "
        "--genomeFastaFiles {input.genome_fasta} "
        "--sjdbGTFfile {input.annotation_gtf} "
        "--sjdbOverhang {params.index_read_length}"
