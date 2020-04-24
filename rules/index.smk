
rule star_index:

    singularity:
        f"{container_dir}/{config['star_image']}"
        
    input:
        genome_fasta = expand_path(reference_dir, config["genome_file"]),
        annotation_gtf = expand_path(reference_dir, config['transcripts_file'])

    output:
        genome_dir = expand_path(reference_dir, "star_genome")

    params:
        index_threads = str(config["star_index_threads"]),
        index_read_length = str(config["read_length"])

    shell:
        "mkdir {output.genome_dir} && "
        "STAR --runThreadN {params.index_threads} "
        "--runMode genomeGenerate "
        "--genomeDir {output.genome_dir} "
        "--genomeFastaFiles {input.genome_fasta} "
        "--sjdbGTFfile {input.annotation_gtf} "
        "--sjdbOverhang {params.index_read_length}"
