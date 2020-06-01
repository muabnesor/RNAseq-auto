
rule fuse_all:
    input:
        expand(f"{fusion_dir}/{{sample}}.tsv", sample=sample_names)

rule fuse:
    input:
        bam_file = f"{align_dir}/{{sample}}.bam",
        genome_fasta = expand_path(reference_dir, config["references"]["genome_file"]),
        annotation_gtf = expand_path(reference_dir, config["references"]["transcripts_file"])

    output:
        fusions = f"{fusion_dir}/{{sample}}.tsv",

    params:
        slurm_log_dir = f"{str(slurm_logdir_arriba)}"

    singularity:
        f"{container_dir}/{config['containers']['arriba_image']}"

    threads: 1

    shell:
        "arriba -x {input.bam_file} "
        "-g {input.annotation_gtf} "
        "-a {input.genome_fasta} "
        "-o {output.fusions}"
