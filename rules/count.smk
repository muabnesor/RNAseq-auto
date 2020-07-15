
rule count_all:
    input:
        expand(f"{count_dir}/{{sample}}_counts.tsv", sample=sample_names)

rule count_genes:
    input:
        bam_file = f"{align_dir}/{{sample}}.bam",
        annotation_gtf = expand_path(reference_dir, config["references"]["transcripts_file"])

    output:
        count_file = f"{count_dir}/{{sample}}_counts.tsv"

    singularity:
        f"{container_dir}/{config['containers']['htseq_image']}"

    params:
        slurm_log_dir = f"{str(slurm_logdir_count)}"

    threads: 1

    shell:
        "htseq-count -s reverse -f bam -r pos {input.bam_file} {input.annotation_gtf} "
        "> {output.count_file}"

rule build_genes:
    input:
        gtf_file = expand_path(reference_dir, config["references"]["transcripts_file"])

    output:
        gene_db = expand_path(reference_dir, "genes.sqllite")

    singularity:
        f"{container_dir}/{config['containers']['R_image']}"

    params:
        slurm_log_dir = f"{str(slurm_logdir_count)}"

    threads: 4

    script: "../scripts/ensembl-db.R"


rule htseq_deseq2:
    input:
        gene_db = expand_path(reference_dir, "genes.sqllite"),
        count_files = expand(f"{count_dir}/{{sample}}_counts.tsv", sample=sample_names),
        sample_data = f"{sample_data}"

    output:
        f"{results_dir}/htseq-deseq2/htseq-deseq2.html"

    singularity: f"{container_dir}/{config['containers']['R_image']}"

    params:
        slurm_log_dir = f"{str(slurm_logdir_count)}",
        count_dir = f"{count_dir}"

    threads: 1

    script: "../scripts/htseq-deseq2.Rmd"
