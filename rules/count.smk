
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

    shell:
        "htseq-count -s reverse -f bam -r pos {input.bam_file} {input.annotation_gtf} "
        "> {output.count_file}"
