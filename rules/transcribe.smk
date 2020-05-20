
rule transcribe:
    input:
        bam_file = f"{align_dir}/{{sample}}.bam",
        annotation_gtf = expand_path(reference_dir, config["references"]["transcripts_file"])

    output:
        gtf_file = f"{transcripts_dir}/transcripts/{{sample}}.gtf"

    singularity:
        f"{container_dir}/{config['containers']['stringtie_image']}"

    threads: config["transcribe"]["stringtie_transcribe_threads"]

    params:
        slurm_log_dir = f"{str(slurm_logdir_transcripts)}",
        gtf_dir = f"{transcripts_dir}/transcripts/"

    shell:
        "mkdir -p {params.gtf_dir} && "
        "stringtie {input.bam_file} "
        "-p {threads} "
        "--rf "
        "-G {input.annotation_gtf} "
        "-o {output.gtf_file} "


rule merge_transcripts:
    input:
        gtf_files = expand(f"{transcripts_dir}/transcripts/{{sample}}.gtf", sample=sample_names),
        annotation_gtf = expand_path(reference_dir, config["references"]["transcripts_file"]),
        genome_fasta = expand_path(reference_dir, config["references"]["genome_file"])

    output:
        merged_gtf = f"{transcripts_dir}/merged.gtf"

    singularity:
        f"{container_dir}/{config['containers']['stringtie_image']}"

    threads: config["transcribe"]["stringtie_merge_threads"]

    params:
        transcripts_list = f"{transcripts_dir}/gtf_list.txt",
        slurm_log_dir = f"{str(slurm_logdir_transcripts)}"

    shell:
        "if [ -e {params.transcripts_list} ] ; then rm {params.transcripts_list} ; fi && "
        "for file_name in {input.gtf_files}; do echo ${{file_name}} >> {params.transcripts_list} ; done && "
        "stringtie --merge "
        "-p {threads} "
        "-o {output.merged_gtf} "
        "-G {input.annotation_gtf} "
        "{params.transcripts_list}"


rule count_transcripts:
    input:
        bam_file = f"{align_dir}/{{sample}}.bam",
        merged_gtf = f"{transcripts_dir}/merged.gtf"

    output:
        count_file = f"{transcripts_dir}/counts/{{sample}}.ctab",
        gtf_file = f"{transcripts_dir}/counts/{{sample}}.gtf"


    singularity:
        f"{container_dir}/{config['containers']['stringtie_image']}"

    threads: config["transcribe"]["stringtie_count_threads"]

    params:
        count_dir = f"{transcripts_dir}/counts/",
        slurm_log_dir = f"{str(slurm_logdir_transcripts)}"

    shell:
        "mkdir -p {params.count_dir} && "
        "stringtie {input.bam_file} "
        "-p {threads} "
        "-e -b {output.count_file} "
        "-o {output.gtf_file} "
        "-G {input.merged_gtf} "

rule count_matrix:
    input:
        count_files = expand(f"{transcripts_dir}/counts/{{sample}}.gtf", sample=sample_names)

    output:
        gene_count_matrix = f"{transcripts_dir}/gene_counts.csv",
        transcript_count_matrix = f"{transcripts_dir}/transcrpt_counts.csv",

    singularity:
        f"{container_dir}/{config['containers']['stringtie_image']}"

    threads: config["transcribe"]["stringtie_matrix_threads"]

    params:
        transcripts_list = f"{transcripts_dir}/count_gtf_list.txt",
        slurm_log_dir = f"{str(slurm_logdir_transcripts)}"

    shell:
        "if [ -e {params.transcripts_list} ] ; then rm {params.transcripts_list} ; fi && "
        "for file_name in {input.count_files}; do echo -e ${{file_name%.gtf}}'\t'${{file_name}} >> {params.transcripts_list} ; done && "
        "prepDE.py -i {params.transcripts_list} "
        "-l 151 "
        "-g {output.gene_count_matrix} "
        "-t {output.transcript_count_matrix}"
