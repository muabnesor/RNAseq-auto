rule download_all:
    input:
        genome_path = expand_path(reference_dir, config["references"]["genome_file"]),
        transcripts_path = expand_path(reference_dir, config["references"]["transcripts_file"])

rule download_genome:
    params:
        genome_url = config["references"]["genome_url"]
    output:
        genome_path = expand_path(reference_dir, config["references"]["genome_file"])
    shell:
        "wget -O {output.genome_path}.gz {params.genome_url} && "
        "gunzip {output.genome_path}.gz"

rule download_transcripts:
    params:
        transcripts_url = config["references"]["transcripts_url"]
    output:
        transcripts_path = expand_path(reference_dir, config["references"]["transcripts_file"])
    shell:
        "wget -O {output.transcripts_path}.gz {params.transcripts_url} && "
        "gunzip {output.transcripts_path}.gz"
