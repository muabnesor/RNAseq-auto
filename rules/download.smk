rule download_all:
    input:
        genome_path = expand_path(reference_dir, config["genome_file"]),
        transcripts_path = expand_path(reference_dir, config["transcripts_file"])

rule download_genome:
    params:
        genome_url = config["genome_url"]
    output:
        genome_path = expand_path(reference_dir, config["genome_file"])
    shell:
        "wget -O {output.genome_path}.gz {params.genome_url} && "
        "gunzip {output.genome_path}.gz"

rule download_transcripts:
    params:
        transcripts_url = config["transcripts_url"]
    output:
        transcripts_path = expand_path(reference_dir, config["transcripts_file"])
    shell:
        "wget -O {output.transcripts_path}.gz {params.transcripts_url} && "
        "gunzip {output.transcripts_path}.gz"
