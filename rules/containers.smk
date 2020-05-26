rule build_all_containers:
    input:
        star_sif = expand_path(container_dir, config["containers"]["star_image"]),
        preprocess_sif = expand_path(container_dir, config["containers"]["preprocess_image"])

rule build_star:
    params:
        image_url = config["containers"]['star_container_url']
    output:
        sif_file = expand_path(container_dir, config["containers"]["star_image"])
    shell:
        "sudo singularity build {output.sif_file} {params.image_url}"

rule build_stringtie:
    params:
        image_url = config["containers"]['stringtie_container_url']
    output:
        sif_file = expand_path(container_dir, config["containers"]["stringtie_image"])
    shell:
        "sudo singularity build {output.sif_file} {params.image_url}"

rule build_preprocess:
    params:
        image_url = config["containers"]['preprocess_container_url']
    output:
        sif_file = expand_path(container_dir, config["containers"]["preprocess_image"])
    shell:
        "sudo singularity build {output.sif_file} {params.image_url}"

rule build_samtools:
    params:
        image_url = config["containers"]['samtools_container_url']
    output:
        sif_file = expand_path(container_dir, config["containers"]["samtools_image"])
    shell:
        "sudo singularity build {output.sif_file} {params.image_url}"

rule build_htseq:
    params:
        image_url = config["containers"]['htseq_container_url']
    output:
        sif_file = expand_path(container_dir, config["containers"]["htseq_image"])
    shell:
        "sudo singularity build {output.sif_file} {params.image_url}"

rule build_salmon:
    params:
        image_url = config["containers"]['salmon_container_url']
    output:
        sif_file = expand_path(container_dir, config["containers"]["salmon_image"])
    shell:
        "sudo singularity build {output.sif_file} {params.image_url}"
