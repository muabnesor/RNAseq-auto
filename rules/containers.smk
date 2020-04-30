rule build_all_containers:
    input:
        star_sif = expand_path(container_dir, config["star_image"]),
        preprocess_sif = expand_path(container_dir, config["preprocess_image"])

rule build_star:
    params:
        image_url = config['star_container_url']
    output:
        sif_file = expand_path(container_dir, config["star_image"])
    shell:
        "sudo singularity build {output.sif_file} {params.image_url}"

rule build_preprocess:
    params:
        image_url = config['preprocess_container_url']
    output:
        sif_file = expand_path(container_dir, config["preprocess_image"])
    shell:
        "sudo singularity build {output.sif_file} {params.image_url}"

rule build_samtools:
    params:
        image_url = config['samtools_container_url']
    output:
        sif_file = expand_path(container_dir, config["samtools_image"])
    shell:
        "sudo singularity build {output.sif_file} {params.image_url}"
