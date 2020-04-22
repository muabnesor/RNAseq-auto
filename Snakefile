from pathlib import Path

configfile: "config.yaml"

working_dir = str(Path(config["working_dir"]))
reference_dir = str(Path(config["reference_dir"]))
containers_dir = str(Path(config["container_dir"]))

def expand_path(dir: str, extension: str):
    return str(Path(dir).joinpath(extension))

genome_fa_path = expand_path(reference_dir, config["genome_file"])
transcripts_gtf_path = expand_path(reference_dir, config["transcript_file"])

align_dir = expand_path(working_dir, "align")
stats_dir = expand_path(working_dir, "stats")

rule all:
    input:
        expand_path(stats_dir, "result.txt")

singularity: "docker://continuumio/miniconda3"

report: "report/workflow.rst"

# Downloading references
#include: "rules/download.smk"

# Building singularity containers
#include: "rules/containers.smk"

# Preprocessing fastq-files
#include: "rules/preprocess"

# Indexing reference genome
include: "rules/index.smk"

# Aligning reads with STAR
include: "rules/align.smk"

# Quality check of aligned reads
include: "rules/qc.smk"

# Build count matrix
include: "rules/count.smk"

# Check splice junctions
include: "rules/splice.smk"

# Perform DE-analysis with DESeq2
include: "rules/de.smk"

# Compile results into report
include: "rules/compile.smk"
