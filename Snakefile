from pathlib import Path

# Include common python functions
include: "rules/common.smk"

configfile: "config.yaml"

working_dir = str(Path(config["working_dir"]))
reference_dir = str(Path(config["reference_dir"]))
containers_dir = str(Path(config["container_dir"]))

align_dir = expand_path(working_dir, "align")
results_dir = expand_path(working_dir, "results")

rule all:
    input:
        expand_path(results_dir, "results.txt")

report: "report/workflow.rst"

# Downloading references
include: "rules/download.smk"

# Building singularity containers
#include: "rules/containers.smk"

# Preprocessing fastq-files
#include: "rules/preprocess"

# Indexing reference genome
include: "rules/index.smk"

# Aligning reads with STAR
#include: "rules/align.smk"

# Quality check of aligned reads
#include: "rules/qc.smk"

# Build count matrix
#include: "rules/count.smk"

# Check splice junctions
#include: "rules/splice.smk"

# Perform DE-analysis with DESeq2
#include: "rules/de.smk"

# Compile results into report
#include: "rules/compile.smk"
