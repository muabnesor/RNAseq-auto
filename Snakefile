from pathlib import Path
import pandas as pd

# Include common python functions
include: "rules/common.smk"

configfile: "config.yaml"

# Set base directories
data_dir = str(Path(config["data_dir"]))
analysis_dir = str(Path(config["analysis_dir"]))
reference_dir = str(Path(config["reference_dir"]))
container_dir = str(Path(config["container_dir"]))
log_dir = str(Path(config["log_dir"]))

# Set analysis subdirectories
trim_galore_dir = expand_path(analysis_dir, "trim_galore")
align_dir = expand_path(analysis_dir, "align")

# get sample data
sample_data = config.get("sample_data") or expand_path(data_dir, "samples.tsv")
coldata = pd.read_table(sample_data)

sample_names = list(coldata.loc[:, "sample"])

# Parse base directory
fastq_dict = get_fastq_dict(base_dir=Path(data_dir),
                            fastq1_suffix=config["fastq1_suffix"],
                            fastq2_suffix=config["fastq2_suffix"],
                            sample_names=sample_names)

report: "report/workflow.rst"

# Downloading references
include: "rules/download.smk"

# Building singularity containers
include: "rules/containers.smk"

# Preprocessing fastq-files
include: "rules/preprocess.smk"

# Indexing reference genome
include: "rules/index.smk"

# Aligning reads with STAR
include: "rules/align.smk"

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
