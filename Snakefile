from pathlib import Path
import pandas as pd

# Include common python functions
include: "rules/common.smk"

configfile: "config.yaml"

# Set base directories
data_dir = str(Path(config["data_dir"]))
analysis_dir = str(Path(config["analysis_dir"]))
reference_dir = str(Path(config["references"]["reference_dir"]))
container_dir = str(Path(config["containers"]["container_dir"]))
log_dir = str(Path(config["log_dir"]))

# Create slurm log dirs
slurm_logdir = Path(log_dir).joinpath("slurm")

slurm_logdir_preprocess = Path(slurm_logdir).joinpath("preprocess")
slurm_logdir_preprocess.mkdir(parents=True, exist_ok=True)

slurm_logdir_index = Path(slurm_logdir).joinpath("index")
slurm_logdir_index.mkdir(parents=True, exist_ok=True)

slurm_logdir_align = Path(slurm_logdir).joinpath("align")
slurm_logdir_align.mkdir(parents=True, exist_ok=True)

slurm_logdir_count = Path(slurm_logdir).joinpath("count")
slurm_logdir_count.mkdir(parents=True, exist_ok=True)

slurm_logdir_transcripts = Path(slurm_logdir).joinpath("transcripts")
slurm_logdir_transcripts.mkdir(parents=True, exist_ok=True)

slurm_logdir_pseudo_map = Path(slurm_logdir).joinpath("pseudo_map")
slurm_logdir_pseudo_map.mkdir(parents=True, exist_ok=True)

slurm_logdir_arriba = Path(slurm_logdir).joinpath("arriba")
slur_logdir_arriba.mkdir(parents=True, exist_ok=True)

# Set analysis subdirectories
trim_galore_dir = expand_path(analysis_dir, "trim_galore")
align_dir = expand_path(analysis_dir, "align")
fusion_dir = expand_path(analysis_dir, "fusion")
count_dir = expand_path(analysis_dir, "count")
transcripts_dir = expand_path(analysis_dir, "transcripts")
pseudo_map_dir = expand_path(analysis_dir, "pseudo_map")

# get sample data
sample_data = config.get("sample_data") or expand_path(data_dir, "samples.tsv")
coldata = pd.read_table(sample_data)

sample_names = list(coldata.loc[:, "sample"])

# Parse base directory
fastq_dict = get_fastq_dict(base_dir=Path(data_dir),
                            fastq1_suffix=config["fastq1_suffix"],
                            fastq2_suffix=config["fastq2_suffix"],
                            coldata=coldata)

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

# Build count matrix
include: "rules/count.smk"

# Find transcripts
include: "rules/transcribe.smk"

# Perform pseudo mapping with salmon
include: "rules/pseudo_map.smk"

# check for fusions with arriba
include: "rules/fuse.smk"

# Perform DE-analysis with DESeq2
#include: "rules/de.smk"

# Compile results into report
#include: "rules/compile.smk"
