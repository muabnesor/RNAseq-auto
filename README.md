# RNAseq-auto
This is a fully automated and scalable pipeline written for the snakemake workflow
Manager.

All the external tools used by the pipeline exists in docker containers handled
by the developer by default. Workflows for downloading containers, and reference
files are included, but the user may specify these in the config.yaml file.

The pipeline supports usage on a HPC with the slurm scheduler and singularity
installed.
