# CHANGELOG
All pushed to the repository must be noted here

## v0.11.0 [2020-07-15]
### Added
- r-markdown for ballgown, deseq2, drimseq

## v0.10.0 [2020-07-08]
### Added
- Get list of genes with frequent fusions in treated samples


## v0.9.1 [2020-06-02]
### Fixed
- STAR command call

## v0.9.0 [2020-06-01]
### Added
- Arriba: Arriba is a tool developed to extract possible gene fusion from STAR alignment
- BEDGraphs from bam-files, for later visualisation.

### Fixed
- Fix the STAR aligement to output the chimeric reads into the output alignment file. This is done so that arriba can be used on the output.

## v0.8.0 [2020-05-27]
### Added
- salmon

## v0.7.0 [2020-05-20]
### Added
- stringtie

## v0.6.0 [2020-05-05]
### Added
- gene count rule
- HTSeq container

### Fixed
- star command
- config indentation
- bam index threads

## v.0.5.0 [2020-04-30]
### Added
- bam index rule
- cluster status script

## v.0.4.0 [2020-04-28]
### Added
- Better structure for slurm logs

## v.0.3.0 [2020-04-28]
### Added
- Added alignment rule
- added multiqc rule for trimmed data

## v.0.2.0 [2020-04-27]
### Added
- trimmed_multiqc rule

## v0.1.1 [2020-04-27]
### Fixed
- pipeline now find the sample names in samples.tsv with first row specifying sample names

## v0.1.0 [2020-04-26]
### Added
- support for slurm configuration

## Singularity fix (v0.0.1) [2020-04-24]
### fixed
- Add singularity images to star_index, and trim rules

### Added
- add bumpversion config file to repo
- add VERSION file

## Add Rules [2020-04-24]
### Added
- create preprocessing rules
- create container rules

## Add rules [2020-04-22]
### Added
- create download rules
- create indexing rule

## Set up structure for  project [2020-04-22]
### Added
- touch rule files
- create config.yaml
- create README.md
- create CHANGELOG.md
- create Snakefile
