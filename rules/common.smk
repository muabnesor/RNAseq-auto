from pathlib import Path

def expand_path(dir: str, extension: str):
    return str(Path(dir).joinpath(extension))

def get_fastq_dict(base_dir: Path, fastq1_suffix: str, fastq2_suffix: str) -> dict:
    fastq_dict = dict()
    for sub_dir in base_dir.iterdir():
        if not sub_dir.is_dir():
            continue
        sample_name = sub_dir.name
        sample_dict = dict()
        for file in sub_dir.iterdir():
            if not file.is_file():
                continue
            if file.name.endswith(fastq1_suffix):
                sample_dict["fastq_first"] = str(file)
            if file.name.endswith(fastq2_suffix):
                sample_dict["fastq_second"] = str(file)
            if len(sample_dict) == 2:
                fastq_dict[sample_name] = sample_dict
                break
    return fastq_dict
