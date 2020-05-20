from pathlib import Path

def expand_path(dir: str, extension: str):
    return str(Path(dir).joinpath(extension))

def get_fastq_dict(base_dir: Path, fastq1_suffix: str, fastq2_suffix: str, coldata) -> dict:
    fastq_dict = dict()
    sample_names = list(coldata.loc[:, "sample"])
    for sub_dir in base_dir.iterdir():
        if not sub_dir.is_dir():
            continue
        sample_name = sub_dir.name.split("_")[-1]
        if sample_name not in sample_names:
            continue
        sample_dict = dict()
        for file in sub_dir.iterdir():
            if not file.is_file():
                continue
            if file.name.split("_")[0] != sample_name:
                continue
            if file.name.endswith(fastq1_suffix):
                sample_dict["fastq_first"] = str(file)
            if file.name.endswith(fastq2_suffix):
                sample_dict["fastq_second"] = str(file)
            if len(sample_dict) == 2:
                fastq_dict[sample_name] = sample_dict
                break

    return fastq_dict
