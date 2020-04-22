from pathlib import Path

def expand_path(dir: str, extension: str):
    return str(Path(dir).joinpath(extension))
