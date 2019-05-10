import gzip
from pathlib import Path

from ruamel.yaml import YAML
yaml = YAML(typ='safe')


def zopen(path: Path, *args):
    if path.suffix == ".gz":
        return gzip.open(path, *args)
    else:
        return path.open(*args)


def write_bytes(filename: str, content: bytes):
    with open(filename, 'wb') as f:
        f.write(content)
