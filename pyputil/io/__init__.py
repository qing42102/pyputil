import gzip

from ruamel.yaml import YAML
yaml = YAML(typ='safe')


def zopen(filename: str, *args):
    if filename.endswith(".gz"):
        return gzip.open(filename, *args)
    else:
        return open(filename, *args)


def write_bytes(filename: str, content: bytes):
    with open(filename, 'wb') as f:
        f.write(content)
