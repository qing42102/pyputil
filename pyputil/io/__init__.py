import gzip

from ruamel.yaml import YAML
yaml = YAML(typ='safe')


def zopen(filename: str, *args):
    if filename.endswith(".gz"):
        return gzip.open(filename, *args)
    else:
        return open(filename, *args)
