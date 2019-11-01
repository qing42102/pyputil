import re
import numpy as np
from pathlib import Path

from pymatgen import Structure, Lattice, Element
from ruamel.yaml import YAML
yaml = YAML(typ='safe')


def zopen(path: Path, *args):
    if path.suffix == ".gz":
        import gzip
        return gzip.open(path, *args)
    else:
        return path.open(*args)


def write_bytes(filename: str, content: bytes):
    with open(filename, 'wb') as f:
        f.write(content)


# note: does not read multiple structures from a single file
def read_xyz(path: str, vacuum: float = 15.0):
    spaces = re.compile(r"\s+")
    lines = Path(path).read_text().split("\n")

    natoms = int(lines[0])
    atoms = [spaces.split(l) for l in lines[2:]]

    types = np.array([a[0] for a in atoms])
    positions = np.array([
        [float(a[1]), float(a[2]), float(a[3])] for a in atoms
    ])
    assert natoms == len(positions)
    assert natoms == len(types)

    minv = positions.min(axis=0)
    maxv = positions.max(axis=0)
    lattice = (maxv - minv) + np.array([vacuum, vacuum, vacuum])

    return Structure(
        lattice=Lattice(np.diag(lattice)),
        species=[Element(t) for t in types],
        coords=positions + lattice / 2,
        coords_are_cartesian=True,
    )
