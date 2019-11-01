import numpy as np
from pymatgen import Structure, Lattice, Element

from pyputil.structure import add_hydrogen
from pyputil.structure.constants import DEFAULT_CC_DIST, DEFAULT_CH_DIST


def graphene_unit_cell(vacuum_sep: float = 10) -> (np.array, np.array):
    lattice = np.array([
        [np.sqrt(3), 0, 0],
        [np.sqrt(3) / 2, 1.5, 0],
        [0, 0, vacuum_sep],
    ])

    coords = np.array([
        [0, 0, vacuum_sep / 2],
        [np.sqrt(3) / 2, 0.5, vacuum_sep / 2],
    ])

    # put coords in the center of the cell
    coords += np.sum(lattice, axis=0) / 2 - np.average(coords, axis=0)
    return lattice, coords


def _generate_graphene_flake_simple(
        radius: float,
        vacuum_sep: float,
):
    ref_lattice, ref_coords = graphene_unit_cell(vacuum_sep)

    visited = {(0, 0)}
    to_visit = [(0, 0)]

    coords = []
    while len(to_visit) != 0:
        (a, b) = to_visit.pop()

        delta = a * ref_lattice[0] + b * ref_lattice[1]
        coords.extend(ref_coords + delta)

        if np.linalg.norm(delta) > radius + 4.0:
            continue

        for da in [-1, 0, 1]:
            for db in [-1, 0, 1]:
                next = (a + da, b + db)
                if next not in visited:
                    visited.add(next)
                    to_visit.append(next)

    # origin centering (middle of hexagon)
    coords_origin = np.row_stack(coords)
    # center centering (in between A and B sites)
    coords_center = np.copy(coords_origin)
    coords_center[:, :2] -= np.average(ref_coords, axis=0)[:2]
    # a centering (center on A site)
    coords_a = np.copy(coords_origin)
    coords_a[:, :2] -= ref_coords[0, :2]

    def build_structure(c):
        # filter coordinates
        c = c[np.linalg.norm(c[:, :2], axis=1) <= radius, :]
        # build lattice
        min_c = c.min(axis=0)
        max_c = c.max(axis=0)
        lattice = np.diag((max_c - min_c) + np.array([vacuum_sep, vacuum_sep, vacuum_sep]))
        # center flake in the unit cell
        c += np.sum(lattice, axis=0) / 2 - np.average(c, axis=0)
        return Structure(
            lattice=Lattice(lattice),
            species=[Element.C] * c.shape[0],
            coords=c,
            coords_are_cartesian=True,
        )

    return (
        build_structure(coords_origin),
        build_structure(coords_center),
        build_structure(coords_a),
    )


def generate_graphene_flakes(
        radius: float,
        bond_dist: float = DEFAULT_CC_DIST,
        vacuum_sep: float = 15,
) -> [Structure]:
    def process_structure(s):
        # kill off any lone carbons or those with only 1 bond
        from pyputil.structure.bonds import calculate_bond_list
        n_bonds = np.array([len(b) for b in calculate_bond_list(s)])
        s.remove_sites(np.where(n_bonds <= 1)[0])
        # add hydrogen atoms on the edges
        add_hydrogen(s, cutoff=1.05, dist=DEFAULT_CH_DIST / DEFAULT_CC_DIST)
        # scale coordinates to bond distance
        s.lattice = Lattice(matrix=s.lattice.matrix * bond_dist)
        return s

    return [process_structure(s) for s in _generate_graphene_flake_simple(
        radius=radius / bond_dist,
        vacuum_sep=vacuum_sep / bond_dist,
    )]
