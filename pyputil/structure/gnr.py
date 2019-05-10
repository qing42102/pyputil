import typing as tp

import numpy as np
from pymatgen import Structure, Lattice, Element

from pyputil.structure.bonds import bonds_to_positions, \
    calculate_bond_list


def calc_hydrogen_loc(
        a: np.array,
        delta_ab: np.array,
        delta_ac: np.array,
        dist: float = 0.75
) -> np.array:
    direction = -(delta_ab + delta_ac) / 2
    direction *= dist / np.linalg.norm(direction)
    return a + direction


def add_hydrogen(structure: Structure, cutoff: float = 1.05):
    all_bonds = calculate_bond_list(structure=structure, cutoff=cutoff)

    lattice = np.copy(structure.lattice.matrix)
    coords = np.copy(structure.cart_coords)

    for bonds in all_bonds:
        # only add hydrogen to atoms with two neighbors
        if len(bonds) != 2:
            continue

        bond_pos = bonds_to_positions(bonds, lattice, coords)
        a, b = bond_pos[0]
        _, c = bond_pos[1]
        structure.append(
            species=Element.H,
            coords=calc_hydrogen_loc(a, b - a, c - a),
            coords_are_cartesian=True,
        )


def add_vacuum_sep(
        structure: Structure,
        vx: tp.Optional[float] = None,
        vy: tp.Optional[float] = None,
        vz: tp.Optional[float] = None,
) -> Structure:
    pos_min = np.min(structure.cart_coords, axis=0)
    pos_max = np.max(structure.cart_coords, axis=0)
    pos_range = pos_max - pos_min

    new_lattice = np.copy(structure.lattice.matrix)
    for i, vsep in enumerate([vx, vy, vz]):
        if vsep is None:
            continue

        # current lattice vector should be orthogonal
        assert new_lattice[i, i] != 0
        assert new_lattice[i, (i + 1) % 3] == 0
        assert new_lattice[i, (i + 2) % 3] == 0

        new_lattice[i, i] = vsep + pos_range[i]

    return Structure(
        lattice=Lattice(matrix=new_lattice),
        species=structure.species,
        coords=structure.cart_coords,
        coords_are_cartesian=True,
    )


def center_structure(structure: Structure) -> Structure:
    structure = structure.copy()
    avg_pos = np.average(structure.frac_coords, axis=0)
    structure.translate_sites(
        np.arange(0, structure.num_sites),
        np.array([0.5, 0.5, 0.5]) - avg_pos,
        frac_coords=True,
        to_unit_cell=True,
    )
    return structure


def agnr_unit_cell(vacuum_sep: float = 10) -> Structure:
    return Structure(
        lattice=Lattice([
            [3, 0, 0],
            [0, np.sqrt(3), 0],
            [0, 0, vacuum_sep],
        ]),
        species=[
            Element.C,
            Element.C,
            Element.C,
            Element.C,
        ],
        coords=np.array([
            [0.5, 0, 0],
            [1.5, 0, 0],
            [0.0, np.sqrt(3) / 2, 0],
            [2.0, np.sqrt(3) / 2, 0],
        ]) + np.array([0.5, np.sqrt(3) / 4, vacuum_sep / 2]),
        coords_are_cartesian=True,
    )


def generate_periodic_agnr(
        width_n: int,
        bond_dist: float = 1.44,
        vacuum_sep: float = 15,
        hydrogen: bool = True,
):
    assert width_n >= 2

    # put vacuum separation in terms of the bond distance
    vacuum_sep /= bond_dist

    cell = agnr_unit_cell(vacuum_sep=vacuum_sep)
    if width_n % 2 == 0:
        cell.make_supercell([1, width_n // 2, 1])
    else:
        cell.make_supercell([1, width_n // 2 + 1, 1])
        # need to remove the two topmost atoms
        cell.remove_sites(
            np.argsort(cell.cart_coords[:, 1])[-2:]
        )

    # add vacuum separation
    cell = add_vacuum_sep(cell, vy=vacuum_sep, vz=vacuum_sep)
    # center new structure
    cell = center_structure(cell)
    if hydrogen:
        add_hydrogen(cell)

    # scale to bond distances
    cell.lattice = Lattice(matrix=cell.lattice.matrix * bond_dist)
    return cell


def _generate_finite_agnr_from_periodic(
        periodic_gnr: Structure,
        length_m: int,
        bond_dist: float,
        vacuum_sep: float,
        hydrogen: bool,
):
    assert length_m >= 2
    cell = periodic_gnr

    if length_m % 2 == 0:
        cell.make_supercell([length_m // 2, 1, 1])
    else:
        num_to_remove = cell.num_sites // 2
        cell.make_supercell([length_m // 2 + 1, 1, 1])
        # remove rightmost atoms
        to_remove = np.argsort(cell.cart_coords[:, 0])[-num_to_remove:]
        cell.remove_sites(to_remove)

    # put vacuum separation in terms of the bond distance
    vacuum_sep /= bond_dist
    cell = add_vacuum_sep(cell, vx=vacuum_sep, vy=vacuum_sep, vz=vacuum_sep)
    # center new structure
    cell = center_structure(cell)

    # kill off any lone carbons or those with only 1 bond
    n_bonds = np.array([len(b) for b in calculate_bond_list(cell, cutoff=1.1)])
    cell.remove_sites(np.where(n_bonds <= 1)[0])

    if hydrogen:
        add_hydrogen(cell)

    # scale to bond distances
    cell.lattice = Lattice(matrix=cell.lattice.matrix * bond_dist)
    return cell


def generate_finite_agnr(
        width_n: int,
        length_m: tp.Union[tp.Iterable[int], int],
        bond_dist: float = 1.44,
        vacuum_sep: float = 15,
        hydrogen: bool = True,
):
    cell = generate_periodic_agnr(
        width_n=width_n, bond_dist=1.0, vacuum_sep=vacuum_sep, hydrogen=False)

    def generate(length: int) -> Structure:
        return _generate_finite_agnr_from_periodic(
            cell.copy(), length, bond_dist=bond_dist, vacuum_sep=vacuum_sep, hydrogen=hydrogen)

    try:
        # length_m is hopefully iterable
        return map(generate, length_m)
    except TypeError:
        return generate(length_m)
