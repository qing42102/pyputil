import typing as tp

import numpy as np
from pymatgen import Structure, Lattice, Element

from pyputil.structure import add_hydrogen, add_vacuum_sep, center_structure
from pyputil.structure.constants import DEFAULT_CC_DIST, DEFAULT_CH_DIST
from pyputil.structure.bonds import calculate_bond_list


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


def zgnr_unit_cell(vacuum_sep: float = 10) -> Structure:
    return Structure(
        lattice=Lattice([
            [np.sqrt(3), 0, 0],
            [0, 3, 0],
            [0, 0, vacuum_sep],
        ]),
        species=[
            Element.C,
            Element.C,
            Element.C,
            Element.C,
        ],
        coords=np.array([
            [0, 0.5, 0],
            [0, 1.5, 0],
            [np.sqrt(3) / 2, 0.0, 0],
            [np.sqrt(3) / 2, 2.0, 0],
        ]) + np.array([np.sqrt(3) / 4, 0.5, vacuum_sep / 2]),
        coords_are_cartesian=True,
    )


def generate_periodic_agnr(
        width_n: int,
        bond_dist: float = DEFAULT_CC_DIST,
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
        add_hydrogen(cell, cutoff=1.05, dist=DEFAULT_CH_DIST / DEFAULT_CC_DIST)

    # scale to bond distances
    cell.lattice = Lattice(matrix=cell.lattice.matrix * bond_dist)
    return cell


def generate_periodic_zgnr(
        width_n: int,
        bond_dist: float = DEFAULT_CC_DIST,
        vacuum_sep: float = 15,
        hydrogen: bool = True,
):
    assert width_n > 0

    # put vacuum separation in terms of the bond distance
    vacuum_sep /= bond_dist

    cell = zgnr_unit_cell(vacuum_sep=vacuum_sep)
    cell.make_supercell([1, width_n // 2 + 1, 1])
    if width_n % 2 == 0:
        # need to remove the two topmost atoms
        cell.remove_sites(
            np.argsort(cell.cart_coords[:, 1])[-2:]
        )

    # add vacuum separation
    cell = add_vacuum_sep(cell, vy=vacuum_sep, vz=vacuum_sep)
    # center new structure
    cell = center_structure(cell)
    if hydrogen:
        add_hydrogen(cell, cutoff=1.05, dist=DEFAULT_CH_DIST / DEFAULT_CC_DIST)

    # scale to bond distances
    return Structure(
        lattice=Lattice(matrix=cell.lattice.matrix * bond_dist),
        species=cell.species,
        coords=cell.frac_coords,
    )


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
        add_hydrogen(cell, cutoff=1.05, dist=DEFAULT_CH_DIST / DEFAULT_CC_DIST)

    # scale to bond distances
    cell.lattice = Lattice(matrix=cell.lattice.matrix * bond_dist)
    return cell


def generate_finite_agnr(
        width_n: int,
        length_m: tp.Union[tp.Iterable[int], int],
        bond_dist: float = DEFAULT_CC_DIST,
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
