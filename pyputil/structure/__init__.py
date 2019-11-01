import typing as tp
import numpy as np
from pymatgen import Structure, Element, Lattice

from pyputil.structure.constants import DEFAULT_CH_DIST, DEFAULT_CC_DIST
from pyputil.structure.bonds import calculate_bond_list, bonds_to_positions


def get_structure_dims(arg: tp.Union[str, Structure]):
    if type(arg) == str:
        structure = Structure.from_file(arg)
    else:
        structure = arg

    cc = structure.cart_coords
    avg_c = np.average(cc, axis=0)
    min_c = cc.min(axis=1)
    max_c = cc.max(axis=1)
    radius = np.linalg.norm(cc - avg_c, axis=0).max()
    return {
        "min": min_c,
        "max": max_c,
        "avg": avg_c,
        "radius": radius,
    }


def calc_hydrogen_loc(
        a: np.array,
        delta_ab: np.array,
        delta_ac: np.array,
        dist: float,
) -> np.array:
    direction = -(delta_ab + delta_ac) / 2
    direction *= dist / np.linalg.norm(direction)
    return a + direction


def add_hydrogen(
        structure: Structure,
        cutoff: float,
        dist: float,
):
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
            coords=calc_hydrogen_loc(a, b - a, c - a, dist),
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


def remove_random_carbon(structure: Structure) -> Structure:
    structure = structure.copy()

    carbon_indices = [idx for idx, s in enumerate(structure.species) if s == Element.C]
    to_remove = np.random.randint(len(carbon_indices))

    structure.remove_sites([to_remove])
    return structure


# note: adds hydrogen to the edges at the end
def add_vacancy_defects(
        pristine: Structure,
        defect_percent: float,
        cutoff: float = 1.1 * DEFAULT_CC_DIST,
) -> Structure:
    structure = pristine.copy()

    def should_add_defect():
        return np.random.rand() < defect_percent

    # basically just go through and try to remove carbon atoms
    structure.remove_sites([
        idx for idx, s in enumerate(structure.species)
        if s == Element.C and should_add_defect()
    ])

    # remove lone atoms (based on # of carbon neighbors)
    while True:
        species = structure.species
        to_remove = []
        for idx, bonds in enumerate(calculate_bond_list(structure=structure, cutoff=cutoff)):
            num_carbon_neighbors = np.sum(
                1 for b in bonds[:, 4]
                if species[b] == Element.C
            )
            if num_carbon_neighbors == 0:
                to_remove.append(idx)

        if len(to_remove) == 0:
            break
        else:
            structure.remove_sites(to_remove)

    # make sure to add hydrogen back to the edges when done though
    add_hydrogen(structure, cutoff, dist=DEFAULT_CH_DIST)

    return structure
