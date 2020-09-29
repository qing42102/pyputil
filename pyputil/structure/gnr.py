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


def edge_types(cell: Structure, atom_list = None, combinations = None):
    #Get all the bonds
    all_bonds = calculate_bond_list(structure=cell)

    #Get the atoms with only 2 bonds
    atoms_2bonds = []
    atoms_2bonds_index = []
    for index, atoms in enumerate(all_bonds):
        if len(atoms) == 2:
            atoms_2bonds.append(atoms)
            atoms_2bonds_index.append(index)

    #To check whether the graphene structure is broken if there is an atom with only 1 bond
    def broken_structure(cell: Structure):
        all_bonds = calculate_bond_list(structure=cell)
        for atoms in all_bonds:
            if len(atoms) == 1:
                return True
        return False

    #Store the combinations of atoms indices list
    if combinations == None:
        combinations = set()

    if atom_list is None:
        atom_list = [i for i in range(len(all_bonds))]

    #Store all possible combinations of the cell with the removed atoms in a list
    structures = []
    for atoms in atoms_2bonds:
        for bonds in atoms:
            #Go through the atoms with only 2 bonds and check if its neighbor has 2 bonds. 
            if bonds[4] in atoms_2bonds_index:
                atom_list_temp = atom_list.copy()
                
                #Remove the atoms from indices list
                if bonds[3] > bonds[4]:
                    del atom_list_temp[bonds[3]]
                    del atom_list_temp[bonds[4]]
                else:
                    del atom_list_temp[bonds[4]]
                    del atom_list_temp[bonds[3]]
                
                #Remove the atoms from structure
                cell_temp = cell.copy()
                cell_temp.remove_sites([bonds[3], bonds[4]])
                
                #Check whether there are duplicates and if the structure is broken
                if tuple(atom_list_temp) not in combinations and broken_structure(cell_temp) == False:
                    combinations.add(tuple(atom_list_temp))
                    structures.append(cell_temp)
                    structures.extend(edge_types(cell_temp, atom_list_temp, combinations))

    return structures

'''
def edge_types(cell: Structure, combinations = None):
    #Get all the bonds
    all_bonds = calculate_bond_list(structure=cell)

    #Get the atoms with only 2 bonds
    atoms_2bonds = []
    atoms_2bonds_index = []
    for index, atoms in enumerate(all_bonds):
        if len(atoms) == 2:
            atoms_2bonds.append(atoms)
            atoms_2bonds_index.append(index)

    #To check whether the graphene structure is broken if there is an atom with only 1 bond
    def broken_structure(cell: Structure):
        all_bonds = calculate_bond_list(structure=cell)
        for atoms in all_bonds:
            if len(atoms) == 1:
                return True
        return False

    #Store all possible combinations of the cell with the removed atoms in a list
    if combinations == None:
        combinations = []

    for atoms in atoms_2bonds:
        for bonds in atoms:
            #Go through the atoms with only 2 bonds and check if its neighbor has 2 bonds. 
            if bonds[4] in atoms_2bonds_index:
                cell_temp = cell.copy()
                cell_temp.remove_sites([bonds[3], bonds[4]])
                #Check whether there are duplicates and if the structure is broken
                if broken_structure(cell_temp) == False and cell_temp not in combinations:
                    combinations.append(cell_temp)
                    edge_types(cell_temp, combinations)
    
    # #If there are no combinations left, end the recursion
    # if len(combinations) == 0:
    #     return combinations
    # else:
    #     #Recurse for every combinations
    #     for all_cell in combinations:
    #         combinations += edge_types(all_cell)
    #         #Remove the duplicates as it recurses
    #         combinations = [i for n, i in enumerate(combinations) if i not in combinations[:n]] 

    return combinations
'''

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
    # if hydrogen:
    #     add_hydrogen(cell, cutoff=1.05, dist=DEFAULT_CH_DIST / DEFAULT_CC_DIST)

    # scale to bond distances
    #cell.lattice = Lattice(matrix=cell.lattice.matrix * bond_dist)
    return cell


def generate_periodic_zgnr(
        width_n: int,
        bond_dist: float = DEFAULT_CC_DIST,
        vacuum_sep: float = 15,
        hydrogen: bool = True,
):
    assert width_n > 0
    # the way we count is a bit off, so subtract one here to fix it
    width_n -= 1

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
