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

def agnr_edge_types(cell: Structure, atom_list = None, combinations = None):
    
    #To check whether the graphene structure is broken if there is an atom with only 1 bond
    def broken_structure(cell: Structure):
        all_bonds = calculate_bond_list(structure=cell)
        for atoms in all_bonds:
            if len(atoms) == 1:
                return True
        return False

    #Get all the bonds
    all_bonds = calculate_bond_list(structure=cell)

    #Get the atoms with only 2 bonds
    atoms_2bonds = []
    atoms_2bonds_index = []
    for index, atoms in enumerate(all_bonds):
        if len(atoms) == 2:
            atoms_2bonds.append(atoms)
            atoms_2bonds_index.append(index)

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

def zgnr_edge_types(cell: Structure, atom_list = None, combinations = None):

    def recurse_bonds(current_bonds: np.array, all_bonds: np.array, previous_atom = -1, traverse_atoms = None):
        ''' 
        Get a set of atoms that can be traversed by starting from a atom

        Return a dictionary, traverse_atoms, in which the keys are the atoms traversed and the 
        values are a set of periodic direction to move to that atom
        '''

        if traverse_atoms == None:
            traverse_atoms = dict()
        
        for bonds in current_bonds:
            # Do not move back to the previous atom
            # Do not move in the negative periodic direction
            if bonds[4] != previous_atom and bonds[0] != -1:  
                current_atom = bonds[3]
                neighbor_atom = bonds[4]

                # Create the key if the atom has not been traversed
                if neighbor_atom not in traverse_atoms:
                    traverse_atoms[neighbor_atom] = {bonds[0]}
                    recurse_bonds(all_bonds[neighbor_atom], all_bonds, current_atom, traverse_atoms)
                # Add the new periodic direction. The values can only be 0 or 1
                else:
                    traverse_atoms[neighbor_atom].add(bonds[0])
        
        return traverse_atoms

    def broken_structure(cell: Structure):
        ''' 
        Check whether the graphene structure is broken if there is an atom with only 1 bond or 
        if the structure is no longer periodic
        '''
        try: 
            all_bonds = calculate_bond_list(structure=cell)
        except IndexError: 
            return True

        for atoms in all_bonds:
            if atoms.shape[0] == 1:
                return True

        # Get the list of atoms traversed by starting from atom 0
        atom0 = all_bonds[0][0][3]
        traverse_atoms = recurse_bonds(all_bonds[0], all_bonds)
        traverse_path = set()
        for i in traverse_atoms.values():
            traverse_path.update(i)

        # If there is no path moving back to atom 0 or if no positive periodic direction is 
        # taken, then the structure is broken
        if atom0 not in traverse_atoms.keys() or 1 not in traverse_path \
            or len(traverse_atoms.keys()) != len(all_bonds):
            return True

        return False

    all_bonds = calculate_bond_list(structure=cell)

    # Get the atoms with only 2 bonds
    atoms_2bonds_index = set()
    for index, atoms in enumerate(all_bonds):
        if atoms.shape[0] == 2:
            atoms_2bonds_index.add(index)

    # Store the set of atoms of the structures for search for duplicates
    if combinations == None:
        combinations = set()

    if atom_list is None:
        atom_list = [i for i in range(len(all_bonds))]

    # Store all possible structures with the removed atoms in a list
    structures = []

    for atoms in all_bonds:
        bond1 = atoms[0]
        bond2 = atoms[1]
        bond1_in_2bonds = bond1[4] in atoms_2bonds_index
        bond2_in_2bonds = bond2[4] in atoms_2bonds_index

        cell_temp = cell.copy()
        atom_list_temp = atom_list.copy()

        # Remove the atoms from both the atom_list and the structure
        bond_list = None
        if bond1_in_2bonds and bond2_in_2bonds:
            bond_list = [bond1[3], bond1[4], bond2[4]]
            cell_temp.remove_sites([bond1[3], bond1[4], bond2[4]])
        elif atoms.shape[0] == 3:
            bond3 = atoms[2]
            bond3_in_2bonds = bond3[4] in atoms_2bonds_index
            
            if bond1_in_2bonds and bond3_in_2bonds:
                bond_list = [bond1[3], bond1[4], bond3[4]]
                cell_temp.remove_sites([bond1[3], bond1[4], bond3[4]])
            elif bond2_in_2bonds and bond3_in_2bonds:
                bond_list = [bond2[3], bond2[4], bond3[4]]
                cell_temp.remove_sites([bond2[3], bond2[4], bond3[4]])
        elif bond1[3] in atoms_2bonds_index and bond1_in_2bonds:
            del atom_list_temp[max([bond1[3], bond1[4]])]
            del atom_list_temp[min([bond1[3], bond1[4]])]
            cell_temp.remove_sites([bond1[3], bond1[4]])
        elif bond1[3] in atoms_2bonds_index and bond2_in_2bonds:
            del atom_list_temp[max([bond1[3], bond2[4]])]
            del atom_list_temp[min([bond1[3], bond2[4]])]
            cell_temp.remove_sites([bond1[3], bond2[4]])


        if bond_list != None:
            del atom_list_temp[max(bond_list)]
            bond_list.remove(max(bond_list))
            del atom_list_temp[max(bond_list)]
            del atom_list_temp[min(bond_list)]

        # Check whether there are duplicates and if the structure is broken
        if cell_temp != cell and atom_list_temp != atom_list:
            if broken_structure(cell_temp) == False and tuple(atom_list_temp) not in combinations:
                combinations.add(tuple(atom_list_temp))
                structures.append(cell_temp)
                structures.extend(edge_types(cell_temp, atom_list_temp, combinations))

    return structures


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
    # if hydrogen:
    #     add_hydrogen(cell, cutoff=1.05, dist=DEFAULT_CH_DIST / DEFAULT_CC_DIST)

    # scale to bond distances
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
