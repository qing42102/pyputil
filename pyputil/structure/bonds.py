import itertools

import numpy as np
import scipy.spatial
from pymatgen import Structure


def bond_to_positions(bond: np.array, lattice, cart_coords):
    # lattice input should be each lattice vector as a row
    offset = np.dot(lattice.T, bond[:3])
    pos_a = cart_coords[bond[3]]
    pos_b = cart_coords[bond[4]] + offset.T
    return pos_a, pos_b


def bonds_to_positions(bonds: np.array, lattice, cart_coords):
    pos_a, pos_b = bond_to_positions(bonds.T, lattice, cart_coords)
    return np.hstack((pos_a, pos_b)).reshape((-1, 2, 3))


def calculate_bonds(structure: Structure, cutoff=1.6):
    # just make a huge system to avoid periodicity edge cases
    lattice = structure.lattice.matrix
    vecs = [lattice[1], lattice[2], lattice[0]]

    def plane_dist(la, lb, lc):
        normal = np.cross(la, lb)
        normal /= np.linalg.norm(normal)
        return abs(np.dot(normal, lc))

    image_ranges = []
    for idx in range(3):
        dist = plane_dist(*np.roll(vecs, idx, axis=0))
        image_ranges.append(
            np.arange(max(int(np.ceil(cutoff / dist)), 2)))

    # build cartesian offsets from the multiples of each lattice vector
    images = np.array(list(itertools.product(*image_ranges)))
    ref_coords = structure.cart_coords

    bonds = []
    for a, b, c in images:
        offset = np.dot(lattice.T, [a, b, c])
        image_coords = ref_coords + offset

        distances = scipy.spatial.distance.cdist(
            ref_coords, image_coords)

        new_bonds = np.transpose((distances < cutoff).nonzero())
        image_offsets = np.repeat([[a, b, c]], new_bonds.shape[0], axis=0)
        new_bonds = np.hstack((image_offsets, new_bonds))

        if a == 0 and b == 0 and c == 0:
            # don't double count
            bonds.append(new_bonds[new_bonds[:, 3] < new_bonds[:, 4]])
        else:
            bonds.append(new_bonds)

    return np.vstack(bonds)


def calculate_bond_list(structure: Structure, cutoff=1.6):
    bonds = calculate_bonds(structure=structure, cutoff=cutoff)

    # double count bonds
    reverse_bonds = np.column_stack((
        -1 * bonds[:, :3],
        bonds[:, 4],
        bonds[:, 3],
    ))
    # don't double count bonds that connect to the same atom though
    reverse_bonds = reverse_bonds[reverse_bonds[:, 4] != reverse_bonds[:, 3]]
    bonds = np.vstack((bonds, reverse_bonds))

    # sort by first atom index in bond
    bonds = bonds[np.argsort(bonds[:, 3])]

    # split bonds where the indices change to make a list of bond arrays
    split_locations = np.diff(bonds[:, 3]).nonzero()[0] + 1
    return np.split(bonds, split_locations)

