import math

import numpy as np
from pymatgen import Structure, Element, Lattice

from pyputil.misc import rotation_matrix


def _calc_cnt_param(n: int, m: int):
    """
    get standardized cnt parameters
    """
    # ensure n > m
    if n < m:
        n, m = m, n

    # n (transverse)
    n_t = -(2 * m + n)
    # m (transverse)
    m_t = +(m + 2 * n)

    gcd = math.gcd(n_t, m_t)
    n_t = n_t // gcd
    m_t = m_t // gcd
    return n, m, n_t, m_t


def _cnt_ref_cell() -> (np.array, np.array):
    lattice = np.array([
        [np.sqrt(3), 0, 0],
        [np.sqrt(3) / 2, 1.5, 0],
        [0, 0, 15],
    ])

    coords = np.array([
        [0, 0, 7.5],
        [np.sqrt(3) / 2, 0.5, 7.5],
    ])

    # put coords in the center of the cell
    coords += np.sum(lattice, axis=0) / 2 - np.average(coords, axis=0)

    return lattice, coords


def _gen_flat_cnt(n: int, m: int, n_t: int, m_t: int):
    ref_lattice, ref_coords = _cnt_ref_cell()

    a1 = ref_lattice[0]
    a2 = ref_lattice[1]

    periodic_vec = n * a1 + m * a2
    transverse_vec = n_t * a1 + m_t * a2

    # vectors should be orthogonal
    assert math.fabs(np.dot(periodic_vec, transverse_vec)) < 1e-12

    ref_coords = ref_coords
    cnt_coords = []

    # cover the whole area but try not to double-count (probably an
    # easier way to do this but whatever, it works)
    for m_current in range(max(m, m_t)):
        start = n_t
        stop = n

        if m_t < m and m_t <= m_current:
            start = n + n_t
        elif m < m_t and m <= m_current:
            stop = n + n_t

        for n_current in range(start, stop):
            translation = n_current * a1 + m_current * a2
            cnt_coords.append(ref_coords + translation)

    cnt_coords = np.row_stack(cnt_coords)
    # rotate positions to align to the x axis
    xy_theta = math.atan2(periodic_vec[1], periodic_vec[0])
    xy_rot = rotation_matrix([0, 0, 1], xy_theta)
    cnt_coords = np.dot(cnt_coords, xy_rot)

    # put coordinates in the xz plane (swap y/z since they were
    # originally in the xy plane)
    cnt_coords[:, [1, 2]] = cnt_coords[:, [2, 1]]

    # sort positions by z and x (not sure why this is so difficult...)
    sort_x = np.argsort(cnt_coords[:, 0])
    cnt_coords = cnt_coords[sort_x]
    sort_z = np.argsort(cnt_coords[:, 2], kind='stable')
    cnt_coords = cnt_coords[sort_z]

    # build lattice matrix
    lattice = np.diag([
        np.linalg.norm(periodic_vec),
        ref_lattice[2, 2],
        np.linalg.norm(transverse_vec),
    ])

    return lattice, cnt_coords


def generate_cnt(n: int, m: int, vacuum_sep: float = 15.0, bond_dist: float = 1.44):
    # get normalized parameters
    n, m, n_t, m_t = _calc_cnt_param(n, m)

    # build a "flat" cnt in the xz plane
    flat_lattice, flat_coords = _gen_flat_cnt(n, m, n_t, m_t)

    # we're mapping the flat cnt to a cylinder so just get the radius
    radius = np.linalg.norm(flat_lattice[0]) / (2 * np.pi)
    # then map them to a cylinder
    cnt_coords = flat_coords
    thetas = cnt_coords[:, 0] / radius

    cnt_coords = np.column_stack((
            radius * np.sin(thetas),
            radius * np.cos(thetas),
            cnt_coords[:, 2],
    ))

    # build the lattice matrix for the new system
    cell_width = 2 * radius + vacuum_sep / bond_dist
    lattice = np.diag([
        cell_width,
        cell_width,
        flat_lattice[2, 2],
    ])

    # center the cnt in the new unit cell
    cnt_coords += (lattice[0] + lattice[1]) / 2

    # finally, build the new structure (note: scale by bond_dist)
    return Structure(
        lattice=Lattice(matrix=lattice * bond_dist),
        species=[Element.C] * len(cnt_coords),
        # note: offset coords to be in the center of the unit cell
        coords=cnt_coords * bond_dist,
        coords_are_cartesian=True,
    )
