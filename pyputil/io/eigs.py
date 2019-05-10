import pathlib

import numpy as np

from pyputil.io import yaml, zopen

# convert eigenvalues to wavenumber (factors from rsp2)
SQRT_EIG_TO_THZ = 15.6333043006705  # sqrt(eV/amu)/Angstrom/(2*pi)/THz
THZ_TO_WAVENUM = 33.3564095198152  # THz / (c / cm)
SQRT_EIG_TO_WAVENUM = SQRT_EIG_TO_THZ * THZ_TO_WAVENUM

SPEED_OF_LIGHT_MS = 299_792_458
THZ_TO_WAVENUMBER = 1e12 / (SPEED_OF_LIGHT_MS * 100)


def solve_dynmat(dynmat):
    evals, evecs = np.linalg.eigh(dynmat.todense())

    # eigenvalues are frequency^2 and in the wrong units, so convert them
    # to cm^-1
    evals = SQRT_EIG_TO_WAVENUM * np.multiply(
        np.sign(evals),
        np.sqrt(np.abs(evals))
    )

    # reshape to 3d vectors
    evals = np.asarray(evals)
    evecs = np.asarray(evecs.T).reshape((len(evals), -1, 3))

    # note: eigenvectors are orthonormal, and to convert to displacements
    # they must be divided by sqrt(mass) for each atom
    return evals, evecs


def from_file(path: str):
    path = pathlib.Path(path)

    if path.suffix == ".hdf5":
        import h5py
        with h5py.File(path) as data:
            frequencies = data['frequency'][0]
            # convert from THz to cm^-1
            frequencies *= THZ_TO_WAVENUMBER

            eigs = data['eigenvector'][0]
            eigs = eigs.reshape((len(frequencies), -1, 3))
    elif path.name.startswith("gamma-dynmat") and path.suffix == ".npz":
        import scipy.sparse
        # looks like rsp2 output, so load the dynamical matrix and solve it
        dynmat = scipy.sparse.load_npz(path)
        frequencies, eigs = solve_dynmat(dynmat)
    elif path.suffix == ".npz":
        # maybe just eigenvalues/eigenvectors in an npz file
        data = np.load(path)
        frequencies = data['frequencies']
        eigs = data['eigenvectors']
    else:
        with zopen(path, 'r') as f:
            data = yaml.load(f)
            band = data['phonon'][0]['band']

            eigs = np.array([mode['eigenvector'] for mode in band])
            # make complex
            eigs = eigs[:, :, :, 0] + 1j * eigs[:, :, :, 1]

            frequencies = np.array([mode['frequency'] for mode in band])
            # convert from THz to cm^-1
            frequencies *= THZ_TO_WAVENUMBER

    n_modes = frequencies.shape[0]
    assert n_modes % 3 == 0

    n_atoms = n_modes // 3
    assert eigs.shape[0] == n_modes
    assert eigs.shape[1] == n_atoms
    assert eigs.shape[2] == 3

    return frequencies, eigs
