import numpy as np
import h5py

from pyputil.io import yaml, zopen

SPEED_OF_LIGHT_MS = 299_792_458
THZ_TO_WAVENUMBER = 1e12 / (SPEED_OF_LIGHT_MS * 100)


def load_eigs_phonopy(filename: str):
    if filename.endswith(".hdf5"):
        with h5py.File(filename) as data:
            frequencies = data['frequency'][0]
            eigs = data['eigenvector'][0]
            eigs = eigs.reshape((len(frequencies), -1, 3))
    else:
        with zopen(filename, 'r') as f:
            data = yaml.load(f)
            band = data['phonon'][0]['band']

            eigs = np.array([
                np.array(mode['eigenvector'])[:, :, 0]
                + 1j * np.array(mode['eigenvector'])[:, :, 1]
                for mode in band
            ])
            frequencies = THZ_TO_WAVENUMBER * np.array([mode['frequency']
                                                        for mode in band])

    return frequencies, eigs
