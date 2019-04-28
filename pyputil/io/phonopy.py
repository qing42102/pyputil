import numpy as np

from pyputil.io import yaml, zopen

SPEED_OF_LIGHT_MS = 299_792_458
THZ_TO_WAVENUMBER = 1e12 / (SPEED_OF_LIGHT_MS * 100)


def load_eigs_band_yaml(filename: str):
    with zopen(filename, 'r') as f:
        data = yaml.load(f)
        band = data['phonon'][0]['band']

        # real only
        eigs = [np.array(mode['eigenvector'])[:, :, 0] for mode in band]
        frequencies = THZ_TO_WAVENUMBER * np.array([mode['frequency'] for mode in band])

    return frequencies, eigs
