#!/usr/bin/env python

from setuptools import setup, find_packages
import os

__author__ = "Colin Daniels"
__maintainer__ = "Colin Daniels"
__email__ = "colin.r.daniels@gmail.com"
__date__ = "Apr 27, 2019"

module_dir = os.path.dirname(os.path.abspath(__file__))

if __name__ == "__main__":
    setup(
        name='pyputil',
        version='0.0.1',
        description='Python physics utility programs',
        long_description='',
        url='https://github.com/colin-daniels/pyputil',
        author='Colin Daniels',
        author_email='colin.r.daniels@gmail.com',
        license='MIT',
        packages=find_packages(),
        package_data={},
        # required due to use of @dataclass
        python_requires='>=3.7',
        zip_safe=False,
        install_requires=[
            # note: CairoSVG actually uses semver
            'CairoSVG~=2.3',
            'h5py~=2.9.0',
            'lxml~=4.3.3',
            'numpy~=1.16.3',
            'phonopy==2.1.3',
            'Pillow~=6.2.0',
            'pymatgen==2019.4.11',
            'ruamel.yaml~=0.15.94',
            'scipy~=1.2.1',
        ],
        classifiers=['Programming Language :: Python',
                     'Development Status :: 3 - Alpha',
                     'Intended Audience :: Science/Research',
                     'License :: OSI Approved :: MIT License',
                     'Operating System :: OS Independent',
                     'Topic :: Other/Nonlisted Topic',
                     'Topic :: Scientific/Engineering'],
        entry_points={
            'console_scripts': [
                'modeplot = pyputil.bin.modeplot:main',
                'structure-gen = pyputil.bin.structure_gen:main',
                'rsp2-util = pyputil.bin.rsp2_util:main',
            ]
        }
    )
