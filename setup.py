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
        version="0.0.1",
        description='Python physics utility programs',
        # long_description=open(os.path.join(module_dir, 'README.md')).read(),
        long_description='',
        url='https://github.com/colin-daniels/pyputil',
        author='Colin Daniels',
        author_email='colin.r.daniels@gmail.com',
        license='MIT',
        packages=find_packages(),
        package_data={},
        zip_safe=False,
        # install_requires=['ruamel.yaml>=0.15.35'],
        classifiers=['Programming Language :: Python',
                     'Development Status :: 3 - Alpha',
                     'Intended Audience :: Science/Research',
                     'Intended Audience :: System Administrators',
                     'Intended Audience :: Information Technology',
                     'Operating System :: OS Independent',
                     'Topic :: Other/Nonlisted Topic',
                     'Topic :: Scientific/Engineering'],
        entry_points={
            'console_scripts': [
                'modeplot = pyputil.bin.modeplot:main',
                'structure-gen = pyputil.bin.structure_gen:main',
            ]
        }
    )
