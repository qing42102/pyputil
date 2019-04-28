import math
import typing as tp
from dataclasses import dataclass, asdict

import numpy as np
from numba import jit
from pymatgen import Structure
from pymatgen.io.vasp import Poscar
import lxml.etree as etree

from pyputil.io import yaml
from pyputil.io.phonopy import load_eigs_band_yaml
from pyputil.misc import default_field
from pyputil.structure.bonds import calculate_bonds, bond_to_positions, bonds_to_positions
from pyputil.structure.element import mass_from_symbol

NSMAP = {'xlink': 'http://www.w3.org/1999/xlink'}


@dataclass(unsafe_hash=True)
class RenderSettings:
    # image size multiplier
    scaling: float = 30.0
    # how much to translate the structure
    translation: tp.List[float] = default_field([0.0, 0.0, 0.0])
    # structure rotation (applied after translation)
    rotation: tp.List[tp.List[float]] = default_field([
        [1, 0, 0],
        [0, 1, 0],
        [0, 0, 1],
    ])

    # extra space to add to fit structure in the image
    padding: float = 2.0

    # text settings for mode info
    draw_info_text: bool = True
    info_text_size: float = padding / 4

    # can be none for transparent
    background_color: tp.Optional[str] = '#FFFFFF'

    bonds: tp.Dict = default_field({
        'stroke': '#7F8285',
        'stroke-width': 0.1251,
        # whether to draw bonds across periodic bounds
        'draw-periodic': False,
    })

    displacements: tp.Dict = default_field({
        'color': '#EB1923',
        'stroke-width': 0.1251,
        'max-length': 1.44 / 2,
        'arrow-width': 3 * 0.1251,
    })

    # per-atom render settings
    atom_types: tp.Dict = default_field({
        # Carbon
        "C": {
            'radius': 0.2607,
            'fill':   '#111417',
            'stroke': '#111417',
            'stroke-width': 0.0834,
        },
        # Hydrogen
        "H": {
            'radius': 0.1669,
            'fill':   '#FFFFFF',
            'stroke': '#111417',
            'stroke-width': 0.0834,
        },
    })

    def to_yaml(self, stream):
        return yaml.dump(asdict(self), stream)

    def to_file(self, filename: str):
        with open(filename, 'w') as f:
            self.to_yaml(f)

    @classmethod
    def from_file(cls, filename: str):
        with open(filename, 'r') as f:
            return cls(**yaml.load(f))

    @classmethod
    def from_file_or_default(cls, filename: tp.Optional[str]):
        if filename is None:
            return cls()

        try:
            return cls.from_file(filename)
        except FileNotFoundError:
            print("Warning, could not open settings file.")
            return cls()


class ModeRenderer:
    def __init__(
            self,
            structure_filename: str,
            eigs_filename: str,
            # defaults to [1, 1, 1]
            supercell: tp.Optional[tp.List[int]] = None,
            settings: RenderSettings = RenderSettings()
    ):
        # read structure
        structure = Structure.from_file(structure_filename)
        # initial translation
        structure.translate_sites(
            np.arange(structure.num_sites), np.array(settings.translation, dtype=float), to_unit_cell=True)

        # make supercell if specified
        if supercell is not None:
            supercell_images = np.prod(supercell)
            structure.make_supercell(supercell)
            Poscar(structure).write_file("supercell.vasp")
        else:
            supercell_images = 1

        self.structure = structure
        self.bonds = calculate_bonds(structure)

        bond_is_periodic = np.any(self.bonds[:, :3], axis=1)
        self.bond_positions = bonds_to_positions(self.bonds, self.structure)
        self.non_periodic_bond_positions = self.bond_positions[np.logical_not(bond_is_periodic)]

        self.frequencies, self.eigenvectors = load_eigs_band_yaml(eigs_filename)
        self.eigenvectors = np.array([
            np.repeat(e, supercell_images, axis=0) for e in self.eigenvectors
        ])

        self.settings = settings

    def render(self, mode_id: tp.Optional[int]):
        rset = self.settings

        # our coordinate transform then consists of a rotation then scaling
        transform = rset.scaling * np.array(rset.rotation).astype(float)

        def coord_transform(coords):
            return np.dot(transform, coords.T).T

        # find structure bounds
        padding = rset.scaling * rset.padding
        transformed_coords = coord_transform(self.structure.cart_coords)
        min_pos = transformed_coords.min(axis=0) - padding
        max_pos = transformed_coords.max(axis=0) + padding
        range_pos = max_pos - min_pos

        svg = etree.Element("svg", nsmap=NSMAP, attrib={
            'width': '{:.1f}'.format(range_pos[0]),
            'height': '{:.1f}'.format(range_pos[1]),
            'version': '1.1',
            'viewBox': '{:.1f} {:.1f} {:.1f} {:.1f}'.format(
                min_pos[0], min_pos[1],
                range_pos[0], range_pos[1],
            ),
            'xmlns': "http://www.w3.org/2000/svg",
        })

        # definitions for use later
        defs = etree.SubElement(svg, 'defs')

        if rset.background_color is not None:
            etree.SubElement(svg, 'rect', attrib={
                'id': 'background',
                'x': '{:.1f}'.format(min_pos[0]),
                'y': '{:.1f}'.format(min_pos[1]),
                'width': '{:.1f}'.format(range_pos[0]),
                'height': '{:.1f}'.format(range_pos[1]),
                'fill': rset.background_color,
            })

        if rset.draw_info_text and mode_id is not None:
            text_size = rset.info_text_size * rset.scaling
            mode_info = etree.SubElement(svg, 'text', attrib={
                'id': 'mode-info',
                'x': '{:.1f}'.format(min_pos[0] + text_size / 3),
                'y': '{:.1f}'.format(max_pos[1] - text_size / 3),
                'font-family': 'monospace',
                'font-size': '{}'.format(text_size),
                'font-weight': 'bold',
            })
            mode_info.text = "id: {} | frequency: {:.4f} cm^-1".format(mode_id + 1, self.frequencies[mode_id])

        self.svg_bonds(svg, coord_transform)

        if mode_id is not None:
            self.svg_displacements(svg, defs, mode_id, coord_transform)

        self.svg_atoms(svg, coord_transform)

        return etree.ElementTree(svg)

    def svg_atoms(self, svg, coord_transform):
        rset = self.settings

        atom_group = etree.SubElement(svg, 'g', attrib={
            'id': 'atoms'
        })

        for key, value in rset.atom_types.items():
            value['group'] = etree.SubElement(atom_group, 'g', attrib={
                'id': '{}-atoms'.format(key.lower()),
                'fill': value['fill'],
                'stroke': value['stroke'],
                'stroke-width': str(value['stroke-width'] * rset.scaling),
            })
        for site in self.structure.sites:
            x, y, z = coord_transform(site.coords)
            info = rset.atom_types[str(site.specie)]

            etree.SubElement(info['group'], 'circle', attrib={
                'r': str(info['radius'] * rset.scaling),
                'cx': '{:.1f}'.format(x),
                'cy': '{:.1f}'.format(y),
            })

    def svg_displacements(self, svg: etree.Element, defs: etree.Element, mode_id: int, coord_transform):
        rset = self.settings

        # displacements (phonon mode)
        etree.SubElement(defs, 'path', attrib={
            'id': 'arrowhead',
            'd': 'M -1 0 h 2 l -1 1.732 z',
            'fill': rset.displacements['color'],
            'stroke-width': '0',
            'transform': 'scale({})'.format(rset.scaling * rset.displacements['arrow-width'] / 2)
        })

        displacement_group = etree.SubElement(svg, 'g', nsmap=NSMAP, attrib={
            'id': 'displacements',
            'stroke': '#EB1923',
            'stroke-width': str(rset.displacements['stroke-width'] * rset.scaling),
        })

        eigs = self.eigenvectors[mode_id]
        # scale by 1 / sqrt(mass) to get displacements
        masses = np.array([mass_from_symbol(s.name) for s in self.structure.species])
        disps = (eigs.T / np.sqrt(masses)).T

        # normalize, then scale to arrow length
        disps = rset.displacements['max-length'] * disps / np.max(np.linalg.norm(disps, axis=1))

        for pos, disp in zip(coord_transform(self.structure.cart_coords), coord_transform(disps)):
            # pos = coord_transform(pos)
            # disp = coord_transform(disp)
            end = pos + disp

            etree.SubElement(displacement_group, 'path', attrib={
                'd': 'M {:.1f} {:.1f} {:.1f} {:.1f}'.format(
                    pos[0], pos[1], end[0], end[1]
                )
            })

            # add arrow head
            etree.SubElement(displacement_group, 'use', nsmap=NSMAP, attrib={
                '{http://www.w3.org/1999/xlink}href': '#arrowhead',
                'x': '{:.1f}'.format(end[0]),
                'y': '{:.1f}'.format(end[1]),
                'transform': 'rotate({:.1f} {:.1f} {:.1f})'.format(
                    math.degrees(math.atan2(disp[1], disp[0]) - np.pi / 2), end[0], end[1]
                ),
            })

    def svg_bonds(self, svg, coord_transform):
        rset = self.settings

        # bonds
        bond_group = etree.SubElement(svg, 'g', attrib={
            'id': 'bonds',
            'stroke': rset.bonds['stroke'],
            'stroke-width': str(rset.bonds['stroke-width'] * rset.scaling),
        })

        bond_pos = self.bond_positions
        if not rset.bonds['draw-periodic']:
            bond_pos = self.non_periodic_bond_positions

        for pos_a, pos_b in zip(coord_transform(bond_pos[:, 0, :]), coord_transform(bond_pos[:, 1, :])):
            etree.SubElement(bond_group, 'path', attrib={
                'd': 'M {:.1f} {:.1f} {:.1f} {:.1f}'.format(
                    pos_a[0], pos_a[1], pos_b[0], pos_b[1]
                )
            })
