import math
import typing as tp
from dataclasses import dataclass, asdict

import numpy as np
from pymatgen import Structure
import lxml.etree as etree

from pyputil.io import yaml
from pyputil.io.phonopy import load_eigs_band_yaml
from pyputil.misc import default_field
from pyputil.structure.bonds import calculate_bonds, bonds_to_positions
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

@dataclass
class ViewBounds:
    x: float
    y: float
    w: float
    h: float


class ModeRenderer:
    def __init__(
            self,
            structure_filename: str,
            eigs_filename: str,
            # defaults to [1, 1, 1]
            supercell: tp.Optional[tp.List[int]] = None,
            settings: RenderSettings = RenderSettings()
    ):
        self.settings = settings

        # read structure
        structure = Structure.from_file(structure_filename)
        # initial translation
        structure.translate_sites(
            np.arange(structure.num_sites),
            np.array(settings.translation, dtype=float),
            to_unit_cell=True)

        # read eigenvectors
        self.frequencies, self.eigenvectors = load_eigs_band_yaml(eigs_filename)

        # make supercell if specified
        if supercell is not None:
            supercell_images = np.prod(supercell)
            structure.make_supercell(supercell)
            # Poscar(structure).write_file("supercell.vasp")
            self.eigenvectors = np.repeat(
                self.eigenvectors, supercell_images, axis=1)
        else:
            supercell_images = 1

        self.structure = structure
        self.bonds = calculate_bonds(structure)
        self.transform = settings.scaling * \
                         np.array(settings.rotation).astype(float)

        # transformed positions
        self.coords = self.transform_coords(
            self.structure.cart_coords)
        self.lattice = self.transform_coords(
            self.structure.lattice.matrix)
        self.bond_coords = self._calc_bond_coords()

        # transformed and normalized displacements
        self.normalized_displacements = self._calc_displacements()

        # view boundary
        self.view_bounds = self._calc_view_bounds()

        self.frequencies, self.eigenvectors = load_eigs_band_yaml(eigs_filename)
        self.eigenvectors = np.array([
            np.repeat(e, supercell_images, axis=0) for e in self.eigenvectors
        ])

    def _calc_displacements(self):
        # scale by 1 / sqrt(mass) to get displacements
        sqrt_masses = np.array([
            math.sqrt(mass_from_symbol(s.name))
            for s in self.structure.species
        ])

        # only look at real components
        real_eigs = np.real_if_close(self.eigenvectors)
        disps = real_eigs / sqrt_masses[np.newaxis, :, np.newaxis]

        # for display purposes, normalize displacements such that the maximum
        # displacement for each mode is of length 1
        displacement_norms = np.linalg.norm(disps, axis=2, keepdims=True)
        disps /= np.max(displacement_norms, axis=1, keepdims=True)

        return self.transform_coords(disps)

    def _calc_bond_coords(self):
        bond_positions = bonds_to_positions(
            self.bonds, self.lattice, self.coords)

        if not self.settings.bonds['draw-periodic']:
            bond_is_periodic = np.any(self.bonds[:, :3], axis=1)
            return bond_positions[np.logical_not(bond_is_periodic)]
        else:
            return bond_positions

    def _calc_view_bounds(self):
        rset = self.settings
        padding = max(rset.padding, rset.displacements['max-length'])
        padding *= rset.scaling

        min_pos = self.coords.min(axis=0) - padding
        max_pos = self.coords.max(axis=0) + padding
        range_pos = max_pos - min_pos

        return ViewBounds(
            x=min_pos[0],
            y=min_pos[1],
            w=range_pos[0],
            h=range_pos[1]
        )

    def transform_coords(self, coords):
        return np.dot(coords, self.transform.T)

    def render(self, mode_id: tp.Optional[int]) -> etree.ElementTree:
        rset = self.settings

        # find structure bounds
        bounds = self.view_bounds
        svg = self.svg_base(bounds)

        # definitions for use later
        defs = etree.SubElement(svg, 'defs')

        if rset.background_color is not None:
            self.svg_background(bounds, svg)

        if rset.draw_info_text and mode_id is not None:
            self.svg_info_text(bounds, mode_id, svg)

        self.svg_bonds(svg, self.bond_coords)

        if mode_id is not None:
            self.svg_displacements(svg, defs, mode_id)

        self.svg_atoms(svg, self.coords)

        return etree.ElementTree(svg)

    def svg_base(self, bounds: ViewBounds):
        return etree.Element("svg", nsmap=NSMAP, attrib={
            'width': '{:.1f}'.format(bounds.w),
            'height': '{:.1f}'.format(bounds.h),
            'version': '1.1',
            'viewBox': '{:.1f} {:.1f} {:.1f} {:.1f}'.format(
                bounds.x, bounds.y,
                bounds.w, bounds.h,
            ),
            'xmlns': "http://www.w3.org/2000/svg",
        })

    def svg_background(self, bounds: ViewBounds, svg):
        return etree.SubElement(svg, 'rect', attrib={
            'id': 'background',
            'x': '{:.1f}'.format(bounds.x),
            'y': '{:.1f}'.format(bounds.y),
            'width': '{:.1f}'.format(bounds.w),
            'height': '{:.1f}'.format(bounds.h),
            'fill': self.settings.background_color,
        })

    def svg_info_text(self, bounds: ViewBounds, mode_id, svg):
        text_size = self.settings.info_text_size * self.settings.scaling
        mode_info = etree.SubElement(svg, 'text', attrib={
            'id': 'mode-info',
            'x': '{:.1f}'.format(bounds.x + text_size / 3),
            'y': '{:.1f}'.format(bounds.y + bounds.h - text_size / 3),
            'font-family': 'monospace',
            'font-size': '{}'.format(text_size),
            'font-weight': 'bold',
        })
        mode_info.text = "id: {} | frequency: {:.4f} cm^-1"\
            .format(mode_id + 1, self.frequencies[mode_id])

    def svg_atoms(self, svg, coords):
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

        for specie, [x, y, _] in zip(self.structure.species, coords):
            info = rset.atom_types[str(specie)]
            etree.SubElement(info['group'], 'circle', attrib={
                'r': str(info['radius'] * rset.scaling),
                'cx': '{:.1f}'.format(x),
                'cy': '{:.1f}'.format(y),
            })

    def svg_displacements(
            self,
            svg: etree.Element,
            defs: etree.Element,
            mode_id: int,
    ):
        rset = self.settings

        # displacements (phonon mode)
        etree.SubElement(defs, 'path', attrib={
            'id': 'arrowhead',
            'd': 'M -1 0 h 2 l -1 1.732 z',
            'fill': rset.displacements['color'],
            'stroke-width': '0',
            'transform': 'scale({})'.format(
                rset.scaling * rset.displacements['arrow-width'] / 2)
        })

        displacement_group = etree.SubElement(svg, 'g', nsmap=NSMAP, attrib={
            'id': 'displacements',
            'stroke': '#EB1923',
            'stroke-width':
                str(rset.displacements['stroke-width'] * rset.scaling),
        })

        disps = self.normalized_displacements[mode_id]
        disps *= rset.displacements['max-length']
        angles = np.degrees(np.arctan2(disps[:, 1], disps[:, 0]) - np.pi / 2)

        for pos, disp, angle in zip(self.coords, disps, angles):
            end = pos + disp
            etree.SubElement(displacement_group, 'path', attrib={
                'd': 'M {:.1f} {:.1f} {:.1f} {:.1f}'.format(
                    pos[0], pos[1], end[0], end[1]
                )
            })

            # add arrow head via href
            etree.SubElement(displacement_group, 'use', nsmap=NSMAP, attrib={
                '{http://www.w3.org/1999/xlink}href': '#arrowhead',
                'x': '{:.1f}'.format(end[0]),
                'y': '{:.1f}'.format(end[1]),
                'transform': 'rotate({:.1f} {:.1f} {:.1f})'.format(
                    angle, end[0], end[1]
                ),
            })

    def svg_bonds(self, svg, bond_coords):
        rset = self.settings
        bond_group = etree.SubElement(svg, 'g', attrib={
            'id': 'bonds',
            'stroke': rset.bonds['stroke'],
            'stroke-width': str(rset.bonds['stroke-width'] * rset.scaling),
        })

        for atom_a, atom_b in bond_coords:
            etree.SubElement(bond_group, 'path', attrib={
                'd': 'M {:.1f} {:.1f} {:.1f} {:.1f}'.format(
                    atom_a[0], atom_a[1], atom_b[0], atom_b[1]
                )
            })
