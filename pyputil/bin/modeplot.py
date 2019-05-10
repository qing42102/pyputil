from multiprocessing.pool import Pool
import sys

from pymatgen import Structure
import numpy as np
import typing as tp

from pyputil.io.phonopy import load_eigs_phonopy
from pyputil.modeplot import RenderSettings, ModeRenderer
import argparse

from pyputil.svg import write_svg_gif


def render_svg(mode_ids: tp.Union[int, tp.List[int]]):
    if type(mode_ids) == list:
        for m in mode_ids:
            render_svg(m)
    else:
        filename = f'mode_{mode_ids + 1}.svg'
        mode_svg = renderer.render(mode_ids)
        mode_svg.write(filename)
        print(filename)


def render_gif(mode_ids: tp.Union[int, tp.List[int]]):
    if type(mode_ids) == list:
        for m in mode_ids:
            render_gif(m)
    else:
        filename = f'mode_{mode_ids + 1}.gif'
        write_svg_gif(
            filename,
            renderer.render_mode_anim(mode_ids),
            duration=renderer.settings.gif['frame-delay'],
        )
        print(filename)


class PrintDefaultsAction(argparse.Action):
    def __call__(self, parser, namespace, values, option_string=None):
        defaults = RenderSettings()
        defaults.to_yaml(sys.stdout)
        parser.exit()


def main():
    global renderer

    parser = argparse.ArgumentParser(
        prog="modeplot",
        description='Generate SVG/GIF phonon mode plots from phonopy output.')

    parser.register('action', 'print_defaults', PrintDefaultsAction)

    parser.add_argument(
        '-c', '--config',
        type=str,
        default=None,
        help='optional render settings config in yaml format')

    parser.add_argument(
        '-s', '--supercell',
        metavar='N',
        nargs=3,
        type=int,
        default=[1, 1, 1])

    parser.add_argument(
        '-i', '--input',
        metavar='STRUCTURE_FILE',
        type=str,
        required=True,
        help='input structure file (e.g. POSCAR, input.cif)')

    parser.add_argument(
        '-e', '--eigs',
        metavar='EIGENVECTOR_FILE',
        type=str,
        required=True,
        help='eigenvalue input file (e.g. band.yaml, qpoints.hdf5)')

    parser.add_argument(
        '-g', '--gif',
        metavar='MODE_ID',
        type=int,
        help='render a single mode gif')

    parser.add_argument(
        '--all-gifs',
        action='store_true',
        help='render all modes as gifs in addition to svgs')

    parser.add_argument(
        '--parallel',
        action='store_true',
        help='try to render modes in parallel')

    parser.add_argument(
        '--print-defaults',
        nargs=0,
        action='print_defaults',
        help='print default yaml settings to stdout and exit')

    args = parser.parse_args()

    # settings
    settings = RenderSettings.from_file_or_default(args.config)

    # read structure
    structure = Structure.from_file(args.input)

    # initial translation
    structure.translate_sites(
        np.arange(structure.num_sites),
        np.array(settings.translation, dtype=float),
        to_unit_cell=True)

    # read eigenvectors
    frequencies, eigs = load_eigs_phonopy(args.eigs)

    # make supercell if specified
    if args.supercell is not None:
        supercell_images = np.prod(args.supercell)
        structure.make_supercell(args.supercell)
        eigs = np.repeat(eigs, supercell_images, axis=1)

    renderer = ModeRenderer(
        structure=structure,
        frequencies=frequencies,
        eigenvectors=eigs,
        settings=settings,
    )

    if args.gif:
        assert 0 < args.gif <= len(renderer.frequencies)
        render_gif(args.gif - 1)
    else:
        mode_ids = list(range(len(frequencies)))
        if args.parallel:
            with Pool() as pool:
                pool.map(render_svg, mode_ids)

                if args.all_gifs:
                    pool.map(render_gif, mode_ids)
        else:
            render_svg(mode_ids)

            if args.all_gifs:
                render_gif(mode_ids)


if __name__ == '__main__':
    main()
