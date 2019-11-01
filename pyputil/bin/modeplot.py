import os
import sys
from pathlib import PurePath

from pymatgen import Structure
import numpy as np
import typing as tp

import pyputil.io.eigs
from pyputil.modeplot import RenderSettings, ModeRenderer
import argparse

from pyputil.svg import write_svg_gif


def render_svg(
        info: tp.Union[tp.Tuple[int, str], tp.List[tp.Tuple[int, str]]]
):
    if type(info) == list:
        for m in info:
            render_svg(m)
    else:
        mode_id, filename = info
        filename = filename + ".svg"

        mode_svg = renderer.render(mode_id)
        mode_svg.write(filename)
        print(filename)


def render_gif(
        info: tp.Union[tp.Tuple[int, str], tp.List[tp.Tuple[int, str]]]
):
    if type(info) == list:
        for m in info:
            render_gif(m)
    else:
        mode_id, filename = info
        filename = filename + ".gif"

        write_svg_gif(
            filename,
            renderer.render_mode_anim(mode_id),
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
        help='eigenvector input file (e.g. band.yaml, qpoints.hdf5, gamma-dynmat*.npz)')

    parser.add_argument(
        '-o', '--output',
        metavar='DIR',
        type=str,
        required=False,
        default=None,
        help='output directory, default is the current directory')

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
        '-j', '--parallel',
        const=os.cpu_count() or 1,
        default=None,
        action='store',
        nargs='?',
        type=int,
        metavar='N',
        help='enable parallel processing on N processes, default is number of cores if not specified')

    parser.add_argument(
        '--print-defaults',
        nargs=0,
        action='print_defaults',
        help='print default yaml settings to stdout and exit')

    args = parser.parse_args()
    run(vars(args))


def run(args):
    global renderer
    # settings
    settings = RenderSettings.from_file_or_default(args["config"])

    # read structure
    structure = Structure.from_file(args["input"])

    # initial translation
    structure.translate_sites(
        np.arange(structure.num_sites),
        np.array(settings.translation, dtype=float),
        to_unit_cell=True)

    # read eigenvectors
    frequencies, eigs = pyputil.io.eigs.from_file(args["eigs"])

    # make supercell if specified
    if args["supercell"] is not None:
        supercell_images = np.prod(args["supercell"])
        structure.make_supercell(args["supercell"])
        eigs = np.repeat(eigs, supercell_images, axis=1)

    renderer = ModeRenderer(
        structure=structure,
        frequencies=frequencies,
        eigenvectors=eigs,
        settings=settings,
    )

    output = PurePath(os.getcwd())
    if args["output"]:
        output = PurePath(args["output"])

    if args["gif"]:
        mode_id = args["gif"]
        assert 0 < mode_id <= len(renderer.frequencies)
        render_gif((
            mode_id - 1,
            str(output / f"mode_{mode_id}")
        ))
    else:
        modes = [
            (idx, str(output / f"mode_{idx + 1}"))
            for idx in range(len(frequencies))
        ]

        if args["parallel"] and args["parallel"] > 1:
            from multiprocessing.pool import Pool

            with Pool(processes=args["parallel"]) as pool:
                pool.map(render_svg, modes)

                if args["all_gifs"]:
                    pool.map(render_gif, modes)
        else:
            render_svg(modes)

            if args["all_gifs"]:
                render_gif(modes)


if __name__ == '__main__':
    main()
