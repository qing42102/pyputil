from pyputil.modeplot import RenderSettings, ModeRenderer
import argparse


def main():
    parser = argparse.ArgumentParser(description='Generate SVG phonon mode plots from phonopy output.')

    parser.add_argument(
        '-c', '--config',
        type=str,
        default=None,
        help='optional render settings config in yaml format')

    parser.add_argument(
        '-s', '--supercell',
        nargs=3,
        type=int,
        default=[1, 1, 1])

    parser.add_argument(
        '-i', '--input',
        type=str,
        required=True,
        help='input structure file (e.g. POSCAR)')

    parser.add_argument(
        '-e', '--eigs',
        type=str,
        required=True,
        help='eigenvalue input file (e.g. band.yaml, eigs.yaml)')

    args = parser.parse_args()

    renderer = ModeRenderer(
        structure_filename=args.input,
        eigs_filename=args.eigs,
        supercell=args.supercell,
        settings=RenderSettings.from_file_or_default(args.config),
    )

    for i in range(len(renderer.frequencies)):
        filename = 'mode_{}.svg'.format(i + 1)

        print(filename)
        renderer.render(i).write(
            filename,
            encoding='utf-8',
            pretty_print=True,
            xml_declaration=True
        )


if __name__ == '__main__':
    main()
