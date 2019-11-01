import argparse
import sys

from pyputil.structure import add_vacancy_defects
from pyputil.structure.gnr import generate_periodic_agnr, generate_finite_agnr, \
    generate_periodic_zgnr
from pyputil.structure.nanotube import generate_cnt


def main():
    # create the top-level parser
    parser = argparse.ArgumentParser(
        prog='structure-gen',
        description="Generate various structures, see subcommands for more information")

    # subcommands
    subparsers = parser.add_subparsers(title="subcommands")
    setup_gnr_opts(subparsers)
    setup_cnt_opts(subparsers)
    setup_flake_opts(subparsers)
    setup_add_defects_opts(subparsers)

    # error if no subcommand
    parser.set_defaults(func=lambda _: parser.error("missing subcommand"))

    # parse and call selected subcommand
    args = parser.parse_args()
    args.func(args)


def setup_gnr_opts(subparsers):
    # gnr options
    gnr_parser = subparsers.add_parser(
        'gnr',
        help='generate graphene nanoribbons in POSCAR format',
        description='Generate finite GNR structures in POSCAR format.')

    gnr_parser.add_argument(
        '-o', '--output',
        type=str,
        required=False,
        help='output filename, default is finite-gnr-<width>x<length>.vasp')

    gnr_parser.add_argument(
        '-w', '--width',
        metavar='N',
        type=int,
        required=True)

    gnr_parser.add_argument(
        '-l', '--length',
        metavar='M',
        type=int,
        default=None)

    gnr_parser.add_argument(
        '--no-hydrogen',
        action='store_true',
        help='don\'t hydrogenate edges')

    gnr_parser.add_argument(
        '-p', '--periodic',
        action='store_true',
        help='generate a periodic GNR instead of finite')

    gnr_parser.add_argument(
        '-z', '--zigzag',
        action='store_true',
        help='generate a ZGNR instead of AGNR, currently limited to periodic')

    # set function to be run if we select this subcommand
    gnr_parser.set_defaults(func=lambda args: run_gnr(vars(args)))


def run_gnr(args):
    if not args["periodic"] and not args["length"]:
        print("missing --length argument for finite gnr", file=sys.stderr)
        sys.exit(1)
    elif args["periodic"]:
        if args["zigzag"]:
            name = "z"
            func = generate_periodic_zgnr
        else:
            name = "a"
            func = generate_periodic_agnr

        if not args["output"]:
            args["output"] = f"periodic-{name}gnr-{args['width']}.vasp"

        gnr = func(args["width"], hydrogen=not args["no_hydrogen"])
        if args["length"]:
            gnr.make_supercell([args["length"], 1, 1])

        gnr.to(filename=args["output"], fmt="poscar")
    else:
        if args["zigzag"]:
            print("only periodic zigzag GNR generation is supported",
                  file=sys.stderr)
            sys.exit(1)

        if not args["output"]:
            args["output"] = f"finite-gnr-{args['width']}x{args['length']}.vasp"

        gnr = generate_finite_agnr(args["width"], args["length"], hydrogen=not args["no_hydrogen"])
        gnr.to(filename=args["output"], fmt="poscar")


def setup_flake_opts(subparsers):
    parser = subparsers.add_parser(
        'flake',
        help='generate circular graphene flakes in POSCAR format',
        description='Generate graphene flake structures in POSCAR format.')

    parser.add_argument(
        '-o', '--output',
        type=str,
        required=True,
        help='output prefix')

    parser.add_argument(
        '-r', '--radius',
        metavar='RADIUS',
        type=float,
        required=True)

    parser.set_defaults(func=lambda args: run_flake(vars(args)))


def run_flake(args):
    from pyputil.structure.flake import generate_graphene_flakes
    [a, b, c] = generate_graphene_flakes(args["radius"])
    a.to(filename=f"{args['output']}-a.vasp", fmt="poscar")
    b.to(filename=f"{args['output']}-b.vasp", fmt="poscar")
    c.to(filename=f"{args['output']}-c.vasp", fmt="poscar")
        

def setup_cnt_opts(subparsers):
    # gnr options
    parser = subparsers.add_parser(
        'cnt',
        help='generate carbon nanotubes in POSCAR format',
        description='Generate carbon nanotube structures in POSCAR format.')

    parser.add_argument(
        '-o', '--output',
        type=str,
        required=False,
        help='output filename, default is cnt-<n>x<m>.vasp')

    parser.add_argument(
        '-n',
        metavar='N',
        type=int,
        help='chirality',
        required=True)

    parser.add_argument(
        '-m',
        metavar='M',
        type=int,
        help='chirality',
        required=True)

    parser.set_defaults(func=lambda args: run_cnt(vars(args)))


def run_cnt(args):
    n = args["n"]
    m = args["m"]
    output = args["output"]

    if not output:
        output = f"cnt-{n}x{m}.vasp"

    cnt = generate_cnt(n, m)
    cnt.to(filename=output, fmt="poscar")


def setup_add_defects_opts(subparsers):
    # gnr options
    parser = subparsers.add_parser(
        'add_defects',
        help='add vacancy defects to structures in POSCAR format',
        description='Add random vacancy defects to carbon structures in POSCAR format.')

    parser.add_argument(
        '-o', '--output-prefix',
        type=str,
        default="defective",
        metavar='PREFIX',
        help='output *prefix*, default is "defective", full paths are "<prefix>-N.vasp".')

    parser.add_argument(
        '-p', "--defect-percent",
        metavar='P',
        type=float,
        help='percent chance for a carbon atom to be removed',
        required=True)

    parser.add_argument(
        '-n', "--num-structures",
        metavar='N',
        type=int,
        help='number of randomly-defected structures to generate from the input',
        required=True)

    parser.add_argument(
        '-i', '--input',
        type=str,
        required=True,
        help='input filename for the structure that will have defects added')

    parser.set_defaults(func=lambda args: run_add_defects(vars(args)))


def run_add_defects(args):
    from pymatgen import Structure

    # read original structure
    pristine = Structure.from_file(args["input"])

    defect_probability = args["defect_percent"] / 100.0
    prefix = args["output_prefix"]

    for n in range(args["num_structures"]):
        output = f"{prefix}-{n}.vasp"
        defective = add_vacancy_defects(pristine, defect_probability)
        defective.to(filename=output, fmt="poscar")


if __name__ == '__main__':
    main()
