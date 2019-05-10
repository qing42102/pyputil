import argparse
import sys

from pyputil.structure.gnr import generate_periodic_agnr, generate_finite_agnr
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
        '-p', '--periodic',
        action='store_true',
        help='generate a periodic GNR instead of finite')

    # set function to be run if we select this subcommand
    gnr_parser.set_defaults(func=run_gnr)


def run_gnr(args):
    if not args.periodic and not args.length:
        print("missing --length argument for finite gnr", file=sys.stderr)
        sys.exit(1)
    elif args.periodic:
        if not args.output:
            args.output = f"periodic-agnr-{args.width}.vasp"

        gnr = generate_periodic_agnr(args.width)
        if args.length:
            gnr.make_supercell([args.length, 1, 1])

        gnr.to(filename=args.output, fmt="poscar")
    else:
        if not args.output:
            args.output = f"finite-gnr-{args.width}x{args.length}.vasp"

        gnr = generate_finite_agnr(args.width, args.length)
        gnr.to(filename=args.output, fmt="poscar")


def setup_cnt_opts(subparsers):
    # gnr options
    gnr_parser = subparsers.add_parser(
        'cnt',
        help='generate carbon nanotubes in POSCAR format',
        description='Generate carbon nanotube structures in POSCAR format.')

    gnr_parser.add_argument(
        '-o', '--output',
        type=str,
        required=False,
        help='output filename, default is cnt-<n>x<m>.vasp')

    gnr_parser.add_argument(
        '-n',
        metavar='N',
        type=int,
        help='chirality',
        required=True)

    gnr_parser.add_argument(
        '-m',
        metavar='M',
        type=int,
        help='chirality',
        required=True)

    gnr_parser.set_defaults(func=run_cnt)


def run_cnt(args):
    if not args.output:
        args.output = f"cnt-{args.n}x{args.m}.vasp"

    cnt = generate_cnt(args.n, args.m)
    cnt.to(filename=args.output, fmt="poscar")


if __name__ == '__main__':
    main()
