import argparse

import numpy as np


def main():
    # create the top-level parser
    parser = argparse.ArgumentParser(
        prog='rsp2-util',
        description="Utilities for working with rsp2 and rsp2 output")

    # subcommands
    subparsers = parser.add_subparsers(title="subcommands")

    setup_solve_opts(subparsers)

    # error if no subcommand
    parser.set_defaults(func=lambda _: parser.error("missing subcommand"))

    # parse and call selected subcommand
    args = parser.parse_args()
    args.func(args)


def setup_solve_opts(subparsers):
    # gnr options
    solve_parser = subparsers.add_parser(
        'solve-dynmat',
        help='solve rsp2 dynamical matrix',
        description='Solve rsp2 dynamical matrix and output eigenvectors/eigenvalues using numpy.linalg.eigh')

    solve_parser.add_argument('DYNMAT', help='dynmat file (sparse npz format, e.g. gamma-dynmat-01.npz)')
    solve_parser.add_argument('--output', '-o', type=str, required=True,
                              help="frequency/eigenvector output filename (npz format, uses numpy.savez_compressed)")

    # set function to be run if we select this subcommand
    solve_parser.set_defaults(func=run_solve)


def run_solve(args):
    from pyputil.io.eigs import solve_dynmat
    import scipy.sparse

    dynmat = scipy.sparse.load_npz(args.DYNMAT)
    evals, evecs = solve_dynmat(dynmat)

    np.savez_compressed(
        args.output,
        frequencies=evals,
        eigenvectors=evecs)


if __name__ == '__main__':
    main()
