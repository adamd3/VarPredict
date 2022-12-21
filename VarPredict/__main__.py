import os, sys
import argparse
import textwrap
from .BactVarMiner import dbv_parser
from .BactVarMiner import dbv_parser

from .__init__ import __version__


class UltimateHelpFormatter(
    argparse.RawTextHelpFormatter, argparse.ArgumentDefaultsHelpFormatter):
    pass


def get_options(args):


def main():
    parser = argparse.ArgumentParser(formatter_class=UltimateHelpFormatter)
    subparsers = parser.add_subparsers(
        dest="command", title="Available commands")

    parser.add_argument(
        "--version", action="version", version="%(prog)s " + __version__
    )

    # add subcommands
    ncv_subparser = subparsers.add_parser("nested-cv",
        help="estimate model accuracy via nested cross-validation")
    ncv_subparser = ncv_parser(ncv_subparser)

    dbv_subparser = subparsers.add_parser("run",
        help="s")
    dbv_subparser = dbv_parser(dbv_subparser)

    # parse arguments and run function
    args = parser.parse_args()
    args.func(args)

    return


if __name__ == "__main__":
    main()
