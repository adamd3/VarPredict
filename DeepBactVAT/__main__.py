import os, sys
import argparse
import textwrap

## TODO: import other modules from the DeepBactVAT here

from .__init__ import __version__


class UltimateHelpFormatter(
    argparse.RawTextHelpFormatter, argparse.ArgumentDefaultsHelpFormatter):
    pass



def main():
    parser = argparse.ArgumentParser(formatter_class=UltimateHelpFormatter)
    subparsers = parser.add_subparsers(dest="command", title="Available commands")

    parser.add_argument(
        "--version", action="version", version="%(prog)s " + __version__
    )

    # add subcommands
    dbv_subparser = subparsers.add_parser("run",
        help="run the DeepBactVAT pipeline")
    dbv_subparser = dbv_parser(dbv_subparser)

    ## TODO: add subparsers for different functionalities
    ##   that can be run independently from the main pipeline

    # parse arguments and run function
    args = parser.parse_args()
    args.func(args)

    return


if __name__ == "__main__":
    main()
