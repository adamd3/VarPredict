import os, sys
import argparse

from .__init__ import __version__

from .rf_model import rf_parser
# from .en_model import en_parser

class UltimateHelpFormatter(
    argparse.RawTextHelpFormatter, argparse.ArgumentDefaultsHelpFormatter):
    pass


def main():
    main_parser = argparse.ArgumentParser(formatter_class=UltimateHelpFormatter)

    # main arguments
    main_parser.add_argument("--filter",
        dest = "filter",
        help = "Filter gene expression data to remove genes with low variance",
        action = "store_true",
        default = False
    )
    main_parser.add_argument(
        "--version", 
        action="version", 
        version="%(prog)s " + __version__
    )

    # argument groups
    io_opts = main_parser.add_argument_group("Input/output files")
    other_opts = main_parser.add_argument_group("Other options")

    # I/O options
    io_opts.add("-o",
        "--out_dir",
        dest = "output_dir",
        required = True,
        help = "Location of output directory, which should already exist",
        type = lambda x: is_valid_dir(main_parser, x)
    )
    io_opts.add("-g",
        "--genotypes",
        dest = "geno_f",
        required = True,
        help = "TSV file containing genotypes across strains (can be " +
        "integers or continuous values)",
        type = argparse.FileType("r"),
        default = None
    )
    io_opts.add("-c",
        "--counts",
        dest = "counts_f",
        required = True,
        help = "TSV file containing normalised expression values per feature",
        type = argparse.FileType("r"),
        default = None
    )
    io_opts.add("-m",
        "--metadata",
        dest = "meta_f",
        required = False,
        help = "TSV file containing metadata for samples. Columns contain " + 
        "model covariates.",
        type = argparse.FileType("r"),
        default = None
    )

    # subparsers
    subparsers = main_parser.add_subparsers(
        dest="command", 
        title="Available commands"
    )

    rf_subparser = subparsers.add_parser(
        "random-forest",
        help="random forest model"
    )
    rf_subparser = rf_parser(rf_subparser)

    # parse arguments and run function
    args = main_parser.parse_args()
    args.func(args)

    return


if __name__ == "__main__":
    main()
