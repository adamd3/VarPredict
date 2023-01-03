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
        help = "filter gene expression data to remove genes with low variance",
        action = "store_true",
        default = False
    )
    main_parser.add_argument(
        "--version", 
        action = "version", 
        version = "%(prog)s " + __version__
    )

    # subparsers
    subparsers = main_parser.add_subparsers(
        dest = "model", help = "specify model type",
        title = "available models"
    )

    # parent subparser for shared args
    # note `add_help=False` and creation via `argparse.`
    # (see comment by Michael @ https://stackoverflow.com/questions/7498595/python-argparse-add-argument-to-multiple-subparsers)

    parent_parser = argparse.ArgumentParser(add_help = False)

    # input/output options
    parent_parser.add_argument("-o",
        "--out_dir",
        dest = "output_dir",
        required = True,
        help = "Location of output directory, which should already exist",
        type = lambda x: is_valid_dir(main_parser, x)
    )
    parent_parser.add_argument("-g",
        "--genotypes",
        dest = "geno_f",
        required = True,
        help = "TSV file containing genotypes across strains (can be " +
        "integers or continuous values)",
        type = argparse.FileType("r"),
        default = None
    )
    parent_parser.add_argument("-c",
        "--counts",
        dest = "counts_f",
        required = True,
        help = "TSV file containing normalised expression values per feature",
        type = argparse.FileType("r"),
        default = None
    )
    parent_parser.add_argument("-m",
        "--metadata",
        dest = "meta_f",
        required = False,
        help = "TSV file containing metadata for samples. Columns contain " + 
        "model covariates.",
        type = argparse.FileType("r"),
        default = None
    )
    

    rf_subparser = subparsers.add_parser(
        "random-forest", 
        parents=[parent_parser], # if not specified, common options not used
        help="random forest model"
    )
    rf_subparser = rf_parser(rf_subparser)

    # parse arguments and run function
    args = main_parser.parse_args()
    args.func(args)

    return

    ## AD example usage:
    ## VarPredict random-forest  -o [output] -g [genotypes] -c [counts] -m [metadata] \
    ##  -d [5,10, None] -f ['sqrt', 'log2', None] -s [2, 4, 8] -e [300]


if __name__ == "__main__":
    main()
