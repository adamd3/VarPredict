import os, sys
import argparse

from .__init__ import __version__

from .rf_model import rf_parser
# from .en_model import en_parser

class UltimateHelpFormatter(
    argparse.RawTextHelpFormatter, argparse.ArgumentDefaultsHelpFormatter):
    pass


def main():

    # parent subparser for shared args
    # (see https://stackoverflow.com/questions/7498595/python-argparse-add-argument-to-multiple-subparsers)

    main_parser = argparse.ArgumentParser(formatter_class=UltimateHelpFormatter)

    # main arguments
    main_parser.add_argument(
        "-f",
        "--filter",
        dest = "filter",
        help = "filter gene expression data to remove genes with low variance",
        action = "store_true",
        default = False
    )
    main_parser.add_argument(
        "-v",
        "--version", 
        action = "version", 
        version = "%(prog)s " + __version__
    )

    # subparsers for ML models
    subparsers = main_parser.add_subparsers(
        dest = "model", help = "specify model type",
        title = "available models", required = True
    )

    rf_subparser = subparsers.add_parser(
        "random-forest", 
        parents = [main_parser], # if not specified, common options not used
        help = "random forest model", 
        add_help = False
    )
    rf_subparser = rf_parser(rf_subparser)

    en_subparser = subparsers.add_parser(
        "elastic-net", 
        parents = [main_parser], # if not specified, common options not used
        help = "elastic net-regularised logistic regression model", 
        add_help = False
    )
    en_subparser = en_parser(en_subparser)


    # parse arguments and run function
    args = main_parser.parse_args()
    args.func(args)

    return



if __name__ == "__main__":
    main()
