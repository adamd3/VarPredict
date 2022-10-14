# Runs BactVarMiner on a set of bacterial genomic variants

import os
import argparse
from .is_valid import *
from .parse_gff import *
from .sift import *

from .__init__ import __version__

def dbv_parser(parser):

    parser.description = "Run BactVarMiner on a set of bacterial genomic variants"

    # Arguments for input and output files
    io_opts = parser.add_argument_group("Input/output")
    io_opts.add_argument("-o",
                        "--out_dir",
                        dest = "output_dir",
                        required = True,
                        help = "Location of output directory, which should already exist",
                        type = lambda x: is_valid_dir(parser, x))
    io_opts.add_argument("-v",
                        "--vars",
                        dest = "vars_f",
                        help = "TSV file containing amino acid variants across strains",
                        type = argparse.FileType("r"),
                        default = None)
    io_opts.add_argument("-g",
                        "--gff",
                        dest = "gff_f",
                        help = "GFF annotation file containing gene coordinates in reference sequence",
                        type = argparse.FileType("r"),
                        default = None)
    io_opts.add_argument("-f",
                        "--fasta",
                        dest = "fasta_f",
                        help = "Multi-FASTA file containing amino acid sequences for the reference strain",
                        type = argparse.FileType("r"),
                        default = None)
    # io_opts.add_argument("-c",
    #                     "--counts",
    #                     dest = "counts_f",
    #                     help = "TSV file containing normalised gene expression values across strains. " +
    #                     "Must not contain missing values",
    #                     type = argparse.FileType("r"),
    #                     default = None)
    # io_opts.add_argument("--example_string",
    #                     dest = "example_string",
    #                     help = "String arg",
    #                     type = str,             < for Strings
    #                     default = None)
    # io_opts.add_argument("--example_bool",
    #                     dest = "example_bool",
    #                     help = "Boolean arg",
    #                     action = "store_true",  < for Booleans
    #                     default = False)
    # Other arguments
    parser.add_argument("-t",
                        "--threads",
                        dest = "threads",
                        help = "Number of threads to use. Default = 8",
                        default = "8")
    parser.add_argument("--version",
                        action = "version",
                        version = "%(prog)s " + __version__)

    parser.set_defaults(func=BactVarMiner)

    return parser



def BactVarMiner(args):

    print("Running SIFT")

    run_sift(fasta_file = args.fasta_f,
        vars_file = args.vars_f,
        gff_file = args.gff_f,
        # counts_file = args.counts_f,
        out_dir = args.output_dir)

    return



def main():
    # set up and parse arguments
    parser = argparse.ArgumentParser()
    parser = dbv_parser(parser)
    args = parser.parse_args()

    # run BactVarMiner
    args.func(args)

    return

if __name__ == "__main__":
    main()
