# Runs VarPredict 

import os
import argparse
from .is_valid import *
from .parse_gff import *
from .sift import *

from .__init__ import __version__

def VarPredict(args):

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
    parser = ncv_parser(parser)
    args = parser.parse_args()

    # run VarPredict
    args.func(args)

    return

if __name__ == "__main__":
    main()
