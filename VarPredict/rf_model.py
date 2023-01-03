# Random Forest model

import os
import argparse

from .__init__ import __version__


def rf_parser(parser):

    parser.description = "Random Forest model"

    # Hyperparameter arguments
    hyper_opts = parser.add_argument_group("Hyperparameter options")
    hyper_opts.add("--max_depth",
            dest = "max_depth",
            nargs = '+',
            default = [5,10, None],
            help = "The number of splits that each decision tree is allowed to make"
            )
    hyper_opts.add("--max_features",
            dest = "max_features",
            nargs = '+',
            default = ['sqrt', 'log2', None],
            help = "The number of features to consider when looking for the best split"
            )
    hyper_opts.add("--min_samples_split",
            dest = "min_samples_split",
            nargs = '+',
            default = [2, 4, 8],
            help = "The minimum number of samples required to split an internal node"
            )
    hyper_opts.add("--n_estimators",
            dest = "n_estimators",
            nargs = '+',
            default = [300],
            help = "The number of trees in the forest"
            )


    parser.set_defaults(func=rf_model)

    return parser

def rf_model(args):
    #  Construct RF models
    print("Running Random Forest modelling")

    ## AD (UPTOHERE): 
    ## Probably don't need a separate script with 
    ##           a `run_sift` (`run_rf`) function; instead just 
    ##          read in the arguments passed as below

    ## ALSO: ADD FILTERING SCRIPT (MAIN OPTION)


    # run_sift(fasta_file = args.fasta_f,
    #     vars_file = args.vars_f,
    #     gff_file = args.gff_f,
    #     # counts_file = args.counts_f,
    #     out_dir = args.output_dir
    # )

    return

def main():
    # set up and parse arguments
    parser = argparse.ArgumentParser()
    parser = rf_parser(parser)
    args = parser.parse_args()

    # run modelling
    args.func(args)

    return

if __name__ == "__main__":
    main()