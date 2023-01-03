

## AD: this script can be deleted


# # Run the VarPredict pipeline

# import os
# import argparse
# from .isvalid import *
# from .rf_model import *

# from .__init__ import __version__


# def vpd_parser(parser):

#     parser.description = "Run VarPredict on a set of variants and covariates"

#     # Argument groups
#     io_opts = parser.add_argument_group("Input/output files")
#     hyper_opts = parser.add_argument_group("Hyperparameter options")
#     other_opts = parser.add_argument_group("Other options")

#     # I/O arguments
#     io_opts.add("-o",
#             "--out_dir",
#             dest = "output_dir",
#             required = True,
#             help = "Location of output directory, which should already exist",
#             type = lambda x: is_valid_dir(parser, x))
#     io_opts.add("-g",
#             "--genotypes",
#             dest = "geno_f",
#             help = "TSV file containing genotypes across strains (can be integers or continuous values)",
#             type = argparse.FileType("r"),
#             default = None)
#     io_opts.add("-c",
#             "--counts",
#             dest = "counts_f",
#             help = "TSV file containing normalised expression values per feature",
#             type = argparse.FileType("r"),
#             default = None)
#     io_opts.add("-m",
#             "--metadata",
#             dest = "meta_f",
#             help = "TSV file containing metadata for samples",
#             type = argparse.FileType("r"),
#             default = None)



#     # Other options
#     other_opts.add("--version",
#             action = "version",
#             version = "%(prog)s " + __version__,
#             help = "Display current version")
            

#     parser.set_defaults(func=rf_model)

#     return parser

