from VarPredict.__init__ import __version__
import configargparse
import os.path


def buildParser():
    '''
    Constructs a configargparse.ArgumentParser object with the chosen parameters
    Returns
    -------
    parser: configargparse.ArgumentParser
        configargparse.ArgumentParser with the chosen parameters
    '''

    parser.description = "Run VarPredict on a set of bacterial genomic variants"


    parser = configargparse.ArgumentParser(
            description='Model gene expression using variants and covariates', 
            add_help=False)
    vpd_dir = os.path.dirname(utilityFunctions.__file__)

    # Default params
    def_vals = [line.strip().split("\t")
              for line in open("%s/defaults.txt" % vpd_dir)]
    defaults = {x[0]: x[1] for x in def_vals}

    # Argument groups
    io_opts = parser.add_argument_group("Input/output files")
    hyper_opts = parser.add_argument_group("Hyperparameter options")
    other_opts = parser.add_argument_group("Other options")

    # I/O arguments
    io_opts.add("-o",
            "--out_dir",
            dest = "output_dir",
            required = True,
            help = "Location of output directory, which should already exist",
            type = lambda x: is_valid_dir(parser, x))
    io_opts.add("-g",
            "--genotypes",
            dest = "geno_f",
            help = "TSV file containing genotypes across strains (can be integers or continuous values)",
            type = argparse.FileType("r"),
            default = None)
    io_opts.add("-c",
            "--counts",
            dest = "counts_f",
            help = "TSV file containing normalised expression values per feature",
            type = argparse.FileType("r"),
            default = None)
    io_opts.add("-m",
            "--metadata",
            dest = "meta_f",
            help = "TSV file containing metadata for samples",
            type = argparse.FileType("r"),
            default = None)

    # Hyperparameter arguments
    hyper_opts.add("--max_depth",
            dest = "max_depth",
            default = defaults['max_depth'],
            type = float,
            help = "The number of splits that each decision tree is allowed to make"
            )
    hyper_opts.add("--max_features",
            dest = "max_features",
            default = defaults['max_features'],
            type = float,
            help = "The number of features to consider when looking for the best split"
            )
    hyper_opts.add("--min_samples_split",
            dest = "min_samples_split",
            default = defaults['min_samples_split'],
            type = float,
            help = "The minimum number of samples required to split an internal node"
            )
    hyper_opts.add("--n_estimators",
            dest = "n_estimators",
            default = defaults['n_estimators'],
            type = float,
            help = "The number of trees in the forest"
            )

    # Other
    other_opts.add("--version",
            action = "version",
            version = "%(prog)s " + __version__,
            help = "Display current version")
    return (parser)
