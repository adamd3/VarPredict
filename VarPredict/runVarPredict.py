#!/usr/bin/env python3

import os
import copy
import numpy as np

try:
    import CIAlign.nested_cv as nested_cv
except ImportError:
    import nested_cv

def run(args):

    arr, nams, typ = setupArrays(args, log)
    functions = getFuncs(args)

    if "matrices" in functions:
        # Make similarity matrices
        runMatrix(args, log, orig_arr, orig_nams, arr, nams)


def getFuncs(args):
    '''
    Make a list to track which groups of functions to run
    Parameters
    ----------
    args: configargparse.ArgumentParser
        ArgumentParser object containing the specified parameters
    Returns
    -------
    which_functions: list
        List of groups of functions to run
    '''
    which_functions = []
    # Cleaning Functions
    if any([args.remove_divergent,
            args.remove_insertions,
            args.crop_ends,
            args.remove_short,
            args.remove_gaponly,
            args.crop_divergent,
            args.clean,
            args.all_options]):
        which_functions.append("cleaning")

    # Similarity Matrix
    if any([args.make_simmatrix_input,
            args.make_simmatrix_output,
            args.interpret,
            args.all_options]):
        which_functions.append("matrices")