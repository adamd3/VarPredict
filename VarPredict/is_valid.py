# Functions to check if files and directories are valid

import os

def is_valid_dir(parser, arg):
    if not os.path.isdir(arg):
        parser.error("The folder %s does not exist, please create the folder then rerun" % arg)
    else:
        return(arg)
