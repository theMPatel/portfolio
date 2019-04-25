###################################################################
#
# Insertion point for the genomics tools. We'll call this file
# and setup the appropriate env contexts, execute the appropriate
# modules and return back to the CLI what our success status was.
# 
# This tool is a sample and distillation of the real application
# hosted at: https://github.com/theMPatel/functional_genomics_tools
#
# Author: Milan Patel
# Contact: https://github.com/theMPatel
# Version 1.0
#
###################################################################

import argparse
from collections import namedtuple 
import importlib
import json
import os
import sys
import tempfile
import traceback

this_file = os.path.realpath(__file__)
base_path = os.path.dirname(this_file)

if base_path not in sys.path: 
    sys.path.append(base_path)

from tools.environment import (
    Environment,
    initialize_logging,
    graceful_shutdown_logging,
    ResultWriter,
    log_message,
    log_progress,
    log_error,
    log_exception,
    log_algo_params
)

from genotyping import mutation_finder

Settings = namedtuple("Settings", [
    'query',
    'version',
    'database',
    'percent_identity',
    'min_relative_coverage'
])


# When you run the python interpreter without the -O flag
# the below __debug__ variable will be true.
if __debug__:
    import pdb

class MyArgumentParser(argparse.ArgumentParser):

    # Override to split space based cmdline args that
    # will come in through a file.
    def convert_arg_line_to_args(self, arg_line):
        return arg_line.split()

# Parse the commandline arguments
def parse_cmdline():

    # create the parser, the below will look for any file
    # that starts with an '@' symbol. It makes sense to do this
    # when your command line arguments are so large that it
    # becomes more useful to manage them inside of a file.
    parser = MyArgumentParser(fromfile_prefix_chars='@')

    parser.add_argument('--nThreads',
        help='Number of threads', type=int, default=2)

    parser.add_argument('--tempdir',
        help='Temporary shared directory', type=str)

    parser.add_argument('--resultsdir',
        help='Results directory', type=str)

    parser.add_argument('--databasedir',
        help='Shared directory', type=str)

    parser.add_argument('--toolsdir',
        default='', help='Tools directory', type=str)

    parser.add_argument("--run", default=False, action="store_true")

    # Get the known arguments of the command line, by default
    # prefer the arguments defined in the settings file
    # over any actual cmd-line arguments
    try:
        args, remaining = parser.parse_known_args(args=['@settings.txt'])

    except:
        # In the case that we don't get a settings file, we can see
        # if the arguments we need are defined on the command line.
        args, remaining = parser.parse_known_args()

    # I will instruct you to use --run to do a test run of the software!
    if args.run:

        args.tempdir = tempfile.mkdtemp()
        args.resultsdir = os.path.join(args.tempdir, "results")
        os.makedirs(args.resultsdir)
        print(args.tempdir)

    #return the arguments object
    return args, remaining

def populate_syspath():
    """
    This function will help us find any of the tools that we
    installed through the setup.py script. It exists because
    I don't want to permanently modify your own $PATH variable
    """
    _tools_dirs = [
        os.path.normpath(os.path.expanduser("~/.tools")),
        os.path.normpath(os.path.expanduser("~/.bin")),
        os.path.normpath(os.path.expanduser("~/.tmp"))
        ]

    path_parts = []
    # We can use shutil.which to get the correct binary we want
    # if we ensure that the path is correct. Unfortunately, I
    # didn't feel safe manipulating your path (esp if you're on
    # windows) so am resorting to this.
    for directory in _tools_dirs:
        for root, dirs, files in os.walk(directory):
            path_parts.append(root)

    os.environ["PATH"] += os.pathsep + os.pathsep.join(path_parts)

def main_throw_args(args, remaining, settings):
    # Actually calls the genotyping algorithm
    
    env = Environment()
    env.setup(vars(args))

    initialize_logging(env.logdir)
    ResultWriter(env.resultsdir)

    log_message('Initializing..')
    log_algo_params(vars(args))
    log_progress(0)

    # Actually, the below is a major refactor of the original
    # project layout. Before we were actually dynamically importing
    # the modules that a client wanted to run. This particular file
    # served as the insertion point for any and all modules the client
    # was aware it could run and which one in particular would come
    # in as a cmdline argument. We would import that class dynamically
    # and run it's main method.
    mutation_finder.main(settings, env)
        
    # And we are done!
    log_progress(100)
    log_message('Done running algorithm: {}!'.format(
        "Mutation Finder"))    

def main_throw():
        
    # Parse cmdline arguments
    args, remaining = parse_cmdline()

    # In the actual project, we had a shared network drive
    # that all nodes on the compute cluster had access to.
    # One of the things on the compute was the default settings
    # for the algorithms to be run, and the below would be
    # placed in the run directory of a job which contained
    # the settings that were requested by a remote client.
    # These settings overrode anything defined in the default
    # settings.
    genotyper_settings_path = os.path.join(
        os.getcwd(), 'genotyper_settings.json')

    mutation_finder_settings = {}
    if os.path.exists(genotyper_settings_path):
        with open(genotyper_settings_path, 'r') as f:
            genotyper_settings = json.load(f)

    settings = Settings(os.path.expanduser("~/Downloads/test_genome.fna"),
                        "1.0.0",
                        os.path.expanduser("~/Documents/github/pointfinder_db/escherichia_coli"),
                        0.9,
                        0.6)
    
    # Update the path to add all the directories that might have
    # our tools.
    populate_syspath()

    # Run the main program with arguments
    main_throw_args(args, remaining, settings)

def _main():

    return_code = 0

    try:
        main_throw()

    except Exception:
        log_exception('')
        return_code = 1

    finally:
        graceful_shutdown_logging()
        sys.exit(return_code)
