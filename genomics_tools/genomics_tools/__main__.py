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
    Environment, initialize_logging,
    graceful_shutdown_logging,
    ResultWriter, log_message,
    log_progress, log_error,
    log_exception, log_algo_params,
    get_stack_len, set_base_depth,
    populate_syspath, find_directory_on_path
)

from genotyping import mutation_finder

MutationFinderSettings = namedtuple("Settings", [
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

    parser = MyArgumentParser()

    parser.add_argument('--nThreads',
        help='Number of threads', type=int, default=2)

    parser.add_argument('--tempdir',
        help='Temporary shared directory', type=str)

    parser.add_argument('--resultsdir',
        help='Results directory', type=str)

    parser.add_argument('--databasedir',
        help='Shared directory', type=str)

    parser.add_argument("--run", default=False, action="store_true")

    args, remaining = parser.parse_known_args()

    # I will instruct you to use --run to do a test run of the software!
    if args.run:

        tempdir = tempfile.mkdtemp()
        args.tempdir = os.path.join(tempdir, "tmp")
        args.resultsdir = os.path.join(tempdir, "results")
        os.makedirs(args.resultsdir)

    #return the arguments object
    return args, remaining

def main_throw_args(args, remaining):
    # Actually calls the genotyping algorithm

    env = Environment()
    env.setup(vars(args))
    initialize_logging(env.logdir)
    ResultWriter(env.resultsdir)
    
    log_message('Initializing..')
    
    path = os.environ["PATH"]
    database_dir = find_directory_on_path("pointfinder_db", path)
    sequence_dir = find_directory_on_path("sequence_data", path)
    ecoli_db_dir = os.path.join(database_dir, "escherichia_coli")
    
    base_settings = MutationFinderSettings(query="", version="1.0.0",
                                            database=ecoli_db_dir,
                                            percent_identity=0.9,
                                            min_relative_coverage=0.6)

    log_message("Using temp directory: {}".format(env.tempdir))
    log_message("Using results directory: {}".format(env.resultsdir))
    log_message("Using database: {}".format(database_dir))
    log_message("Using sequences from: {}".format(sequence_dir))
    log_progress(0)

    # Set's the base depth for logging so that we can get tabbed log
    # files that mimic the execution flow of the program
    set_base_depth(-(get_stack_len()))

    for file in os.listdir(sequence_dir):
        query_path = os.path.join(sequence_dir, file)

        settings = MutationFinderSettings(
            query=query_path,

            database=base_settings.database,
            version=base_settings.version,
            percent_identity=base_settings.percent_identity,
            min_relative_coverage=base_settings.min_relative_coverage
        )
    
        # Actually, the below is a major refactor of the original
        # project layout. Before we were actually dynamically importing
        # the modules that a client wanted to run. This particular file
        # served as the insertion point for any and all modules the client
        # was aware it could run and which one in particular would come
        # in as a cmdline argument. We would import that class dynamically
        # and run it's main method.
        mutation_finder.main(settings, env)
        log_message("")
        
    # And we are done!
    log_progress(100)
    log_message('Done running algorithm: {}!'.format(
        "Mutation Finder"))    

_tools_dirs = [
    os.path.normpath(os.path.expanduser("~/.tools")),
    os.path.normpath(os.path.expanduser("~/.bin")),
    os.path.normpath(os.path.expanduser("~/.tmp"))
]
def main_throw():
    """
    Staging grounds for the application. The intention is to
    do some ancillary tasks prior to running the program.
    """
        
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
    
    # Update the path to add all the directories that might have
    # our tools.
    populate_syspath(_tools_dirs)

    # Run the main program with arguments
    main_throw_args(args, remaining)

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
