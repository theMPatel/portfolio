###################################################################
#
# The genotyping algorithm that calls all the scripts
# 
# Author: Milan Patel
# Contact: mpatel5@cdc.gov
# Version 1.0
#
###################################################################

import os
import sys
import logging
import argparse
import importlib
import traceback

from tools.tools import (
    unzip_file,
    is_fasta,
    fasta_iterator,
    parse_fasta,
    process_seq_file,
    process_read_files
)

import tools.environment

from tools.environment import (
    log_message,
    log_error,
    log_exception,
    log_progress,
    log_algo_version,
    log_algo_params,
    sanitize_path,
    get_stack_len,
    set_base_depth
)

from tools.config import Config
from tools.custom_parser import CustomParser

def parse_settings(args, remaining):

    parser = argparse.ArgumentParser()
   
    # Algorithm parameters
    parser.add_argument('--organism', help = 'Organism', type=str)
    parser.add_argument('--query', help = 'Query genome', type=str) 
    
    # Holds the settings for each of the algorithms
    # as well as the paths to their respective database
    parser.add_argument('--configfile', 
        default = r'[SHAREDDIR]/Genotyping/config.json',
        help='Config file', type=str)

    # If the user requests a specific genotyper to run
    parser.add_argument('--genotypers', 
        help="The specific genotyping algorithms you want to run",
        action='append', type=str, default=[])

    parser.add_argument('--query-reads',
        help='Read files', action='append', type=str, default=[])

    # Get the most up-to-date stuffs
    update_args, remaining = parser.parse_known_args(remaining)

    # Get the two most necessary arguments
    update_args.query = sanitize_path(update_args.query)
    update_args.configfile = sanitize_path(update_args.configfile)
    update_args.query_reads = map(sanitize_path, update_args.query_reads)

    return update_args

def run_genotyper(module_name, settings, env):

    # Create the new module name
    module_full_name = 'genotyping.' + module_name

    # Import the modules
    module = importlib.import_module(module_full_name)

    # Get a new copy of the global environment
    module_env = env.copy()

    # Make a new tmp folder just for this special module
    module_env.localdir = os.path.join(module_env.localdir, module_full_name)

    # Set the base message_depth
    set_base_depth(-(get_stack_len()))
    
    # Run the module!
    try:
        module.main(settings, module_env)

    except:
        log_exception('Error running genotyper: {}'.format(module_name))

def setup_genotyper(genotyper, module_name, organism_config, env, data):

    # Check to make sure we actually got a module name
    if module_name is None:
        raise RuntimeError('Requested module does not'
            ' exist. This should not have occurred')

    # Get the settings of the genotyper
    genotyper_settings = organism_config.genotypers[genotyper]

    # Merge the custom args with the client requested arguments
    CustomParser.update(genotyper, genotyper_settings)

    # Add the query path to the settings
    genotyper_settings.query = data.get('query', None)

    # Add the reads to the query genotyper
    genotyper_settings.query_reads = data.get('query_reads', [])

    # Add the cached_query to the settings
    genotyper_settings.cached_query = data.get('cached_query', None)

    run_genotyper(module_name, genotyper_settings, env)

def main(args, remaining, env, module_settings):

    log_message('Initializing genotyping algorithm')
    
    # Log the algorithm version
    log_algo_version(
        algo_version = None,
        settings = module_settings,
        env = env
    )

    # Get the arguments that we need
    specific_args = parse_settings(args, remaining)

    log_algo_params(vars(specific_args))

    # Make sure the path is a __realpath__
    config_filepath = env.get_sharedpath(specific_args.configfile)

    # Get the global config
    global_config = Config(config_filepath)

    # Get the specfic organism config
    organism_config = global_config.organism_config(specific_args.organism)

    if organism_config is None:
        raise RuntimeError('Missing organism config'
            ' for organism: {}'.format(specific_args.organism))

    # Make sure the query genome is in fasta format
    # before running blast
    query_filename = specific_args.query
    
    log_message('Checking query file...')

    # Just so that logging is nice in this section
    set_base_depth(-(get_stack_len()))

    query_filename, cached_query = process_seq_file(query_filename, load=True)

    log_message('Checking read files...', extra=-1)
    
    unpacked_reads = process_read_files(specific_args.query_reads, load=False)
    specific_args.query_reads = [read[0] for read in unpacked_reads]

    # Turn off pretty logging
    set_base_depth(0)

    # The query is good!
    log_message('Query is ready to be analyzed')

    # Ready to rock n' roll!
    log_message('Performing genotyping analysis')

    # Determine the genotypers to run, if specifics were
    # selected
    all_genotypers = set(organism_config.genotypers.keys())

    # Requested genotypers
    activated_genotypers = set()
    client_args = CustomParser.args()

    # Below is a desgin choice, only run the things that the client is aware
    # it can run. Rather than running everything that it didn't request
    # This allows us the flexibility of adding modules without worrying that
    # the user will run it without having properly tested it.
    if client_args is not None:
        activated_genotypers.update(key for key in client_args.keys() \
            if client_args[key]['activated'])

    # The ones we will eventually run
    # Get the intersection of the genotypers
    genotypers_to_run = all_genotypers & activated_genotypers

    # For salmonella serotyping, we need to make sure
    # that insilico pcr has been run prior to doing serotype work
    # this is a cheap way to make sure that insilico pcr is always run
    # by the time serotyping for salmonella starts running
    genotypers_to_run = list(genotypers_to_run)
    genotypers_to_run.sort()

    data = {
        'query' : query_filename,
        'query_reads' : specific_args.query_reads, 
        'cached_query' : cached_query
    }

    for genotyper in genotypers_to_run:

        # Get the actual file name of the module
        module_name = global_config['modules'][genotyper]

        setup_genotyper(
            genotyper,
            module_name,
            organism_config,
            env,
            data
        )

    # Get the modules that we always need to run no matter what
    always_run = organism_config.always_run

    for genotyper in always_run:

        module_name = global_config['modules'][genotyper]

        setup_genotyper(
            genotyper,
            module_name,
            organism_config,
            env,
            data
        )