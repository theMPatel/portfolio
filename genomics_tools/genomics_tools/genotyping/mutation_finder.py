###################################################################
#
# Caller script for the mutation detection workflow
# 
# Author: Milan Patel
# Contact: mpatel5@cdc.gov
# Version 1.0
#
###################################################################

# Set the version here. It makes more sense to set it as a variable
# here because if you edit the file, you can edit the version at the
# the same time

from tools.environment import (
    log_message,
    log_error,
    log_progress,
    log_algo_version,
    write_results
)

from tools.dbinfo import (
    DbInfo,
    SequenceInfo
)

from tools.tools import (
    is_fasta,
    parse_fasta
)

from .ab_detection import (
    mutation_detector
)

import os
import json
from functools import partial
from collections import namedtuple, defaultdict

MutationTarget = namedtuple('MutationTarget', [
    'gene_id', 
    'gene_name', 
    'num_mutations_needed',
    'codon_position',
    'reference_codon', 
    'reference_aa', 
    'resistance_aa', 
    'resistance', 
    'pm_ids',
    'coding_gene'
])

def sequence_parser(header, sequence, sep='_'):
    return SequenceInfo(
        locus = header,
        allele = header,
        accession = '',
        sequence = sequence,
        other = '')

def main(settings, env):

    # Log the initial message
    log_message('Starting running mutation finder algorithm')

    # Write the version number of the database and algorithm
    version_path = settings['version']

    # Set the initial version information
    log_algo_version(
        algo_version = None,
        settings = settings,
        env = env
    )

    # Get the database path
    database_path = env.get_sharedpath(settings.database)

    # Log it
    log_message('Database path found at: {}'.format(
        database_path))

    # Let's load it all
    log_message('Loading resistance sequences and associated'
        ' information')

    sequence_database = DbInfo(
        database_path, seq_parser = sequence_parser)

    # Load the mutation targets
    sequence_database.load_extras()

    log_message('Successfully loaded sequences and metadata!')
    log_message('Running mutation finder pipeline...')

    # The results will come back without being filtered
    results = mutation_detector(
        sequence_database,
        settings.query,
        settings.percent_identity,
        settings.min_relative_coverage,
        env
    )

    final_results, antibios_out = sequence_database.results_parser(results, f=results_parser)

    log_message('Predicting {} mutations of interest'.format(
        sum(final_results['results'].values())))

    log_message('Writing results out...')

    write_results('resistance.point.json', json.dumps(final_results))

    # Success!
    log_message('Successfully ran mutation finder algorithm!')

    return antibios_out

class DbInfo(DbInfo):

    def load_database(self, dirpath, seq_parser, note_parser):

        # Get all of the sequences in the directory
        for seq_file in os.listdir(dirpath):

            # Create the file path
            file_path = os.path.join(dirpath, seq_file)

            # Check to make sure that this is a valid fasta
            # file
            if not is_fasta(file_path):
                continue

            # Load the sequences themselves from the files
            sequences = parse_fasta(file_path)

            for seq_id, sequence in sequences.iteritems():
                # Create a named tuple that contains the 
                # header split out into its different
                # attributes plus the sequence
                seq_info = seq_parser(seq_id, sequence)

                self._sequences[seq_id] = seq_info

        # Load notes if they exist
        file_path = os.path.join(dirpath, 'notes.txt')

        if os.path.exists(file_path):

            with open(file_path, 'r') as f:

                for line in f:

                    line = line.strip()

                    if not line or line[0] == '#':
                        continue

                    notes_info = note_parser(line)

                    self._notes[notes_info.locus] = notes_info

    def load_extras(self):

        self._targets = defaultdict(list)
        self._rna_genes = set()

        # Load the RNA gene file if it exists
        file_path = os.path.join(self._dirpath, 'RNA_genes.txt')

        if os.path.exists(file_path):

            with open(file_path, 'r') as f:

                for line in f:

                    line = line.strip()

                    self._rna_genes.add(line)

        file_path = os.path.join(self._dirpath, 'resistens-overview.txt')

        if not os.path.exists(file_path):
            raise RuntimeError('Missing mutations file for interpretations..')

        with open(file_path, 'r') as f:

            for line in f:

                line = line.strip()

                # Blank or comment line
                if not line or line[0] == '#':
                    continue

                parts = line.split('\t')

                gene_id = parts[0]
                gene_name = parts[1]

                num_mutations_needed = int(parts[2])

                # You need to multiply this number by
                # 3... codon duh
                codon_position = int(parts[3])

                # These don't need splitting per se
                # however you might as well since
                # you never know how the file will look
                # in the future
                reference_codon = parts[4].split(',')
                reference_aa = parts[5].split(',')
                resistance_aa = parts[6].split(',')
                resistance = parts[7].replace('resistance', '')
                resistance = map(str.strip, resistance.split(','))
                pm_ids = parts[8].split(',')

                coding_gene = gene_id not in self._rna_genes and \
                    'promoter' not in gene_id.lower()

                self._targets[gene_id].append(
                    MutationTarget(
                        gene_id = gene_id,
                        gene_name = gene_name,
                        num_mutations_needed = num_mutations_needed,
                        codon_position = codon_position,
                        reference_codon = reference_codon,
                        reference_aa = reference_aa,
                        resistance_aa = resistance_aa,
                        resistance = resistance,
                        pm_ids = pm_ids,
                        coding_gene = coding_gene
                    )
                )

    @property
    def targets(self):
        return self._targets

def results_parser(dbinfo, interpretations):
    # Filter out the ones that don't meet basic coverage requirements
    final_results = {
        'results': {},
        'extra': []
    }

    notes_out = {
        'results' : {},
        'extra' : []
    }

    log_message('Determining optimal gene coverages...')
    for gene, mutation in interpretations.iteritems():

        for mutation_info in mutation:

            gene_name = gene + '@' + str(mutation_info['position'])

            final_results['results'][gene_name] = True

            hit = mutation_info['hit']

            hit_info = {
                        'locus': hit.reference_id,
                        'identity': hit.identity,
                        'coverage': hit.relative_len,
                        'contig_id': hit.query_id,
                        'query_start': hit.query_start,
                        'query_stop': hit.query_stop,
                        'reference_start': hit.reference_start,
                        'reference_stop': hit.reference_stop,
                        'aa_mutation': mutation_info['aa_mutation'],
                        'query_codon': mutation_info['query_codon'],
                        'iscoding' : mutation_info['iscoding'],
                        'resistance' : mutation_info['resistance'],
                        'position' : mutation_info['position']
                }

            final_results['extra'].append(hit_info)
            notes_out['extra'].append(hit_info)

            for r in mutation_info['resistance']:
                notes_out['results'][r] = True

    
    for gene, mutation_targets in dbinfo.targets.iteritems():

        for mutation_target in mutation_targets:

            gene_name = '{gene}@{location}'

            final_name = gene_name.format(
                gene = mutation_target.gene_id,
                location = mutation_target.codon_position
            )

            final_results['results'][final_name] = final_results['results'].get(final_name, False)

            for r in mutation_target.resistance:
                notes_out['results'][r] = notes_out['results'].get(r, False)


    return final_results, notes_out