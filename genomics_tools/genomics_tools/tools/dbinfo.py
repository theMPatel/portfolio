###################################################################
#
# Methods and classes for reference databases
# 
# Author: Milan Patel
# Contact: mpatel5@cdc.gov
# Version 1.0
#
###################################################################

import os
from collections import namedtuple, defaultdict

from .tools import (
    parse_fasta,
    is_fasta,
)

from .environment import (
    check_dir,
    valid_dir
)

# 'Structs' for datastorage
SequenceInfo = namedtuple('SequenceInfo', [
    'locus', 'allele', 'accession', 'sequence', 'other'])

LocusInfo = namedtuple('LocusInfo', [
    'locus', 'note', 'antibiotic', 'other'])

def sequence_parser(header, sequence, sep = ':'):
    # Default parser for the database files
    parts = header.split(sep)
    parts = list(map(str.strip, parts))

    while len(parts) < 4:
        parts.append('')

    return SequenceInfo(
        locus = parts[0],
        allele = parts[1],
        accession = parts[2],
        sequence = sequence,
        other = parts[3]
    )

def notes_parser(line, sep=':'):
    # Default parser for the notes files
    parts = list(map(str.strip, line.split(sep)))

    while len(parts) < 3:
        parts.append('')

    antibiotic = parts[1].replace('resistance', '')
    antibiotic = map(str.strip, antibiotic.split(','))

    return LocusInfo(
        locus = parts[0],
        note = parts[1],
        antibiotic = antibiotic,
        other = parts[2]
    )

class DbInfo(object):
    # Class that will hold the db information
    def __init__(self, dirpath, seq_parser = sequence_parser,
        note_parser = notes_parser):

        self._notes = {}
        self._sequences = {}
        self._dirpath = dirpath
        self._separator = None

        if dirpath is None or not check_dir(dirpath):
            raise RuntimeError('Invalid path provided for '
                ' database fastas: {}'.format(str(dirpath)))

        self.load_database(dirpath, seq_parser, note_parser)

    def load_database(self, dirpath, seq_parser, note_parser):

        sequence_counts = defaultdict(dict)
        allele_id_template = '{}_{}'
        allele_id_template_i = '{}_{}-{}'
        self._separator = '-'

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

                # Get the allele information
                allele_id = allele_id_template.format(
                    seq_info.locus, seq_info.allele)

                if allele_id in sequence_counts:
                    # There are a lot of alleles with the same name
                    # new is created to make each sequence unique
                    new_id = allele_id_template_i.format(
                        seq_info.locus,
                        seq_info.allele, 
                        len(sequence_counts[allele_id])
                    )

                    # Set the seq_info, seq_counts just
                    # functions to store the counts
                    sequence_counts[allele_id][new_id] = True
                    self._sequences[new_id] = seq_info

                else:

                    # # Create the new id, 0-indexed
                    # new_id = allele_id_template_i.format(
                    #     seq_info.locus,
                    #     seq_info.allele, 
                    #     0
                    # )

                    # Set the seq_info, seq_counts just
                    # functions to store the counts
                    sequence_counts[allele_id][allele_id] = True
                    self._sequences[allele_id] = seq_info

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

    @property
    def sequences(self):
        return self._sequences

    @property
    def notes(self):
        return self._notes

    def export_sequences(self, filepath):

        # Make sure the directory exists
        valid_dir(os.path.dirname(filepath))

        # Open
        with open(filepath, 'w') as f:

            for seq_id, seq_info in self._sequences.iteritems():
                
                # The fasta file should look like:
                # 
                # >allele_id|1234
                # ACGTACGTACGTACGTACGTACGTACGTACGT
                # ACGTACGTACGTACGTACGTACGTACGTACGT

                ostr = '>{}|{}\n{}\n'

                f.write(ostr.format(
                    seq_id,
                    len(seq_info.sequence),
                    seq_info.sequence
                    )
                )

    def get_refseq(self, ref):
        if ref in self._sequences:
            return self._sequences[ref].sequence

        raise KeyError('Missing reference: {}'.format(ref))

    def results_parser(self, results, f=None):

        if f is not None and callable(f):
            return f(self, results)

        sequences = self.sequences
        notes = self.notes

        if not len(sequences):
            raise RuntimeError('Trying to assimilate results for'
                ' something with no sequences')

        results_out = {
            'results': {},
            'extra': []
        }

        for result, geno_regions in results.iteritems():

            # The sequence information for this result
            sequence_info = sequences[result]
            
            # The locus of this sequence
            locus = sequence_info.locus

            # The allele of the found genotype, if there is one
            allele = sequence_info.allele

            # gene_name = '_'.join([locus, allele])
            # Changing the output of the character since only the locus
            # is important *for surveillance*
            gene_name = locus

            # Add the locus to the results
            results_out['results'][gene_name] = True

            # This will be the extra information that
            # will get dumped into the entry logs
            # no need to make parsing the results
            # more difficult than it has to be
            hit_information = [
                {
                        'locus': locus,
                        'identity': geno_region.identity,
                        'coverage': geno_region.coverage,
                        'allele': allele,
                        'hits': [
                                    {
                                    'contig_id': hit.query_id,
                                    'query_start': hit.query_start,
                                    'query_stop': hit.query_stop,
                                    'reference_start': hit.reference_start,
                                    'reference_stop': hit.reference_stop,
                                    'full_match' : hit.full_match,
                                    } for hit in geno_region.locations
                                ]
                } for geno_region in geno_regions
            ]

            # Add it to the results out
            results_out['extra'].extend(hit_information)

        for sequence_info in sequences.itervalues():

            locus = sequence_info.locus
            allele = sequence_info.allele

            # gene_name = '_'.join([locus, allele])
            gene_name = locus

            if gene_name not in results_out['results']:
                results_out['results'][gene_name] = False

            # results_out['results'][gene_name] = gene_name in results_out['results']

        return results_out
    


