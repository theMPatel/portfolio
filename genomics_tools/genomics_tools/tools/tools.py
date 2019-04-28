###################################################################
#
# Tools for the genotyping algorithm.
# 
# This tool is a sample and distillation of the real application
# hosted at: https://github.com/theMPatel/functional_genomics_tools
#
# Author: Milan Patel
# Contact: https://github.com/theMPatel
# Version 1.0
#
###################################################################

import base64
import binascii
import gzip
import io
import json
import os
import shutil
import subprocess as sp
import sys

from collections import namedtuple
from itertools import combinations

from .environment import (
    full_path, log_message,
    log_warning, log_exception,
    log_error
)

# Get the codons for any particular amino acid
_AA_BACK_TRANSLATE = {
    'A' : ['GCT', 'GCC', 'GCA', 'GCG'],
    'C' : ['TGT', 'TGC'],
    'D' : ['GAT', 'GAC'],
    'E' : ['GAA', 'GAG'],
    'F' : ['TTT', 'TTC'],
    'G' : ['GGT', 'GGC', 'GGA', 'GGG'],
    'H' : ['CAT', 'CAC'],
    'I' : ['ATT', 'ATC', 'ATA'],
    'K' : ['AAA', 'AAG'],
    'L' : ['TTA', 'TTG', 'CTT', 'CTC', 'CTA', 'CTG'],
    'M' : ['ATG'],
    'N' : ['AAT', 'AAC'],
    'P' : ['CCT', 'CCC', 'CCA', 'CCG'],
    'Q' : ['CAA', 'CAG'],
    'R' : ['CGT', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG'],   
    'S' : ['TCT', 'TCC', 'TCA', 'TCG', 'AGT', 'AGC'],
    'T' : ['ACT', 'ACC', 'ACA', 'ACG'],
    'V' : ['GTT', 'GTC', 'GTA', 'GTG'],
    'W' : ['TGG'],
    'Y' : ['TAT', 'TAC'],
}

# Get the amino acid for any particular codon
_AA_TRANSLATE = {}

# Create the dictionary on the fly
for amino, codons in _AA_BACK_TRANSLATE.items():
    for codon in codons:
        _AA_TRANSLATE[codon] = amino

# All potential sequence based nucleotides
# and their pairings
_NUC_SIBLINGS = {
   'A': "A",
   'C': "C",
   'G': "G",
   'T': "T",
   'R': "AG",
   'Y': "CT",
   'S': "GC",
   'W': "AT",
   'K': "GT",
   'M': "AC",
   'B': "CGT",
   'D': "AGT",
   'H': "ACT",
   'V': "ACG",
   'N': "ACGT"
}

_SET_TO_NON_ACTG = {}

# Convert values to set for easy membership
# checking
for nuc, pairs in _NUC_SIBLINGS.items():
    fset = frozenset(pairs)
    _NUC_SIBLINGS[nuc] = fset
    _SET_TO_NON_ACTG[_NUC_SIBLINGS.get(nuc)] = nuc

def get_non_iupac(f_set):
    return _SET_TO_NON_ACTG.get(f_set, None)

_GZIP_START = b'1f8b'
def check_gzipped(file_path):

    # Read the first two bytes to see if it is the gzip
    # header
    with open(file_path, 'rb') as f:
        head = f.read(2)

    return binascii.hexlify(head) == _GZIP_START

def check_b64encoded(string):
    """
    Determines whether a string is base64 encoded or not.
    A base64 encoded string includes [a-zA-Z0-9] and
    +, / and = for padding. Another thing to note, if the
    length of your string modulo 4 is zero, and it only
    contains the characters above, it is base64 encoded.
    
    :param string: The string to check if base64 encoded.
    """
    try:
        base64.b64decode(string)
    
    except binascii.Error:
        return False

    else:
        return True

def get_all_file_exts(path):
    """
    Returns back all the file extensions of a file.
    Sometimes files like some_file.tar.gz have multiple
    extensions and it might be useful to check a file
    has a particular extension

    :param path: The string path to get all extensions
    """
    
    root, ext = os.path.splitext(path)

    if not ext:
        return root, []

    all_exts = [ext]

    real_root, exts = get_all_file_exts(root)
    exts += all_exts

    return real_root, exts

_FASTAEXTS = set(['.fna', '.fasta', '.fsa'])
def is_fasta(path):
    """
    Checks to see if a particular path could be
    pointing to a fasta file based on it's extension.

    :param path: The path to the file
    """
    if not path:
        return False

    root, exts = get_all_file_exts(path)

    return bool(set(exts) & _FASTAEXTS)

_GENBANKEXTS = set(['.gb', '.gbk'])
def is_genbank(path):
    """
    Checks to see if a particular path could be
    pointing to a genbank file based on it's
    extension.

    :param path: The path to the file
    """
    if not path:
        return False

    root, exts = get_all_file_exts(path)

    return bool(set(exts) & _GENBANKEXTS)

def fasta_iterator_path(path_to_file):
    """
    Interface to load a fasta file from an external
    path.

    :param path_to_file: The path to load
    :raises: OSError when the file does not exist
    """

    with open(path_to_file, 'r') as f:
        yield from fasta_iterator(f)

def fasta_iterator(fl_obj):
    """
    Files usually look like this:
    
    >fasta_id
    ACTGACTGACTGACTGACTGACTGACTG\n
    ACTGACTGACTGACTGACTGACTGACTG\n
    ACTGACTGACTGACTGACTGACTGACTG\n
    
    Parse accordingly. 
    
    :param fl_obj: A file object to load fasta
        entries from
    """

    # Stores the sequences
    sequence_parts = []

    # Key for a sequence
    key = ''

    for line in fl_obj:

        # Get rid of the newline
        line = line.strip()

        if not line:
            continue

        if line[0] == '>':
            if key:
                # Start of a new fasta record
                # return the previous record
                
                # join the sequence
                full_seqence = ''.join(sequence_parts).upper()
                
                yield (key, full_seqence)

            # Start of the file
            key = line[1:].split()[0]
            sequence_parts = []

        else:
            sequence_parts.append(line)

    if key:
        # The last sequence in the file
        full_seqence = ''.join(sequence_parts).upper()
        
        yield (key, full_seqence)

def parse_fasta(flname, rename=False):
    """
    Function provides an interface to parse
    a fasta file and optionally provide a
    way to rename the headers if requested

    :param flname: The path to the file
    :param rename: Whether to rename the header lines
    """
    
    if is_fasta(flname):
        fasta_sequences = {}

        if rename:
            for i, (name, sequence) in enumerate(
                    fasta_iterator_path(flname),1):

                new_name = 'contig_' + str(i)
                fasta_sequences[new_name] = sequence

        else:
            fasta_sequences.update(
                fasta_iterator_path(flname))

        return fasta_sequences

    raise RuntimeError('Requested file path is not of type fasta')

def codon_translation(codon):
    """
    Translates the codon that was requested to its
    amnio counterpart

    :param codon: The 3 letter DNA sequence to translate
    """

    if len(codon) != 3:
        raise ValueError('Cannot translate a codon: {} of '
            'non-standard size: {}'.format(
                codon, len(codon)))

    return _AA_TRANSLATE.get(codon, 'X')

_FWD = 'ATGCRYSWKMBDHVN'
_REV = 'TACGYRSWMKVHDBN'
_COMPLEMENT = str.maketrans(_FWD, _REV)
def reverse_complement(sequence):
    """
    Returns the reverse complement of a DNA
    sequence.

    :param sequence: The sequence to translate
    """
    translated = sequence.translate(_COMPLEMENT)
    return ''.join(reversed(translated))

def check_mismatches(seq1, seq2):
    mismatches = 0

    for nuc1, nuc2 in zip(seq1, seq2):

        if nuc1 == nuc2:
            continue

        if nuc2 not in _NUC_SIBLINGS[nuc1] \
            and nuc1 not in _NUC_SIBLINGS[nuc2]:
            
            mismatches += 1

    return mismatches

def chunked_file_reader(file_obj, chunk_size=io.DEFAULT_BUFFER_SIZE):
    while True:
        data = file_obj.read(chunk_size)

        if not data:
            break

        yield data
