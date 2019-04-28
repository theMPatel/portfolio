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

_FWD = 'ATGCRYSWKMBDHVN'
_REV = 'TACGYRSWMKVHDBN'
_COMPLEMENT = str.maketrans(_FWD, _REV)
_GZIP_START = b'1f8b'

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

def unzip_file(fl_in, fl_out):
    # Assumes that the file is a gzipped
    # file

    # Open the in file
    if not hasattr(fl_in, 'read'):
        fl_in = gzip.open(fl_in, 'rb')

    # Open the outfile
    if not hasattr(fl_out, 'write'):
        fl_out = open(fl_out, 'wb')

    # If the files are already open
    # it will jump to this step
    shutil.copyfileobj(fl_in, fl_out)

    # Close the files to be nice
    fl_in.close()
    fl_out.close()

    # Make sure that the file is in proper dos2unix format:
    dos2unix = ['dos2unix', None]

    # If a file handle was passed, extract the file path
    if hasattr(fl_out, 'write'):
        dos2unix[1] = fl_out.name

    else:
        dos2unix[1] = fl_out

    success = popen(dos2unix)

    if success[0]:
        raise RuntimeError('Error converting file to dos: {}\n{}\n'.format(
            success[1], success[2]))

_FASTAEXTS = set(['.fna', '.fasta', '.fsa'])
def is_fasta(path):
    # Check if the file is a fasta format
    if not path:
        return False

    root, exts = get_all_file_exts(path)

    return bool(set(exts) & _FASTAEXTS)

_GENBANKEXTS = set(['.gb', '.gbk'])
def is_genbank(path):
    # Check if the file is a genbank format
    if not path:
        return False

    root, exts = get_all_file_exts(path)

    return bool(set(exts) & _GENBANKEXTS)

def fasta_iterator(flname):
    # Files look like this:
    #
    # >fasta_id
    # ACTGACTGACTGACTGACTGACTGACTG\n
    # ACTGACTGACTGACTGACTGACTGACTG\n
    # ACTGACTGACTGACTGACTGACTGACTG\n
    # 
    # Parse accordingly

    with open(flname, 'r') as f:
        # f is an interable

        # Stores the sequences
        sequence_parts = []

        # Key for a sequence
        key = ''

        for line in f:

            # Get rid of the newline
            line = line.strip()

            # Empty line
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
    
    if is_fasta(flname):

        fasta_sequences = {}

        if rename:

            for i, (name, sequence) in enumerate(fasta_iterator(flname),1):

                new_name = 'contig_' + str(i)

                fasta_sequences[new_name] = sequence

            return fasta_sequences

        else:

            fasta_sequences.update(fasta_iterator(flname))

            return fasta_sequences

    else:
        raise RuntimeError('Requested file path is not of type fasta')

def codon_translation(codon):
    # Translates the codon that was requested
    if len(codon) != 3:
        raise ValueError('Cannot translate a codon: {} of '
            'non-standard size: {}'.format(
                codon, len(codon)))

    return _AA_TRANSLATE.get(codon, 'X')

def reverse_complement(sequence):
    # Reverses the sequence at hand
    translated = sequence.translate(_COMPLEMENT)
    return ''.join(reversed(translated))

def binary_search(arr, value, lo, hi):
    # Use this function to do a binary search of your list.
    # This works only on a list of non-duplicates, which is to
    # say that if you have duplicates, you will get the index
    # for only one of them. It should be easy enough for you
    # to search around these values *if* you need to.
    # The purpose of this function is to find gap indices
    # for which there can be no duplicates so we should be fine.
    # Returns:
    #   -Either the index of the search value ***OR***
    #   -The index of the largest element
    #    smaller than the requested search value
    #   - If nothing was found then -1 is returned

    if lo > hi:
        return hi

    mid = (hi+lo) // 2

    if value == arr[mid]:
        return mid

    elif value < arr[mid]:
        return binary_search(arr, value, lo, mid-1)

    else:
        return binary_search(arr, value, mid+1, hi)

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
