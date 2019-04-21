###################################################################
#
# Tools for the genotyping algorithm
# 
# Author: Milan Patel
# Contact: mpatel5@cdc.gov
# Version 1.0
#
###################################################################

import base64
import os
import sys
import gzip
import json
import shutil
import binascii
import subprocess as sp
from string import maketrans
from itertools import izip, combinations
from collections import OrderedDict, namedtuple

from Bio import SeqIO

from .environment import (
    full_path,
    log_message,
    log_warning,
    log_exception,
    log_error
)

_FASTAEXTS = ['.fna', '.fasta', '.fsa']
_GENBANKEXTS = ['.gb', '.gbk']
_FWD = 'ATGCRYSWKMBDHVN'
_REV = 'TACGYRSWMKVHDBN'
_COMPLEMENT = maketrans(_FWD, _REV)
_GZIP_START = b'1f8b'

# Get the codons for any particular amino acid
_AA_BACK_TRANSLATE = {
    'A' : ['GCT', 'GCC', 'GCA', 'GCG'],
    'L' : ['TTA', 'TTG', 'CTT', 'CTC', 'CTA', 'CTG'],
    'R' : ['CGT', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG'],   
    'K' : ['AAA', 'AAG'],
    'N' : ['AAT', 'AAC'],
    'M' : ['ATG'],
    'D' : ['GAT', 'GAC'],
    'F' : ['TTT', 'TTC'],
    'C' : ['TGT', 'TGC'],
    'P' : ['CCT', 'CCC', 'CCA', 'CCG'],
    'Q' : ['CAA', 'CAG'],
    'S' : ['TCT', 'TCC', 'TCA', 'TCG', 'AGT', 'AGC'],
    'E' : ['GAA', 'GAG'],
    'T' : ['ACT', 'ACC', 'ACA', 'ACG'],
    'G' : ['GGT', 'GGC', 'GGA', 'GGG'],
    'W' : ['TGG'],
    'H' : ['CAT', 'CAC'],
    'Y' : ['TAT', 'TAC'],
    'I' : ['ATT', 'ATC', 'ATA'],
    'V' : ['GTT', 'GTC', 'GTA', 'GTG']
}

# Get the amino acid for any particular codon
_AA_TRANSLATE = {}

# Create the dictionary on the fly
for amino, codons in _AA_BACK_TRANSLATE.iteritems():
    for codon in codons:
        _AA_TRANSLATE[codon] = amino

# All potential sequence based nucleotides
# and their parings
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
for nuc, pairs in _NUC_SIBLINGS.iteritems():
    _NUC_SIBLINGS[nuc] = frozenset(pairs)
    _SET_TO_NON_ACTG.__setitem__(_NUC_SIBLINGS.get(nuc), nuc)

def get_non_iupac(f_set):
    return _SET_TO_NON_ACTG.get(f_set, None)

def add_cmdline_args(executable_name, cmd_args, args):
    """
    Modifies a list used for commandline args in place while
    doing some basic validations
    """

    if not isinstance(cmd_args, list):
        raise TypeError('Commandline args destination '
                        'must be of type list, got'
                        ' {} instead'.format(type(cmd_args))
                        )

    if not isinstance(args, (tuple, list)):
        raise TypeError('Args to add must be in list or tuple format'
                        ' got {} instead'.format(type(args))
                        )
    for arg in args:

        if not isinstance(arg, basestring) or not \
            arg.startswith('-'):

            log_warning('Invalid argument provided for {}:'
                        ' {}, skipping..'.format(executable_name, arg))
            continue

        cmd_args.append(arg)

def add_cmdline_kwargs(executable_name, cmd_args, kwargs):
    """
    Modifies a list used for commandline args in place while
    doing some basic validations
    """

    if not isinstance(cmd_args, list):
        raise TypeError('Commandline args destination '
                        'must be of type list, got'
                        ' {} instead'.format(type(cmd_args))
                        )

    if not isinstance(kwargs, dict):
        raise TypeError('Key word args to add must be in dict format'
                        ' got {} instead'.format(type(kwargs))
                        )

    for arg, value in kwargs.iteritems():
        if not isinstance(arg, basestring) or not \
            arg.startswith('-'):

            log_warning('Invalid argument provided for {}:'
                        ' {} -. {}, skipping..'.format(
                                                    executable_name, 
                                                    arg,
                                                    value
                                                )
                        )

            continue

        cmd_args.extend([arg, str(value)])

def popen(args, stdout=None, stderr=None, cwd=None, shell=False):

    if not isinstance(args, list) or not args:
        raise RuntimeError('Provided arguments must be of type list')

    if not stderr:
        stderr = sp.PIPE

    if not stdout:
        stdout = sp.PIPE

    if os.path.isfile(args[0]):
        args[0] = full_path(args[0])

    child = sp.Popen(args, stdout=stdout, stderr=stderr, cwd=cwd, shell=shell)

    out, err = child.communicate()

    return child.returncode, out, err

def check_gzipped(file_path):

    # Read the first two bytes to see if it is the gzip
    # header
    with open(file_path, 'rb') as f:
        head = f.read(2)

    return binascii.hexlify(head) == _GZIP_START

def check_b64encoded(string):
    try:
        base64.decodestring(string)
    
    except binascii.Error:
        return False

    else:
        return True

def get_all_file_exts(path):
    root, ext = os.path.splitext(path)

    if not ext:
        return root, []

    all_exts = [ext]

    real_root, exts = get_all_file_exts(root)
    exts += all_exts

    return real_root, exts

def process_seq_file(file_path, load=True):

    log_message('Checking provided '
                'sequence file: {}'.format(file_path))

    if not os.path.exists(file_path):
        raise RuntimeError('No file exists for path: {}'
                            ''.format(str(path)))

    out_path = file_path

    # First check to see if the data is gzipped:
    if check_gzipped(file_path) or file_path.endswith('.gz'):

        log_message('Found gzipped file!', extra=1)

        root, exts = get_all_file_exts(file_path)

        if exts and exts[0] != '.gz':
            out_path = root + exts[0]
        else:
            out_path = root + '.fasta'

        log_message('Unzipping...', extra=1)
        
        try:
            unzip_file(file_path, out_path)

        except:
            log_error('Unable to unzip query file!')
            raise

    if load:
        return out_path, parse_fasta(out_path)

    else:
        return out_path, None

def process_read_files(reads, load=False):

    if not isinstance(reads, list):
        raise TypeError('Read files must be provided as a'
            ' list of paths')

    to_return = []
    for read in reads:
        to_return.append(process_seq_file(read, load=load))

    return to_return

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

def is_fasta(path):
    # Check if the file is a fasta format
    if not os.path.exists(path):
        return False

    file_extension = os.path.splitext(path)[1]

    return file_extension in _FASTAEXTS

def is_genbank(path):
    # Check if the file is a genbank format
    if not os.path.exists(path):
        return False

    file_extension = os.path.splitext(path)[1]

    return file_extension in _GENBANKEXTS

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

    for nuc1, nuc2 in izip(seq1, seq2):

        if nuc1 == nuc2:
            continue

        if nuc2 not in _NUC_SIBLINGS[nuc1] \
            and nuc1 not in _NUC_SIBLINGS[nuc2]:
            
            mismatches += 1

    return mismatches

# TODO rewrite to fit with code-base
def parse_paired_files(readFiles):
    """
    Written by not me as you can see
    """
    
    tokens = [['_1.', '_R1_', '_R1'], ['_2.', '_R2_', '_R2']]

    singleFiles = []
    pairedFiles = []
    pairedFilesDict = {}
    for fl in readFiles:
        isFirst = any(token in fl for token in tokens[0])
        isSecond = any(token in fl for token in tokens[1])
        if not isFirst and not isSecond:
            singleFiles.append(fl)
        else:
            if isFirst:
                if fl not in pairedFilesDict:
                    pairedFilesDict[fl] = [fl, '']
                else:
                    pairedFilesDict[fl][0] = fl
            if isSecond:
                firstFl = fl
                for t1, t2 in izip(tokens[0], tokens[1]):
                    firstFl = firstFl.replace(t2, t1)
                if firstFl not in pairedFilesDict:
                    pairedFilesDict[firstFl] = ['', fl]
                else:
                    pairedFilesDict[firstFl][1] = fl
    for flPair in pairedFilesDict.itervalues():
        if flPair[0] == '' and flPair[1] != '':
            singleFiles.append(flPair[1])
            continue
        if flPair[0] != '' and flPair[1] == '':
            singleFiles.append(flPair[0])
            continue
        pairedFiles.append(flPair)

    return singleFiles, pairedFiles

def counter(l):
    c = {}

    for i in l:

        if i in c:
            c[i] += 1
        else:
            c[i] = 1

    return c

# TODO rewrite to fit with code-base
def combine_fastas(inputFastaList, outputFasta):
    fastaRecords = list()

    for fastaFL in inputFastaList:
        handle=None
        if os.path.splitext(fastaFL)[1] == ".gz":
            handle = gzip.open(fastaFL, 'rU')
        else:
            handle = open(fastaFL, 'rU')
        for record in SeqIO.parse(handle, "fasta"):
            fastaRecords.append(record)
        handle.close()

    fastaIDs = [record.id for record in fastaRecords]

    assert len(set(fastaIDs)) == len(fastaIDs),"There are some non-unique reference names!"

    SeqIO.write(fastaRecords, outputFasta, "fasta")

    return outputFasta
