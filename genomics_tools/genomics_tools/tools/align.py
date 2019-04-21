###################################################################
#
# All things alignment for the genotyping algorithm
# 
# Author: Milan Patel
# Contact: mpatel5@cdc.gov
# Version 1.0
#
###################################################################

import os
import subprocess as sp
from collections import namedtuple
from environment import (
    log_message,
    log_error,
    check_dir,
    valid_dir
)

BLASTSettings = namedtuple('BLASTSettings', [
    'task', 'identity', 
    'relative_minlen', 'absolute_minlen',
    'include_sequences'
    ])

_platform = os.name

# TODO 
def create_blastdb(fastaflname, dbpath, env):

    toolspath = env.toolsdir

    if not check_dir(toolspath):
        raise RuntimeError('Invalid tools path: {}'.format(
            toolspath))

    # Create the directories if they don't exist
    valid_dir(os.path.dirname(dbpath))

    # This is where the tools should be
    # if that is not the case, either fix it in the filesystem
    # or fix the path below, whichever is easier
    makeblastdb = os.path.join(toolspath, 'all_tools/makeblastdb')

    if _platform == 'nt':
        makeblastdb += '.exe'

    # Check to make sure the tool exists
    if not os.path.exists(makeblastdb):
        raise RuntimeError('Missing ncbi->makeblastdb. Path given:'
            ' {}'.format(toolspath))

    # Make sure that the export of the reference sequences exists
    if not os.path.exists(fastaflname):
        raise RuntimeError('Missing fastafile for blastdb creation'
            ' path given: {}'.format(fastaflname))

    # These are the arguments for blast
    blastdb_args = [

        makeblastdb,
        "-in", fastaflname,
        "-dbtype", "nucl",
        "-out", dbpath
    ]

    # Tell the user that you are doing this
    log_message('BLASTDatabase: Running command: {}'.format(
        ' '.join(blastdb_args)))

    # Call the subprocess
    child = sp.Popen(blastdb_args, stdout=sp.PIPE, stderr=sp.PIPE)

    # Check the output of BLASTMakeDB
    # This will wait for the subprocess to finish before
    # continuing
    stdout, stderr = child.communicate()

    # Check the exit code:
    exit_code = child.returncode
    
    # Log what we got out of the blastn
    for line in stdout.strip().split('\n'):
        log_message(line, extra=1)

    if exit_code:
        log_error(stderr.strip())
        raise RuntimeError('Error making BLASTDatabase')

    log_message('Done creating BLASTDatabase!')

def align_blast_nodb(query, subject, settings, env):
    # There are differences in results between using
    # a formated blastdb, verses just using a
    # subject sequence

    blast_format = [
        '7',
        'qseqid',
        'sseqid',
        'pident',
        'length',
        'mismatch',
        'gapopen',
        'qstart',
        'qend',
        'sstart',
        'send',
        'evalue',
        'bitscore'
    ]

    if settings.include_sequences:
        blast_format.extend(['qseq', 'sseq'])

    # Create the format string
    blast_formatstr = ' '.join(blast_format)

    # Path for the output file
    outputfile = os.path.join(env.localdir, 'blastout.txt')

    # blastn path
    blastn = os.path.join(env.toolsdir, 'all_tools/blastn')

    if _platform == 'nt':
        blastn += '.exe'

    if not os.path.exists(blastn):
        raise RuntimeError('Missing blastn: {}'.format(blastn))

    if not os.path.exists(subject):
        raise RuntimeError('Path to subject sequence does'
            ' not exist {}'.format(subject))

    # BLAST command
    blastn_args = [
        blastn,
       '-task', settings.task,
       '-subject', subject,
       '-query', query,
       # Subtract 1 to inlcude room for the parent that
       # calls this
       '-num_threads', str(max(1, min(4, env.threads-1))),
       '-out', outputfile,
       '-perc_identity', str(int(100.0*settings.identity)),
       '-outfmt',  '{}'.format(blast_formatstr),
       '-max_target_seqs', '1000000',
       '-dust', 'no'
    ]

    log_message('BLASTn running command: {}'.format(
    ' '.join(blastn_args)))

    # Run the blast command
    child = sp.Popen(blastn_args, stdout=sp.PIPE, stderr=sp.PIPE)

    # Check the output of the blastn
    # this will wait until the process has finished
    stdout, stderr = child.communicate()

    # Check the exit code
    exit_code = child.returncode

    # Log what we got out of the blastn
    for line in stdout.strip().split('\n'):
        log_message(line, extra=1)

    if exit_code:
        log_error(stderr.strip())
        raise RuntimeError('Error running BLASTn')

    log_message('Done running BLASTn!')

    # Return the results as a GenotypeResults object
    return GenotypeResults(outputfile, settings, 'blast')

def align_blast(query, blastdb, settings, env):

    # This is the output format that blastn will output
    # 7 is a specific type of predefined header combination
    # providing those headers after the 7 allows you to specify
    # extras that may not be inlcuded in the 7 predefinition
    blast_format = [
        '7',
        'qseqid',
        'sseqid',
        'pident',
        'length',
        'mismatch',
        'gapopen',
        'qstart',
        'qend',
        'sstart',
        'send',
        'evalue',
        'bitscore'
    ]

    # If we want the sequences from the alignment
    # Good for mutation finder and stx subtyper
    if settings.include_sequences:
        blast_format.extend(['qseq', 'sseq'])

    # Create the format string
    blast_formatstr = ' '.join(blast_format)

    # Path for the output file
    outputfile = os.path.join(env.localdir, 'blastout.txt')

    # blastn path
    blastn = os.path.join(env.toolsdir, 'all_tools/blastn')

    if _platform == 'nt':
        blastn += '.exe'

    if not os.path.exists(blastn):
        raise RuntimeError('Missing blastn: {}'.format(blastn))    

    # BLAST command
    blastn_args = [
        blastn,
       '-task', settings.task,
       '-query', query,
       '-db', blastdb,
       # Subtract 1 to inlcude room for the parent that
       # calls this
       '-num_threads', str(max(1, min(4, env.threads-1))),
       '-out', outputfile,
       '-perc_identity', str(int(100.0*settings.identity)),
       '-outfmt',  '{}'.format(blast_formatstr),
       '-max_target_seqs', '1000000',
       '-dust', 'no'
    ]

    log_message('BLASTn running command: {}'.format(
        ' '.join(blastn_args)))

    # Run the blast command
    child = sp.Popen(blastn_args, stdout=sp.PIPE, stderr=sp.PIPE)

    # Check the output of the blastn
    # this will wait until the process has finished
    stdout, stderr = child.communicate()

    # Check the exit code
    exit_code = child.returncode

    # Log what we got out of the blastn
    for line in stdout.strip().split('\n'):
        log_message(line, extra=1)

    if exit_code:
        log_error(stderr.strip())
        raise RuntimeError('Error running BLASTn')

    log_message('Done running BLASTn!')

    # Return the results as a GenotypeResults object
    return GenotypeResults(outputfile, settings, 'blast')

class GenotypeHit(object):

    def __init__(self):
        self.reference_id = ''
        self.query_id = ''
        self._query_start = 0
        self._query_stop = 0
        self._reference_start = 0
        self._reference_stop = 0
        self.forward = True
        self.identity = 0.
        self.absolute_len = 0
        self.relative_len = 0
        self.reference_len = 0
        self.num_mismatches = 0
        self.num_gap_opens = 0
        self.query_seq = None
        self.reference_seq = None
        self.bitscore = 0
        self.evalue = 0
        self.full_match = False

    # ************* IMPORTANT *************
    # 
    # The following methods are used to cope
    # with using 1 indexed BLAST indices.
    # If you are going to use the setters
    # make sure you set it with the 0 indexed
    # value, which is reasonable since you are getting
    # the 0 indexed values from the below
    # property methods to begin with.
    # We did this so I can have good feels about this
    # datastructure (thanks).
    # 

    @property
    def query_start(self):
        return self._query_start - 1

    @query_start.setter
    def query_start(self, val):

        self._query_start = val + 1

    @property
    def query_stop(self):
        return self._query_stop - 1

    @query_stop.setter
    def query_stop(self, val):
        self._query_stop = val + 1

    @property
    def reference_start(self):
        return self._reference_start - 1

    @reference_start.setter
    def reference_start(self, val):
        self._reference_start = val + 1

    @property
    def reference_stop(self):
        return self._reference_stop - 1

    @reference_stop.setter
    def reference_stop(self, val):
        self._reference_stop = val + 1

    @staticmethod
    def from_blast(line):
        # query_id  ref_id|len  iden    alignmentlen    mismatches
        #   gapopens    qstart  qstop   refstart    refstop evalue
        #   bitscore    queryseq    refseq 

        # Create the object
        hit = GenotypeHit()

        # Split the line
        parts = line.split('\t')

        # Parse the results
        hit.query_id = parts[0]
        hit.reference_id, ref_len = parts[1].split('|')
        hit.identity = float(parts[2])
        hit.absolute_len = int(parts[3])
        hit.num_mismatches = int(parts[4])
        hit.num_gap_opens = int(parts[5])

        # NOTE:
        # The alignments come in 1-indexed
        hit._query_start = int(parts[6])
        hit._query_stop = int(parts[7])
        hit._reference_start = int(parts[8])
        hit._reference_stop = int(parts[9])

        hit.evalue = float(parts[10])
        hit.bitscore = float(parts[11])

        # Not all of the genotyping needs to 
        # have the sequences returned
        if len(parts) == 14:
            hit.query_seq = parts[12]
            hit.reference_seq = parts[13]

        # Check to see if the len was provided
        if ref_len:
            hit.reference_len = int(ref_len)
        else:
            hit.reference_len = 0

        # Check to see if the sequence has been reversed
        # If for whatever reason the hit is one nucleotide
        # long, this will be false
        hit.forward = hit._reference_start < hit._reference_stop

        if not hit.forward:
            hit._reference_start, hit._reference_stop = \
                hit._reference_stop, hit._reference_start

        # Change the percentages into decimals,
        # because why not
        hit.identity = hit.identity / 100.

        # Make sure to avoid a divide by zero error
        if hit.reference_len:
            hit.relative_len = float(hit.absolute_len) / float(hit.reference_len)

        if not hit.num_gap_opens and hit.identity == 1.0 and \
            hit.relative_len == 1.0:
            hit.full_match = True

        return hit

    @staticmethod
    def from_mummer(line):
        # PLACEHOLDER FOR IMPLEMENTATION
        raise NotImplementedError('Mummer analysis not yet implemented')

class GenotypeResults(object):

    hit_handlers = {
        'blast' : GenotypeHit.from_blast,
        'mummer': GenotypeHit.from_mummer
    }

    def __init__(self, filename, settings, aligner='blast'):
        self._hits = []
        self._aligner = aligner
        self._settings = settings

        self.load_hits(filename, settings, aligner)

    def read_file(self, flobj):

        found_start = False

        for line in flobj:

            line = line.strip()

            # Empty line
            if not len(line):
                continue

            # Comment line
            if line[0] == '#' or line[0] == '=':
                found_start = True
                continue

            # If you haven't found the start of the file
            if not found_start:
                continue

            yield line

    def load_hits(self, filename, settings, aligner):

        if aligner in GenotypeResults.hit_handlers:
            handler = GenotypeResults.hit_handlers[aligner]
        
        else:
            raise RuntimeError('Requested non-existant aligner')

        # Make sure the settings provided has the correct settings in it
        required_settings = ['relative_minlen', 'absolute_minlen']

        if not all(hasattr(settings, x) for x in required_settings):
            raise RuntimeError('Missing settings for handling'
                ' genotyping hits')

        # If we are given a file path rather than a file obj
        if isinstance(filename, basestring):

            if os.path.exists(filename):

                with open(filename, 'r') as f:

                    for line in self.read_file(f):

                        self._hits.append(handler(line))

            else:

                raise RuntimeError('Provide file path is not a real'
                    ' path for alignment results parsing')

        else:
            # We were given a file handle
            for line in self.read_file(filename):

                self._hits.append(handler(line))

    @property
    def hits(self):
        return self._hits
