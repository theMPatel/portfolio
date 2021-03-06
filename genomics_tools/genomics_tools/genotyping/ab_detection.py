###################################################################
#
# Tools for assembly based genotyping prediction
# 
# This tool is a sample and distillation of the real application
# hosted at: https://github.com/theMPatel/functional_genomics_tools
#
# Author: Milan Patel
# Contact: https://github.com/theMPatel
# Version 1.0
#
###################################################################

from tools.environment import (
    log_message,log_error,
    log_progress,valid_dir
)

from tools.align import (
    BLASTSettings, create_blastdb,
    align_blast_nodb, GenotypeHit
)

from tools.tools import (
    reverse_complement, codon_translation
)

import os
from itertools import combinations, product
from tools.fancy_tools import Disjointset, binary_search
from collections import defaultdict, namedtuple

GenotypeRegion = namedtuple('GenotypeRegion', ['coverage', 'identity', 'locations'])

def mutation_detector(sequence_database, query_path, percent_identity,
    min_relative_coverage, env):
    """
    The primary dispatcher and external interface for the mutation 
    detection pipeline.

    :param sequence_database: The reference sequences to use.
    :param query_path: The path to the query_file to search mutations in.
    :param percent_identity: The minimum percent identity for alignment
        matches.
    :param min_relative_coverage: The minimum coverage in alignment for
        a gene.
    :param env: The env object to retrieve information from
    """
    
    log_message('Exporting references...')
    reference_dir = os.path.join(env.tempdir, 'blastdb')
    valid_dir(reference_dir)
    reference_path = os.path.join(reference_dir, 'references.fasta')
    sequence_database.export_sequences(reference_path)
    log_message('Successfully exported reference database...')

    # Create the path to the blast database
    blast_db_path = os.path.join(env.tempdir, 'blastdb', 'db.fasta')
    log_message('Creating blast database...')

    # Log that we were successful
    log_message('Successfully created blast database!')

    # Return the sequences so that we can search for
    # point mutations
    blast_settings = BLASTSettings(
        task = 'blastn',
        identity = percent_identity,
        relative_minlen = 0,
        absolute_minlen = 0,
        include_sequences = True
        )

    log_message('BLASTing query genome against reference database')
    results = align_blast_nodb(
        query_path,
        reference_path,
        blast_settings,
        env
    )

    log_message('Successfully BLASTed query genome against reference database')
    log_message('Searching for mutations...')
    interpretations = find_mutations(
        sequence_database,
        results,
        min_relative_coverage)

    log_message('Retained {} gene regions after gene analysis'.format(
        len(interpretations)))

    return interpretations

def find_mutations(sequence_database, results, min_relative_coverage):
    """
    Primary function that will search for mutations in BLAST hits

    :param sequence_database: The database of reference sequences
    :param results: The BLAST hit results
    :param min_relative_coverage: The minimum coverage in alignment for
        a gene.
    """

    regions = defaultdict(list)

    for hit in results.hits:
        regions[hit.reference_id].append(hit)

    log_message('Found {} potential regions of interest'.format(
        str(len(regions))))

    # Store the found resistance:
    mutation_results = defaultdict(list)
    for reference, hits in regions.items():
        targets = sequence_database.targets[reference]

        for hit in hits:
            if hit.relative_len < min_relative_coverage:
                continue

            # Store the gap positions for each of the sequences.
            # They will naturally be sorted, we can do binary search
            # for figuring out if we need to offset the string indices
            gaps = hit.num_gap_opens
            ref_gaps = []
            query_gaps = []

            if gaps:
                # Deletions
                query_gaps = [i for i, s in enumerate(hit.query_seq) if s == '-']
                # Insertions
                ref_gaps = [i for i, s in enumerate(hit.reference_seq) if s == '-']

            for target in targets:
                if target.coding_gene:
                    results = validate_coding_gene(hit, target, ref_gaps,
                                                    query_gaps)
                    if results:
                        mutation_results[hit.reference_id].append(results)

                else:
                    results = validate_noncoding_gene(hit, target, ref_gaps,
                                                        query_gaps)

                    if results:
                        mutation_results[hit.reference_id].append(results)

    return mutation_results

def validate_coding_gene(hit, target, ref_gaps, query_gaps, **kwargs):
    """
    Cross references the hit against the target including any insertions
    or deletions in the query and returns back a prediction.

    :param hit: The hit result for a BLAST hit
    :param target: The database target that maps to the hit
    :param ref_gaps: A list of gap indices in the reference sequence
    :param query_gaps: A list of gap indices in the query sequence
    """

    default_return = {}
    ref_start = hit.reference_start
    ref_stop = hit.reference_stop
    query_start = hit.query_start
    query_stop = hit.query_stop
    hit_ref_seq = hit.reference_seq
    hit_query_seq = hit.query_seq

    codon_position = target.codon_position
    reference_codons = target.reference_codon
    reference_aas = target.reference_aa
    resistance_aas = target.resistance_aa
    resistances = target.resistance

    # The sequences come back from blast as the reverse
    # complement of the reference, thus we need to reverse
    # them back to get a sequence in the same orientation
    # as the reference sequence ***IF*** the alignment is a
    # reverse alignment
    if not hit.forward:
        hit_ref_seq = reverse_complement(hit_ref_seq)
        hit_query_seq = reverse_complement(hit_query_seq)

    # The positional information for the coding genes
    # is stored as the position of the codon
    # i.e. codon 415 is nucleotide position 415*3
    # Get the positional information for this target
    # These are the (1-indexed) values which we
    # transform into the (0-indexed) position.
    # You could leave the end position without
    # subtracting 1 to make list slicing easier
    # Instead, I've left this here for the next
    # person so they can understand what is really
    # happening
    start = (codon_position * 3) - 3
    end = (codon_position * 3 ) - 1

    # If the codon does not exist within the
    # range of the hit, obviously it can't
    # exist here. NOTE: despite any
    # insertions into the reference sequence
    # the start and stop will be the absolute
    # start and stop within the reference sequence
    if start < ref_start or end > ref_stop:
        return {}

    # Create the indices for the returned
    # query string
    string_start = start - ref_start
    string_end = end - ref_start

    # Determine how many gaps there are in the
    # returned reference sequence (insertions)
    max_ins = binary_search(
        ref_gaps, string_start, 0, len(ref_gaps)-1)

    ins_offset = max_ins + 1

    # The case where the mutation is the result of an
    # insertion
    if '-' not in reference_codons:
        while hit_ref_seq[string_start+ins_offset] == '-' and \
            string_end+ins_offset+1 < len(hit_ref_seq):
            
            ins_offset += 1

    if len(hit_ref_seq) - 1 < string_end + ins_offset + 1:
        return default_return

    # Check to see that the reference sequence
    # that is returned from blast
    # does not have any gaps that would
    # cause a frame-shift
    hit_ref_codon = hit_ref_seq[
        string_start + ins_offset : string_end + ins_offset + 1]

    # Sanity check
    assert len(hit_ref_codon) == 3

    if hit_ref_codon not in reference_codons:
        # If this happens, it should be because
        # there are insertions into your query
        # genome that caused there to be gaps
        # in your reference sequence within the
        # codon itself
        return default_return

    # Determine how many gaps there are in the returned
    # query (deletions)
    max_dels = binary_search(
        query_gaps, string_start, 0, len(query_gaps)-1)

    # binary_search will return -1 if there are no dels
    del_offset = max_dels + 1

    if '-' not in resistance_aas:
        while hit_query_seq[string_start + del_offset] == '-' and \
            string_end + del_offset + 1 < len(hit_query_seq):

            del_offset += 1

    if len(hit_query_seq) - 1 < string_end + del_offset + 1:
        return default_return

    hit_query_codon = hit_query_seq[
        string_start + del_offset : string_end + del_offset + 1]

    assert len(hit_query_codon) == 3

    ref_slice = create_sequence_slice(hit_ref_seq, string_start, 
                                        ins_offset, 10)

    query_slice = create_sequence_slice(hit_query_seq, string_start,
                                        del_offset, 10)

    # Get the translation
    query_translation = codon_translation(hit_query_codon)

    if query_translation in resistance_aas:

        return {
            'locus' : target.gene_id,
            'identity' : hit.identity,
            'contig_id' : hit.query_id,
            'query_codon': hit_query_codon,
            'query_aa' : query_translation,
            'position': codon_position,
            'aa_mutation' : '{}->{}'.format(
                reference_aas[0], query_translation),
            'resistance' : resistances,
            'hit': hit,
            'iscoding' : target.coding_gene,
            'reference': hit_ref_seq[ref_slice],
            'query' : hit_query_seq[query_slice]
        }

    return default_return

def validate_noncoding_gene(hit, target, ref_gaps, query_gaps, **kwargs):
    """
    Cross references the hit against the target including any insertions
    or deletions in the query and returns back a prediction.

    :param hit: The hit result for a BLAST hit
    :param target: The database target that maps to the hit
    :param ref_gaps: A list of gap indices in the reference sequence
    :param query_gaps: A list of gap indices in the query sequence
    """

    default_return = {}
    ref_start = hit.reference_start
    ref_stop = hit.reference_stop
    query_start = hit.query_start
    query_stop = hit.query_stop
    hit_ref_seq = hit.reference_seq
    hit_query_seq = hit.query_seq

    # The RNA genes and ampC promoter mutation positions are stored
    # as actual positions of the point mutations since these do not 
    # actually code for genes.
    codon_position  = target.codon_position
    reference_codons = target.reference_codon
    reference_aas = target.reference_aa
    resistance_aas = target.resistance_aa
    resistances = target.resistance

    # The ampC promoter mutations are stored as
    # negative indexes from the back of the
    # promoter region. The actual value of the
    # position is going to be 53 + (-index)
    # The promoter is a 54 bp promoter, however
    # It's the first 53 base pairs that make up the
    # negative indices. Here is an example of negative 
    # indices from the actual ampC sequence:
    # 
    #   -11 -10 -9  -8  -7  -6  -5  -4  -3  -2  -1   1
    #    C   A   A   T   C   T   A   A   C   G   C   A
    # 
    if codon_position < 0:
        start = codon_position + 53

    # If it's not negative, then its a 16s or 23s
    # RNA gene
    else:
        start = codon_position - 1

    # If the nucleotide does not exist within the
    # range of the hit, obviously it can't
    # exist here. NOTE: despite any
    # insertions into the reference sequence
    # the start and stop will be the absolute
    # start and stop within the reference sequence
    if not ref_start <= start <= ref_stop:
        return default_return

    # Create the indices for the returned
    # query string
    string_start = start - ref_start

    # Determine how many gaps there are in the
    # returned reference sequence
    max_ins = binary_search(
        ref_gaps, string_start, 0, len(ref_gaps)-1)

    # If max_ins is -1 then it means there are no
    # gaps at indices smaller than the index of the
    # first nucleotide in the codon
    ins_offset = max_ins + 1

    if '-' not in reference_codons:
        while hit_ref_seq[string_start+ins_offset] == '-' and \
            string_start+ins_offset+1 < len(hit_ref_seq):
            
            ins_offset += 1

    if len(hit_ref_seq) - 1 < string_start + ins_offset:
        return default_return

    hit_ref_nucleotide = hit_ref_seq[string_start+ins_offset]

    if hit_ref_nucleotide not in reference_codons:
        # This should only occur in the case that
        # your query sequence has insertions which
        # caused the existence of gaps in your 
        # reference sequence
        return default_return

    max_dels = binary_search(
        query_gaps, string_start, 0, len(query_gaps)-1)

    del_offset = max_dels + 1

    if '-' not in resistance_aas:
        while hit_query_seq[string_start + del_offset] == '-' and \
            string_start + del_offset + 1 < len(hit_query_seq):

            del_offset += 1

    if len(hit_query_seq) - 1 < string_start + del_offset:
        return default_return

    # Get the query string nucleotide
    hit_query_nucleotide = hit_query_seq[string_start+del_offset]

    ref_slice = create_sequence_slice(hit_ref_seq, string_start, 
                                        ins_offset, 10)

    query_slice = create_sequence_slice(hit_query_seq, string_start,
                                        del_offset, 10)

    # Figure out if the hit is the same as a resistance
    # nuc that is known
    if hit_query_nucleotide in resistance_aas:

        return {
            'locus' : target.gene_id,
            'identity' : hit.identity,
            'contig_id' : hit.query_id,
            'query_codon': hit_query_nucleotide,
            'position': codon_position,
            'reference_nuc' : hit_ref_nucleotide,
            'resistance' : resistances,
            'aa_mutation' : '',
            'hit': hit,
            'iscoding' : target.coding_gene,
            'reference': hit_ref_seq[ref_slice],
            'query' : hit_query_seq[query_slice]
        }

    return default_return

def create_sequence_slice(seq, start, offset, size):
    """
    Creates and returns an appropriate slice object that we
    can use to create a localized view around a mutation
    in the query and reference sequence

    :param seq: The actual sequence
    :param start: The start of the mutation
    :param: offset: Any offset that might be needed
    :param size: The size of the slice
    :returns: A constructed slice object
    :rtype: `slice`
    """

    start_pos = max(0, start+offset-(size//2))
    end_pos = min(len(seq), start_pos+size+1)

    return slice(start_pos, end_pos, 1)

def encompassed(hit1, hit2, min_coverage):
    """
    Calculate whether one hit captures another:
    
            max(start)     min(end)
    hit    .___________________.
    hit  ._________________________.
    
    hit1 and hit2 come in as lists
    in case one is a fragment
    
    :param hit1: The first hit list to check
    :param hit2: The second hit list to check
    :param min_coverage: The minimum coverage needed between
        the two hits
    """

    hit1_length = sum(
        hit.query_stop - hit.query_start + 1 for hit in hit1)

    hit2_length = sum(
        hit.query_stop - hit.query_start + 1 for hit in hit2)

    total_overlap = 0

    # Get all the combinations of the two hits
    # will only iterate once unless we have
    # fragments
    for combo in product(hit1, hit2):

        # If they are not on the same contig, obviously
        # no overlap is possible
        if combo[0].query_id != combo[1].query_id:
            continue

        # Calculate the overlap
        overlap = min(combo[0].query_stop, combo[1].query_stop) - \
            max(combo[0].query_start, combo[1].query_start) + 1

        # If the calculated overlap is negative
        overlap = max(0, overlap)

        # Add it to the total_overlap
        total_overlap += overlap

    return total_overlap >= min_coverage * min(
        hit1_length, hit2_length)

def eliminate_overlap(regions, min_merge_overlap):
    """
    This function uses a disjoinset set to merge all hits
    that the aligner resolved to basically the same area.
    Usually this happens when there is more than more one
    reasonable way to align a reference to a contig

    :param regions: A mapping of references to their hits in
        the query
    :param min_merge_overlap: The threshold at which two hits 
        that cover each other are considered the same
    """

    # Create a flat list of all of the hits
    regions = [hit for region in regions.itervalues() for \
        hit in region.predicted]

    dset = Disjointset(len(regions))

    for i, j in combinations(range(len(regions)), 2):
        hit1 = regions[i].locations
        hit2 = regions[j].locations

        if encompassed(hit1, hit2, min_merge_overlap):
            dset.merge(i,j)

    best_hits = defaultdict(set)
    for i in range(len(regions)):

        parent = dset.get_parent(i)

        # Add it to the set of children
        #   NOTE:
        # A parent is a child of itself,
        # so the children list will include all
        # of the parent's children plus itself
        best_hits[parent].add(i)

    # For each of the best parents, see if the children are better
    to_remove = []
    for parent, children in best_hits.items():
        hits_here = list(children)

        # Sort the regions based on identity:
        hits_here.sort(key=lambda x: -regions[x].identity)

        # Get rid of all the worst ones:
        to_remove.extend(hits_here[1:])

    # Sort the indices
    to_remove.sort()
    for index in reversed(to_remove):
        del regions[index]

    accepted = defaultdict(list)
    for region in regions:
        references_here = set()

        for hit in region.locations:
            references_here.add(hit.reference_id)

        if len(references_here) > 1:
            raise RuntimeError('This should never happen: more than one reference '
                ' associated with a region!')

        accepted[references_here.pop()].append(region)

    return accepted

class Genotype(object):

    def __init__(self, seq_id, reference_len, hits):
        self._reference_id = seq_id
        self._reference_len = reference_len
        self._hits = hits
        self._predicted = []
        self._coverage = None
        self._identity = None

    def validate(self, percent_identity, min_relative_coverage, 
        contig_sizes, search_fragments):
        # Validate all of the hits to make sure they are indeed good hits

        def at_edge(hit):
            """
            The reason for calculating whether or not a hit
            is at the edge is in case you have references
            that are split over two contigs:
            
                  contig1             contig2
                           <25bp> <25bp>
            --------------(------|------)-----------------
                              ^      ^
                              |      |
               ends within here      starts within here
            
            """

            # Get the length of the contig
            contig_end = contig_sizes[hit.query_id]-1
            return hit.query_start <= 25 or contig_end-hit.query_stop <= 25

        best_hits = []
        to_check = []

        for hit in self._hits:

            # Get the start stop of the hit
            stop = hit.query_stop
            start = hit.query_start

            # Calculate the relative coverage
            coverage = float(stop - start + 1) / float(self._reference_len)

            # If there are insertions, then start stop of the query
            # can be larger than the length of the reference sequence
            coverage = min(1.0, coverage)

            # If we have good coverage and the identity is solid,
            # lets keep it
            if coverage >= min_relative_coverage and \
                hit.identity >= percent_identity:

                # Keep track of the hit that passed
                geno_region = GenotypeRegion(
                    coverage = coverage,
                    identity = hit.identity,
                    locations = [hit]
                )

                best_hits.append(geno_region)

            # We might have a fragmented sequence so check it out
            else:
                to_check.append(hit)

        if search_fragments:
            # Get a list of all of the edge hits
            to_check = list(filter(at_edge, to_check))

            # Find the best pair of hits
            # since references broken by contigs should be pairs
            for hits in combinations(to_check, 2):

                # Create a mask of the reference to see what
                # amount of the reference is covered
                mask = [0.]*self._reference_len
                
                # For each of the hits set the values
                for hit in hits:

                    # Create the slice
                    mask_slice = slice(
                        hit.reference_start,
                        # This needs to be plus 1
                        # since list slicing is 
                        # exclusive e.g. [1, 10)
                        hit.reference_stop+1,
                        None
                        )

                    # Size of the slice
                    size = hit.reference_stop - hit.reference_start + 1

                    # Set the mask slice
                    mask[mask_slice] = [hit.identity] * size


                # Get the coverage of the two fragments
                coverage = float(len(filter(None, mask))) / float(
                    self._reference_len)

                # Get the real coverage, in case the calculated
                # coverage is greater than 1
                coverage = min(1.0, coverage) 

                # If the coverage of the paired hits is good
                # then check the identity
                if coverage >= min_relative_coverage:

                    # The identity of the fragments is the best
                    # identities of both of the fragments
                    # divided by its relative coverage
                    # which comes out to be just the number of
                    # items in the mask greater than 0
                    identity = sum(mask) / coverage * float(
                        self._reference_len)

                    if identity >= percent_identity:

                        # Keep track of the hit that passed
                        geno_region = GenotypeRegion(
                            coverage = coverage,
                            identity = identity,
                            locations = hits
                        )

                        best_hits.append(geno_region)

        # If there are no best hits, then this
        # genotype is likely to not be present in 
        # our sample
        if not len(best_hits):
            return False

        # Sort the best hits
        best_hits.sort(key=lambda x: float(x.coverage)*float(x.identity))

        # Keep the best one
        self._predicted.append(best_hits.pop())

        # Check to make sure that we've removed all hits
        # for this reference that are contained
        # within the same locations
        while len(best_hits):

            to_check = best_hits.pop()

            for accepted in self._predicted:

                if not encompassed(
                    accepted.locations,
                    to_check.locations,
                    min_relative_coverage
                ):

                    self._predicted.append(to_check)

        # For this genotype save the coverage and 
        # identity
        self._coverage = self._predicted[0].coverage
        self._identity = self._predicted[0].identity
        return True

    @staticmethod
    def find_regions(results, contig_sizes, sequence_database):
        # Find all the regions in the results
        hit_regions = defaultdict(list)

        # Group together all the hits that are for the same
        # reference
        for hit in results.hits:
            hit_regions[hit.reference_id].append(hit)

        # Dictionary for all of the genotypes
        genotypes = {}
        for reference, hits in hit_regions.items():

            # Get the length of the reference
            reference_len = len(sequence_database.get_refseq(reference))

            # Create the genotype object that will hold all hits for a
            # particular reference
            genotypes[reference] = Genotype(reference, reference_len, hits)

        return genotypes

    @property
    def predicted(self):
        return self._predicted

    @property
    def coverage(self):
        return self._coverage

    @property
    def identity(self):
        return self._identity

    @property
    def reference_id(self):
        return self._reference_id
