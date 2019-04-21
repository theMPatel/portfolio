###################################################################
#
# Tools for assembly based genotyping prediction
# 
# Author: Milan Patel
# Contact: mpatel5@cdc.gov
# Version 1.0
#
###################################################################

from tools.environment import (
    log_message,
    log_error,
    log_progress,
    valid_dir
)

from tools.align import (
    BLASTSettings,
    create_blastdb,
    align_blast,
    align_blast_nodb,
    GenotypeHit
)

from tools.tools import (
    reverse_complement,
    codon_translation,
    binary_search
)

import os
from itertools import combinations, product, izip
from tools.fancy_tools import Disjointset
from collections import defaultdict, namedtuple


GenotypeRegion = namedtuple('GenotypeRegion', ['coverage', 'identity', 'locations'])

def presence_detector(sequence_database, query_path, cached_query, percent_identity,
    min_relative_coverage, min_merge_overlap, search_fragments, env):
    
    # Let's export the references
    log_message('Exporting references...')

    # Make the path
    reference_dir = os.path.join(env.localdir, 'blastdb')
    
    # Check to make sure that its a real dir
    valid_dir(reference_dir)

    reference_path = os.path.join(reference_dir, 'references.fasta')

    # Export the reference sequences for blast database creation
    sequence_database.export_sequences(reference_path)

    log_message('Successfully exported reference database...')

    # Create the path to the blast database
    blast_db_path = os.path.join(env.localdir, 'blastdb', 'db.fasta')

    log_message('Creating blast database...')

    # Create the blast database
    create_blastdb(reference_path, blast_db_path, env)

    # Log that we were successful
    log_message('Successfully created blast database!')

    # Create the blast settings so that we can run the thing!
    blast_settings = BLASTSettings(
        task = 'dc-megablast',
        identity = percent_identity,
        relative_minlen = 0,
        absolute_minlen = 0,
        include_sequences = False
        )

    log_message('BLASTing query genome against reference database')
    
    # Run the alignment
    results = align_blast(
        query_path,
        blast_db_path,
        blast_settings,
        env
    )

    # Determine the size of the contigs that we are working with
    contig_sizes = {contig:len(sequence) for \
        contig, sequence in cached_query.iteritems()}

    log_message('Determining genotype coverages...')
    
    # Create the hit objects
    regions = Genotype.find_regions(results, contig_sizes, sequence_database)

    log_message('Found {} potential genotypes!'.format(len(regions)))

    # Figure out which ones are actually worth keeping
    log_message('Determining accetable genotypes...')

    # Remove those that do not have a predicted genotype
    to_remove = set()
    for region in regions.itervalues():
        if not region.validate(
            percent_identity,
            min_relative_coverage,
            contig_sizes,
            search_fragments
        ):
            to_remove.add(region.reference_id)

    # Remove the genotypes that don't have anything
    for rm in to_remove:
        del regions[rm]

    log_message('After filtering, {} genotypes were retained'.format(
        len(regions)))

    # Eliminate any overlap between genotypes of different references
    # and return the accepted genotypes
    log_message('Determining optimal genotype(s)', 2)
    accepted = eliminate_overlap(regions, min_merge_overlap)

    log_message('After overlap analysis, {} genotypes were retained!'.format(
        len(accepted)))

    return accepted

def mutation_detector(sequence_database, query_path, percent_identity,
    min_relative_coverage, env):
    
    # Let's export the references
    log_message('Exporting references...')

    # Make the path
    reference_dir = os.path.join(env.localdir, 'blastdb')
    
    # Check to make sure that its a real dir
    valid_dir(reference_dir)

    reference_path = os.path.join(reference_dir, 'references.fasta')

    # Export the reference sequences for blast database creation
    sequence_database.export_sequences(reference_path)

    log_message('Successfully exported reference database...')

    # Create the path to the blast database
    blast_db_path = os.path.join(env.localdir, 'blastdb', 'db.fasta')

    log_message('Creating blast database...')

    # Create the blast database
    #create_blastdb(reference_path, blast_db_path, env)

    # Log that we were successful
    log_message('Successfully created blast database!')

    # Create the blast settings so that we can run the thing!
    # Return the sequences so that we can see if we've found
    # point mutations
    blast_settings = BLASTSettings(
        task = 'blastn',
        identity = percent_identity,
        relative_minlen = 0,
        absolute_minlen = 0,
        include_sequences = True
        )

    log_message('BLASTing query genome against reference database')
    
    # Run the alignment
    results = align_blast_nodb(
        query_path,
        reference_path,
        blast_settings,
        env
    )

    log_message('Successfully BLASTed query genome against reference database')

    log_message('Searching for mutations...')

    # Figure out the mutations that exist in the sequences
    interpretations = find_mutations(
        sequence_database,
        results,
        min_relative_coverage)

    log_message('Retained {} gene regions after gene analysis'.format(
        len(interpretations)))

    return interpretations

def find_mutations(sequence_database, results, min_relative_coverage):

    regions = defaultdict(list)

    for hit in results.hits:
        regions[hit.reference_id].append(hit)

    log_message('Found {} potential regions of interest'.format(
        str(len(regions))))

    # Store the found resistance:
    mutation_results = defaultdict(list)

    for reference, hits in regions.iteritems():

        targets = sequence_database.targets[reference]
        
        for hit in hits:

            if hit.relative_len < min_relative_coverage:
                continue

            ref_start = hit.reference_start
            ref_stop = hit.reference_stop
            query_start = hit.query_start
            query_stop = hit.query_stop

            hit_ref_seq = hit.reference_seq
            hit_query_seq = hit.query_seq

            # The sequences come back from blast as the reverse
            # complement of the reference, thus we need to reverse
            # them back to get a sequence in the same orientation
            # as the reference sequence ***IF*** the alignment is a
            # reverse alignment
            if not hit.forward:
                hit_ref_seq = reverse_complement(hit_ref_seq)
                hit_query_seq = reverse_complement(hit_query_seq)


            # Store the gap positions for each of the sequences.
            # They will naturally be sorted, we can do binary search
            # for figuring out if we need to offset the string indices
            gaps = hit.num_gap_opens
            ref_gaps = []
            query_gaps = []

            if gaps:
                # Deletions
                query_gaps = [i for i, s in enumerate(hit_query_seq) if s == '-']
                # Insertions
                ref_gaps = [i for i, s in enumerate(hit_ref_seq) if s == '-']


            for target in targets:

                codon_position = target.codon_position
                num_mutations = target.num_mutations_needed
                reference_codons = target.reference_codon
                reference_aas = target.reference_aa
                resistance_aas = target.resistance_aa
                resistances = target.resistance

                # The positional information for the coding genes
                # is stored as the position of the codon
                # i.e. codon 415 is nucleotide position 415*3
                if target.coding_gene:
                    # Get the postional information for this target
                    # These are the (1-indexed) values which we
                    # transform into the (0-indexed) position.
                    # You could leave the end postion without
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
                        continue

                    # Create the indices for the returned
                    # query string
                    string_start = start - ref_start
                    string_end = end - ref_start

                    # Determine how many gaps there are in the
                    # returned reference sequence (insertions)
                    max_ins = binary_search(
                        ref_gaps, string_start, 0, len(ref_gaps)-1)

                    ins_offset = max_ins + 1

                    if '-' not in reference_codons:
                        while hit_ref_seq[string_start+ins_offset] == '-' and \
                            string_end+ins_offset+1 < len(hit_ref_seq):
                            
                            ins_offset += 1

                    if len(hit_ref_seq) - 1 < string_end + ins_offset + 1:
                        continue

                    # Check to see that the reference sequence
                    # that is returned from blast
                    # does not have any gaps that would
                    # cause a frameshift
                    hit_ref_codon = hit_ref_seq[
                        string_start + ins_offset : string_end + ins_offset + 1]

                    # just in case
                    assert len(hit_ref_codon) == 3

                    if hit_ref_codon not in reference_codons:
                        # If this happens, it should be because
                        # there are insertions into your query
                        # genome that caused there to be gaps
                        # in your reference sequence within the
                        # codon itself
                        continue

                    # Determine how many gaps there are in the returned
                    # query (deletions)
                    max_dels = binary_search(
                        query_gaps, string_start, 0, len(query_gaps)-1)

                    del_offset = max_dels + 1

                    if '-' not in resistance_aas:
                        while hit_query_seq[string_start + del_offset] == '-' and \
                            string_end + del_offset + 1 < len(hit_query_seq):

                            del_offset += 1

                    if len(hit_query_seq) - 1 < string_end + del_offset + 1:
                        continue

                    hit_query_codon = hit_query_seq[
                        string_start + del_offset : string_end + del_offset + 1]

                    assert len(hit_query_codon) == 3

                    # Get the translation
                    query_translation = codon_translation(hit_query_codon)

                    if query_translation in resistance_aas:

                        results = {
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
                            'iscoding' : target.coding_gene
                        }


                        mutation_results[hit.reference_id].append(results)


                # The RNA genes and ampC promoter mutation positions are stored
                # as actual positions of the point mutations
                # since these do not actually code for genes
                else:

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
                        continue

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
                        continue

                    hit_ref_nucleotide = hit_ref_seq[string_start+ins_offset]

                    if hit_ref_nucleotide not in reference_codons:
                        # This should only occur in the case that
                        # your query sequence has insertions which
                        # caused the existence of gaps in your 
                        # reference sequence
                        continue

                    max_dels = binary_search(
                        query_gaps, string_start, 0, len(query_gaps)-1)

                    del_offset = max_dels + 1

                    if '-' not in resistance_aas:
                        while hit_query_seq[string_start + del_offset] == '-' and \
                            string_start + del_offset + 1 < len(hit_query_seq):

                            del_offset += 1

                    if len(hit_query_seq) - 1 < string_start + del_offset:
                        continue

                    # Get the query string nucleotide
                    hit_query_nucleotide = hit_query_seq[string_start+del_offset]

                    # Figure out if the hit is the same as a resistance
                    # nuc that is known
                    if hit_query_nucleotide in resistance_aas:

                        results = {
                            'locus' : target.gene_id,
                            'identity' : hit.identity,
                            'contig_id' : hit.query_id,
                            'query_codon': hit_query_nucleotide,
                            'position': codon_position,
                            'reference_nuc' : hit_ref_nucleotide,
                            'resistance' : resistances,
                            'aa_mutation' : '',
                            'hit': hit,
                            'iscoding' : target.coding_gene
                        }

                        mutation_results[hit.reference_id].append(results)

    return mutation_results

def encompassed(hit1, hit2, min_coverage):
    # Calculate whether one hit captures another:
    # 
    #         max(start)     min(end)
    # hit    .___________________.
    # hit  ._________________________.
    # 
    # hit1 and hit2 come in as lists
    # in case one is a fragment

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

    # Create a flat list of all of the hits
    regions = [hit for region in regions.itervalues() for \
        hit in region.predicted]

    # Create the disjoint sets object
    dset = Disjointset(len(regions))

    # Create the indices of the hits list
    for i, j in combinations(xrange(len(regions)), 2):

        hit1 = regions[i].locations
        hit2 = regions[j].locations

        if encompassed(hit1, hit2, min_merge_overlap):
            # If the two hits are overlapping, link them
            dset.merge(i,j)

    # Create a dictionary of the parents and the best hits:
    best_hits = defaultdict(set)

    for i in xrange(len(regions)):

        # Get the parent of hit at this index
        parent = dset.get_parent(i)

        # Add it to the set of children
        #   NOTE:
        # A parent is a child of itself,
        # so the childrens list will include all
        # of the parent's children plus itself
        best_hits[parent].add(i)

    # For each of the best parents, see if the children are better
    to_remove = []
    for parent, children in best_hits.iteritems():

        # Parent and children will be part of
        # childrens list
        hits_here = list(children)

        # Sort the regions based on identity:
        hits_here.sort(key=lambda x: -regions[x].identity)

        # Get rid of all the worst ones:
        to_remove.extend(hits_here[1:])

    # Sort the indices
    to_remove.sort()

    # Remove the genotypes
    for index in reversed(to_remove):
        del regions[index]

    # Return a dictionary of references and it's regions
    accepted = defaultdict(list)

    for region in regions:

        references_here = set()

        for hit in region.locations:
            references_here.add(hit.reference_id)

        if len(references_here) > 1:
            raise RuntimeError('This should never happen!')

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
            # The reason for calculating whether or not a hit
            # is at the edge is in case you have references
            # that are split over two contigs:
            # 
            #       contig1             contig2
            #                <25bp> <25bp>
            # --------------(------|------)-----------------
            #                   ^      ^
            #                   |      |
            #    ends within here      starts within here
            # 
            # 

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

                # Get the real coverage, incase calculated
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

        for reference, hits in hit_regions.iteritems():

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
