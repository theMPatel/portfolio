
def validate_coding_gene(hit, target, ref_gaps, query_gaps, **kwargs):

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
        return default_return

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
        return default_return

    hit_query_codon = hit_query_seq[
        string_start + del_offset : string_end + del_offset + 1]

    assert len(hit_query_codon) == 3

    ref_start_index_log = max(0, string_start + ins_offset - 10)
    ref_end_index_log = min(len(hit_ref_seq),
                            ref_start_index_log + 21)
    ref_slice = slice(ref_start_index_log, ref_end_index_log, 1)

    query_start_index_log = max(0, string_start + del_offset -10)
    query_end_index_log = min(len(hit_query_seq),
                            query_start_index_log + 21)
    query_slice = slice(query_start_index_log, query_end_index_log, 1)

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
    
    ref_start_index_log = max(0, string_start + ins_offset - 10)
    ref_end_index_log = min(len(hit_ref_seq),
                            ref_start_index_log + 21)
    ref_slice = slice(ref_start_index_log, ref_end_index_log, 1)

    query_start_index_log = max(0, string_start + del_offset -10)
    query_end_index_log = min(len(hit_query_seq),
                            query_start_index_log + 21)
    query_slice = slice(query_start_index_log, query_end_index_log, 1)

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

def find_mutations(sequence_database, results, min_relative_coverage):

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
                query_gaps = [i for i, s in enumerate(hit_query_seq) if s == '-']
                # Insertions
                ref_gaps = [i for i, s in enumerate(hit_ref_seq) if s == '-']

            for target in targets:

                if target.coding_gene:

                    results = validate_coding_gene(hit, target, ref_gaps,
                                                    query_gaps)
                    if results:
                        mutation_results[hit.reference_id].append(results)

                # The RNA genes and ampC promoter mutation positions are stored
                # as actual positions of the point mutations
                # since these do not actually code for genes
                else:

                    results = validaet_noncoding_gene(hit, target, ref_gaps,
                                                        query_gaps)

                    if results:
                        mutation_results[hit.reference_id].append(results)

    return mutation_results