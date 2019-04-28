###################################################################
#
# Tests for the align module
#
# Author: Milan Patel
# Contact: https://github.com/theMPatel
# Version 1.0
#
###################################################################

import io
import pytest

from genomics_tools.tools.align import GenotypeHit
from genomics_tools.tools.align import GenotypeResults

class TestAlign:

    def test_genotype_api(self):
        GenotypeHit.from_blast
        
        genotype_hit = GenotypeHit()
        genotype_hit.reference_id
        genotype_hit.query_id
        genotype_hit.query_start
        genotype_hit.query_stop
        genotype_hit.reference_start
        genotype_hit.reference_stop
        genotype_hit.forward
        genotype_hit.identity
        genotype_hit.absolute_len
        genotype_hit.relative_len
        genotype_hit.reference_len
        genotype_hit.num_mismatches
        genotype_hit.num_gap_opens
        genotype_hit.query_seq
        genotype_hit.reference_seq
        genotype_hit.bitscore
        genotype_hit.evalue
        genotype_hit.full_match

    def test_blast_hit_parsing(self):
        line = "CU928145.2\trpoB|4029\t92.000\t25\t2\t0\t4808419\t4808443\t735\t711\t0.68\t37.4\tGCCTTCCAGCACCAGTTCCATCTGC\tGCGTTCCGGCACCAGTTCCATCTGC"

        hit_obj = GenotypeHit.from_blast(line)

        assert hit_obj.query_id == "CU928145.2"
        assert hit_obj.identity == 0.92
        assert hit_obj.absolute_len == 25
        assert hit_obj.num_mismatches == 2
        assert hit_obj.query_seq == "GCCTTCCAGCACCAGTTCCATCTGC"
        assert hit_obj.reference_seq == "GCGTTCCGGCACCAGTTCCATCTGC"

    def test_genotype_results(self):
        line = "#BLASTN 2.9.0+\n# Query: CU928145.2 Escherichia coli 55989 chromosome, complete genome\nCU928145.2\trpoB|4029\t92.000\t25\t2\t0\t4808419\t4808443\t735\t711\t0.68\t37.4\tGCCTTCCAGCACCAGTTCCATCTGC\tGCGTTCCGGCACCAGTTCCATCTGC"
        f = io.StringIO(line)

        genotype_object = GenotypeResults().load_hits(f, 'blast')
        assert len(genotype_object.hits) == 1