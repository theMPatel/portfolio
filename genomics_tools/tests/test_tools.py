###################################################################
#
# Tests for the tools module
#
# Author: Milan Patel
# Contact: https://github.com/theMPatel
# Version 1.0
#
###################################################################

import io
import pytest

from genomics_tools.tools.tools import check_b64encoded
from genomics_tools.tools.tools import check_mismatches
from genomics_tools.tools.tools import codon_translation
from genomics_tools.tools.tools import fasta_iterator
from genomics_tools.tools.tools import get_all_file_exts
from genomics_tools.tools.tools import get_non_iupac
from genomics_tools.tools.tools import is_fasta
from genomics_tools.tools.tools import reverse_complement

class TestGenericTools:

    @pytest.mark.parametrize(
        "path, count", (
            ("/some/path/that/exists", 0),
            ("some.file.with.lots.o.exts", 5),
            ("/var/root/archive.tar.gz", 2),
            ("/etc/init.d/test.gz", 1)
        )
    )
    def test_get_all_extensions(self, path, count):
        root, exts = get_all_file_exts(path)
        assert len(exts) == count

    @pytest.mark.parametrize(
        "nucleotides, iupac_code",(
            ("A", 'A'),
            ("C", 'C'),
            ("G", 'G'),
            ("T", 'T'),
            ("AG", 'R'),
            ("CT", 'Y'),
            ("GC", 'S'),
            ("AT", 'W'),
            ("GT", 'K'),
            ("AC", 'M'),
            ("CA", 'M'),
            ("YN", None),
            ("N",  None),
            ("CGT", 'B'),
            ("AGT", 'D'),
            ("ACT", 'H'),
            ("ACG", 'V'),
            ("ACGT", 'N'),
        )
    )
    def test_get_correct_iupac(self, nucleotides, iupac_code):
        fset = frozenset(nucleotides)
        code = get_non_iupac(fset)
        assert code == iupac_code

    @pytest.mark.parametrize(
        "string, expected", (
            (b'YWxrc2RqZmxranNk\n', True),
            (b'CgoKCgoKXA==\n', True),
            (b'bHNsc2xzbHNsc2xzbHM=\n', True),
            (b'hahahaha', True),
            (b'obvislyhouldwork', True),
            (b'ancdd', False)
        )
    )
    def test_check_b64_encoded(self, string, expected):
        is_encoded = check_b64encoded(string)
        assert is_encoded == expected

    @pytest.mark.parametrize(
        "string, expected", (
            ("test.fna", True),
            ("test.fna.fasta", True),
            ("test.fsa", True),
            ("test.fsa.gz", True),
            ("test.tar.gz", False)
        )
    )
    def test_check_fast_filename(self, string, expected):
        assert is_fasta(string) == expected

    def test_fasta_iterator(self):
        lines = (
                ">SEQUENCE_1\n"
                "MTEITAAMVKELRESTGAGMMDCKNALSETNGDFDKAVQLLREKGLGKAAKKADRLAAEG\n"
                "LVSVKVSDDFTIAAMRPSYLSYEDLDMTFVENEYKALVAELEKENEERRRLKDPNKPEHK\n"
                "IPQFASRKQLSDAILKEAEEKIKEELKAQGKPEKIWDNIIPGKMNSFIADNSQLDSKLTL\n"
                "MGQFYVMDDKKTVEQVIAEKEKEFGGKIKIVEFICFEVGEGLEKKTEDFAAEVAAQL\n"
                ">SEQUENCE_2\n"
                "SATVSEINSETDFVAKNDQFIALTKDTTAHIQSNSLQSVEELHSSTINGVKFEEYLKSQI\n"
                "ATIGENLVVRRFATLKAGANGVVNGYIHTNGRVGVVIAAACDSAEVASKSRDLLRQICMH\n"
            )
        f = io.StringIO(lines)
        records = dict(fasta_iterator(f))
        assert len(records) == 2

    @pytest.mark.parametrize(
        "codon, expected", (
            ("GCT", "A"),
            ("TTT", "F"),
            ("ATG", "M"),
            ("TGG", "W"),
            ("TAT", "Y"),
            ("TAC", "Y"),
            ("GGT", "G"),
            ("CTA", "L"),
            ("AAC", "N")
        )
    )
    def test_codon_translate(self, codon, expected):
        translated = codon_translation(codon)
        assert translated == expected

    @pytest.mark.parametrize(
        "sequence, expected", (
            ("ACTGS", "SCAGT"),
            ("AAAAAAA", "TTTTTTT"),
            ("AVMKYAC", "GTRMKBT")
        )
    )
    def test_reverse_complement(self, sequence, expected):
        reversed_seq = reverse_complement(sequence)
        assert reversed_seq == expected

    @pytest.mark.parametrize(
        "seq1, seq2, expected", (
            ("ACTGS", "ACTGS", 0),
            ("AATAAG", "ATTAAG", 1),
            ("AAAAA", "TTTTT", 5)
        )
    )
    def test_sequence_mismatches(self, seq1, seq2, expected):
        mismatches = check_mismatches(seq1, seq2)
        assert mismatches == expected