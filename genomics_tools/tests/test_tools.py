###################################################################
#
# Tests for the tools module
#
# Author: Milan Patel
# Contact: https://github.com/theMPatel
# Version 1.0
#
###################################################################

import pytest

from genomics_tools.tools.tools import get_all_file_exts
from genomics_tools.tools.tools import get_non_iupac
from genomics_tools.tools.tools import check_b64encoded
from genomics_tools.tools.tools import is_fasta

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