###################################################################
#
# Tests for the fancy tools.
#
# Author: Milan Patel
# Contact: https://github.com/theMPatel
# Version 1.0
#
###################################################################

import pytest
from genomics_tools.tools.fancy_tools import edit_distance
from genomics_tools.tools.fancy_tools import Disjointset
from genomics_tools.tools.fancy_tools import binary_search

class TestFancyTools:

    @pytest.mark.parametrize(
        "stringA, stringB, expected", (
            ("ABC", "DEF", 3),
            ("ABC", "DBC", 1),
            ("IelsElsf", "", 8),
            ("ACTGGATTAC", "ACGGATTAC", 1)
        )
    )
    def test_valid_edit_distance(self, stringA, stringB, expected):
        distance = edit_distance(stringA, stringB)
        assert(distance == expected)

    def test_disjoint_set_api(self):
        Disjointset.get_parent
        Disjointset.merge
        Disjointset.size

    def test_disjoint_set(self):
        size = 50

        test_data = list(range(size))
        dset = Disjointset(len(test_data))

        for left in range(dset.size//2):
            right = -(left+1)+50

            dset.merge(left, right)

        for left in range(dset.size//2):
            right = -(left+1)+50

            assert dset.get_parent(left) == dset.get_parent(right)

        assert dset.get_parent(size//4) != dset.get_parent((size//4)+1)

    @pytest.mark.parametrize(
        "arr, target, expected", (
            ([8,10,18,21,22], 19, 2),
            ([0,1,2,3,4,5], 10, 5),
            ([87,99,129,347,900], 1, -1)
        )
    )
    def test_binary_search(self, arr, target, expected):
        index = binary_search(arr, target, 0, len(arr)-1)
        assert index == expected