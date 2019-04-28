###################################################################
#
# Fancy tools!
# 
# This tool is a sample and distillation of the real application
# hosted at: https://github.com/theMPatel/functional_genomics_tools
#
# Author: Milan Patel
# Contact: https://github.com/theMPatel
# Version 1.0
#
###################################################################

class Disjointset(object):
    # Disjoint sets object. 
    # Balances by:
    #   -> Rank heuristic
    #   -> Path compression

    def __init__(self, size):
        self._size = size
        self._rank = [0] * self._size
        self._parents = list(range(size))

    def get_parent(self, index):
        # This will flatten the disjoint set so that all
        # children point to the true parent
        if self._parents[index] != index:
            self._parents[index] = self.get_parent(self._parents[index])

        return self._parents[index]

    def merge(self, source, destination):
        # Merges based on rank heuristic
        real_source = self.get_parent(source)
        real_destination = self.get_parent(destination)

        if real_source == real_destination:
            return

        if self._rank[real_destination] > self._rank[real_source]:
            self._parents[real_source] = real_destination

        else:
            self._parents[real_destination] = real_source

            if self._rank[real_destination] == self._rank[real_source]:
                self._rank[real_source] += 1

    @property
    def size(self):
        return self._size

def binary_search(arr, value, lo, hi):
    """
    Use this function to do a binary search of your list.
    This works only on a list of non-duplicates, which is to
    say that if you have duplicates, you will get the index
    for only one of them. 

    Returns:
      - The index of the search value
      - The index of the largest element
            smaller than the requested search value
      - If nothing was found then -1 is returned
    
    :param arr: The array to search
    :param value: Desired value
    :param lo: Lower index bound
    :param hi: Higher index bound
    """

    if lo > hi:
        return hi

    mid = (hi+lo) // 2

    if value == arr[mid]:
        return mid

    elif value < arr[mid]:
        return binary_search(arr, value, lo, mid-1)

    else:
        return binary_search(arr, value, mid+1, hi)
        
def edit_distance(string_a, string_b, matrix=False):
    """
    Returns the minimum edit distance between two strings.
    :param string_a: The first string
    :param string_b: The second string
    :param matrix: Whether you would like the calculated matrix
        of distances.
    :returns: The minimum distance
    :returns: The calculated matrix
    :rtype: int | list
    """

    string_a_len = len(string_a)
    string_b_len = len(string_b)

    if not string_a_len:
        return string_b_len

    if not string_b_len:
        return string_a_len

    m = [[0]*(string_b_len+1) for _ in range(string_a_len+1)]

    for i in range(string_a_len+1):
        m[i][0] = i

    for j in range(string_b_len+1):
        m[0][j] = j

    for i in range(string_a_len):
        for j in range(string_b_len):

            if string_a[i] == string_b[j]:
                m[i+1][j+1] = m[i][j]

            else:

                m[i+1][j+1] = min(
                    m[i][j],
                    m[i+1][j],
                    m[i][j+1]
                ) + 1

    if matrix:
        return m
    else:
        return m[string_a_len][string_b_len]

def get_alignment(seq1, seq2):
    """
    Aligns two sequences together via edit distance. This is 
    primarily for *global* alignment, if you want local alignment
    you should use something different.

    :param seq1: The first sequence
    :param seq2: The second sequence
    """

    if not seq1 or not seq2:
        return
    
    m = edit_distance(seq1, seq2, matrix=True)

    i = len(m) - 1
    j = len(m[0]) - 1

    seq1_aln = []
    seq2_aln = []

    # Walk the returned matrix and build up the sequence.
    while i >= 1 and j >= 1:

        v = m[i][j]

        up_left = m[i - 1][j - 1]
        left = m[i][j - 1]
        up = m[i - 1][j]

        if all(up_left == x for x in [left, up]):
            i -= 1
            j -= 1
            seq1_aln.append(seq1[i])
            seq2_aln.append(seq2[j])
        
        else:
            m_val = min([up_left, left, up])

            if up_left == m_val:
                i -= 1
                j -= 1
                seq1_aln.append(seq1[i])
                seq2_aln.append(seq2[j])

            # Use dashes for insertions/deletions
            elif left == m_val:
                j -= 1
                seq2_aln.append(seq2[j])
                seq1_aln.append('-')
            else:
                i -= 1
                seq1_aln.append(seq1[i])
                seq2_aln.append('-')

    while i >= 1:
        i -= 1
        seq1_aln.append(seq1[i])
        seq2_aln.append('-')

    while j >= 1:
        j -= 1
        seq1_aln.append('-')
        seq2_aln.append(seq2[j])

    return list(reversed(seq1_aln)), list(reversed(seq2_aln)), m[-1][-1]

def pretty_aln(a, b, wrap=70):
    """
    Takes two strings that are presumed to be the same length
    and already aligned and creates a pretty output of the
    sequences that's easier for a human to read.
    Output would look like below:

    >>> pretty_aln(string_a, string_b)
    CCATGGTGACTCGGCGGTCTA
    |||||||||| ||||||| ||
    CCATGGTGACGCGGCGGTTTA
    
    :param a: String to put on top.
    :param b: String to put on bottom.
    """
    
    pipe_str = []
    for x, y in zip(a,b):
        
        if x==y:
            pipe_str.append('|')
        else:
            pipe_str.append(' ')
    
    x = [a[i:i+wrap] for i in range(0, len(a), wrap)]
    p = [pipe_str[i:i+wrap] for i in range(0, len(pipe_str), wrap)]
    y = [b[i:i+wrap] for i in range(0, len(b), wrap)]
    
    final = []
    
    for combo in zip(x, p, y):
        
        new_c = list(map(''.join, combo))
        new_c = '\n'.join(new_c)
        final.append(new_c)
    
    return '\n\n'.join(final)
