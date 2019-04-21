###################################################################
#
# Fancy tools!
# 
# Author: Milan Patel
# Contact: mpatel5@cdc.gov
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

class DecisionTree(object):

    def __init__(self):
        self._attr = None
        self._value = None
        self._tree = {}

    def update(self, branch):
        # Attrs should come in transposed and sorted by column header
        # values of attrs can only be True or False

        if not self._attr:
            self._attr = branch[0][0]

        else:
            assert self._attr == branch[0][0]

        if isinstance(branch[0][1], basestring):
            self._value = branch[0][1]

        elif branch[0][1] in self._tree:
            self._tree[branch[0][1]].update(branch[1:])

        else:
            self._tree[branch[0][1]] = DecisionTree()
            self._tree[branch[0][1]].update(branch[1:])

    def recurse(self, branch):

        if self._value:
            return self._value

        elif not branch[0][0] == self._attr:
            return

        else:
            return self._tree[branch[0][1]].recurse(branch[1:])

def edit_distance(c, n, matrix=False):
    c_l = len(c)
    n_l = len(n)

    if not c_l:
        return n_l

    if not n_l:
        return c_l

    m = [[0]*(n_l+1) for _ in range(c_l+1)]

    for i in range(c_l+1):
        m[i][0] = i

    for j in range(n_l+1):
        m[0][j] = j

    for i in range(c_l):

        for j in range(n_l):

            if c[i] == n[j]:

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
        return m[c_l][n_l]

def get_alignment(seq1, seq2):
    if not seq1 or not seq2:
        return
    
    m = edit_distance(seq1, seq2, matrix=True)

    i = len(m) - 1
    j = len(m[0]) - 1

    seq1_aln = []
    seq2_aln = []

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
