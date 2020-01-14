"""Linear algebra, including integral linear algebra."""

from fractions import Fraction
from utility import fracToInt, memorize

class RowSystem(object):
    """Manage a list of row vectors of integers. Find both integer and rational
    linear combinations of these vectors that sum to zero or another row
    vector.

    """
    def __init__(self, vecs):
        """vecs is a nonempty list of vectors. Each element of vecs is a vector
        of integers (rational inputs are possible, but they must have integer
        value, and will be converted to integers.

        """
        assert len(vecs) > 0
        self.ori_vecs = [list(vec) for vec in vecs]
        self.ori_vecs = [[fracToInt(n) for n in vec]
                         for vec in self.ori_vecs]
        self.num_row = len(self.ori_vecs)
        self.num_col = len(self.ori_vecs[0])
        self._rowReduce()
        self.num_reduced_row = len(self.reduced_vecs)
        assert self.num_reduced_row == len(self.pivot)

    def __str__(self):
        result = "Row system with rows:"
        for vec in self.ori_vecs:
            result += "\n"+"\t".join(["%d" % n for n in vec])
        return result

    def __repr__(self):
        return str(self)

    def _rowReduce(self):
        """Perform integeral row reduction. Produces three lists of vectors for
        this object:
        *. self.reduced_vecs: list of reduced vectors. These forms a basis for
        the subspace spanned by the original vectors, in the echelon form.
        *. self.reduced_comb: each row corresponds to a reduced vector, it
        expresses the reduced vector as a linear combination of the original
        vectors (note this is not necessarily unique).
        *. self.zero_comb: each row corresponds to a (integer) linear relation
        between the original vectors.
        Also produces self.pivot, which maps row number in reduced_vecs to the
        position of pivot in that row (first nonzero entry).

        """
        reduced_vecs = [list(vec) for vec in self.ori_vecs]
        # Expresses each row of reduced_vecs as a linear combination of rows
        # in ori_vecs. Initially this is just the identity.
        combs = [[0]*self.num_row for i in range(self.num_row)]
        for i in range(self.num_row):
            combs[i][i] = 1

        def swap_row(a, b):
            """Swaps row #a and #b of reduced_vecs."""
            for i in range(self.num_col):
                reduced_vecs[a][i], reduced_vecs[b][i] = \
                    reduced_vecs[b][i], reduced_vecs[a][i]
            for i in range(self.num_row):
                combs[a][i], combs[b][i] = combs[b][i], combs[a][i]

        def multiply_row(a, factor):
            """Multiply row #a by factor."""
            for i in range(self.num_col):
                reduced_vecs[a][i] *= factor
            for i in range(self.num_row):
                combs[a][i] *= factor

        def add_multiple(add_from, add_to, factor):
            """Add factor*(row #add_from) onto row #add_to."""
            for i in range(self.num_col):
                reduced_vecs[add_to][i] += (factor * reduced_vecs[add_from][i])
            for i in range(self.num_row):
                combs[add_to][i] += (factor * combs[add_from][i])

        cur_row, cur_col = 0, 0
        self.pivot = []
        while cur_row < self.num_row and cur_col < self.num_col:
            # Find an entry with minimum zero absolute value in this column, at
            # or below cur_row
            min_val, min_row = 0, -1
            for row in range(cur_row, self.num_row):
                cur_val = reduced_vecs[row][cur_col]
                if cur_val != 0 and (min_val == 0 or abs(cur_val) < min_val):
                    min_val, min_row = abs(cur_val), row
            # If all entries in this column at or below cur_row are zero, move
            # to next column
            if min_val == 0:
                cur_col += 1
                continue
            # Otherwise, swap minimum row onto cur_row, and correct for signs
            swap_row(cur_row, min_row)
            if reduced_vecs[cur_row][cur_col] < 0:
                multiply_row(cur_row, -1)
            assert min_val == reduced_vecs[cur_row][cur_col]
            # Reduce the remaining rows
            all_gcd = True
            for row in range(cur_row+1, self.num_row):
                factor = -reduced_vecs[row][cur_col] // min_val
                add_multiple(cur_row, row, factor)
                if reduced_vecs[row][cur_col] != 0:
                    all_gcd = False
            # If min_val is the gcd of remaining values in the column (so all
            # other values in the column are now zero), we are done. Otherwise,
            # we need to repeat this with a smaller min_val.
            if all_gcd:
                self.pivot.append(cur_col)
                cur_row += 1

        # At this point, all rows of self.reduced_vecs at or below cur_row
        # should be zero.
        self.reduced_vecs = reduced_vecs[0:cur_row]
        self.reduced_comb = combs[0:cur_row]
        self.zero_comb = combs[cur_row:]

    def getZeroComb(self):
        """Returns a list of linear relations between the original vectors."""
        return self.zero_comb

    def vecReduce(self, vec, use_rational = False):
        """Reduce the given vector to a standard form. If use_rational is set
        to True, then rational multiples of the original vectors are allowed.
        Otherwise only integer multiples are allowed to reduce the vector.
        If use_rational is set to false, only integer values are allowed
        in vec.

        Returns a tuple (comb, reduced_vec), where comb is a vector of length
        num_row, specifying the linear combination of original vectors that,
        when subtracted, will yield the standard form reduced_vec.

        """
        if not use_rational:
            vec = [fracToInt(n) for n in vec]
        comb = [0] * self.num_row
        cur_row = 0
        for col in range(self.num_col):
            if vec[col] == 0:
                continue
            while cur_row < self.num_reduced_row and self.pivot[cur_row] < col:
                cur_row += 1
            if cur_row >= self.num_reduced_row or self.pivot[cur_row] != col:
                continue
            pivot_val = self.reduced_vecs[cur_row][col]
            assert pivot_val > 0
            if use_rational:
                factor = Fraction(vec[col]) / pivot_val
            else:
                factor = vec[col] // pivot_val # integer division
            for i in range(self.num_col):
                vec[i] -= self.reduced_vecs[cur_row][i] * factor
            for i in range(self.num_row):
                comb[i] += self.reduced_comb[cur_row][i] * factor
            if use_rational:
                assert vec[col] == 0
            else:
                # Python definition of integer division means
                assert vec[col] >= 0 and vec[col] < pivot_val
        return (comb, vec)

    @memorize
    def reduceProfile(self, use_rational = False):
        """Returns a vector of length num_col, indicating at each position, the
        range of possible values returned by vec_reduce. A value of 0 means any
        integer is possible. Otherwise, a value of n > 0 indicates the range
        [0, n).

        """
        profile = [0] * self.num_col
        for row in range(self.num_reduced_row):
            col = self.pivot[row]
            profile[col] = self.reduced_vecs[row][col]
            assert profile[col] > 0
            if use_rational:
                profile[col] = 1
        return profile

    def shortForm(self, vec, use_rational = False):
        """Call vecReduce, but collect only those entries in the standard form
        that contain information (not always zero).

        """
        comb, reduced_vec = self.vecReduce(vec, use_rational)
        profile = self.reduceProfile(use_rational)
        result = [reduced_vec[i]
                  for i in range(self.num_col) if profile[i] != 1]
        return result

    def getComb(self, vec, use_rational = False):
        """Write the given vec as a linear combination of the original vectors.
        Note the answer is not necessarily uniquely specified. Return None if
        this is impossible.

        """
        comb, reduced_vec = self.vecReduce(vec, use_rational)
        if all([n == 0 for n in reduced_vec]):
            return comb
        else:
            return None

class F2RowSystem(object):
    """Linear algebra over field F2. Currently using just 0 and 1's (so the F2
    from utility.py is not involved.

    """
    def __init__(self, vecs):
        """Initialize with vecs, a matrix with 0/1 values."""
        assert len(vecs) > 0
        self.ori_vecs = [list(vec) for vec in vecs]
        self.num_row = len(self.ori_vecs)
        self.num_col = len(self.ori_vecs[0])
        self._rowReduce()
        self.num_reduced_row = len(self.reduced_vecs)
        assert self.num_reduced_row == len(self.pivot)

    def __str__(self):
        result = "Row system with rows:"
        for vec in self.ori_vecs:
            result += "\n"+"\t".join(["%d" % n for n in vec])
        return result

    def __repr__(self):
        return str(self)

    def _rowReduce(self):
        """Performs row reduction. Similar to the function of the same name in
        RowSystem.

        """
        reduced_vecs = [list(vec) for vec in self.ori_vecs]
        combs = [[0]*self.num_row for i in range(self.num_row)]
        for i in range(self.num_row):
            combs[i][i] = 1

        def swap_row(a, b):
            """Swaps row #a and #b of reduced_vecs."""
            for i in range(self.num_col):
                reduced_vecs[a][i], reduced_vecs[b][i] = \
                    reduced_vecs[b][i], reduced_vecs[a][i]
            for i in range(self.num_row):
                combs[a][i], combs[b][i] = combs[b][i], combs[a][i]

        def add_multiple(add_from, add_to):
            """Add row #add_from onto row #add_to."""
            for i in range(self.num_col):
                reduced_vecs[add_to][i] += reduced_vecs[add_from][i]
                reduced_vecs[add_to][i] %= 2
            for i in range(self.num_row):
                combs[add_to][i] += combs[add_from][i]
                combs[add_to][i] %= 2

        cur_row, cur_col = 0, 0
        self.pivot = []
        while cur_row < self.num_row and cur_col < self.num_col:
            # Find an entry with non-zero absolute value in this column, at or
            # below cur_row
            pivot_row = -1
            for row in range(cur_row, self.num_row):
                if reduced_vecs[row][cur_col] == 1:
                    pivot_row = row
                    break
            # If all entries in this column at or below cur_row are zero, move
            # to next column
            if pivot_row == -1:
                cur_col += 1
                continue
            # Otherwise, swap pivot row onto cur_row
            swap_row(cur_row, pivot_row)
            # Reduce the remaining rows
            for row in range(cur_row+1, self.num_row):
                if reduced_vecs[row][cur_col] == 1:
                    add_multiple(cur_row, row)
            self.pivot.append(cur_col)
            cur_row += 1
            cur_col += 1

        # At this point, all rows of self.reduced_vecs at or below cur_row
        # should be zero.
        self.reduced_vecs = reduced_vecs[0:cur_row]
        self.reduced_comb = combs[0:cur_row]
        self.zero_comb = combs[cur_row:]

    def getZeroComb(self):
        """Returns a list of linear relations between the original vectors."""
        return self.zero_comb

    def vecReduce(self, vec):
        """Reduce the given vector to a standard form.

        Returns a tuple (comb, reduced_vec), where comb is a vector of length
        num_row, specifying the linear combination of original vectors that,
        when subtracted, will yield the standard form reduced_vec.

        """
        comb = [0] * self.num_row
        cur_row = 0
        for col in range(self.num_col):
            if vec[col] == 0:
                continue
            while cur_row < self.num_reduced_row and self.pivot[cur_row] < col:
                cur_row += 1
            if cur_row >= self.num_reduced_row or self.pivot[cur_row] != col:
                continue
            pivot_val = self.reduced_vecs[cur_row][col]
            assert pivot_val == 1
            for i in range(self.num_col):
                vec[i] += self.reduced_vecs[cur_row][i]
                vec[i] %= 2
            for i in range(self.num_row):
                comb[i] += self.reduced_comb[cur_row][i]
                comb[i] %= 2
        return (comb, vec)

    def getComb(self, vec):
        """Write the given vec as a linear combination of the original vectors.
        Note the answer is not necessarily uniquely specified. Return None if
        this is impossible.

        """
        comb, reduced_vec = self.vecReduce(vec)
        if all([n == 0 for n in reduced_vec]):
            return comb
        else:
            return None
