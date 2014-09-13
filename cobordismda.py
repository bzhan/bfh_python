"""Producing type DA structures for cobordisms, using local actions."""

from autocompleteda import autoCompleteDA
from cobordism import Cobordism
from cobordism import LEFT, RIGHT
from dastructure import DAStructure, SimpleDAGenerator
from extendbyid import ExtendedDAStructure, LocalDAStructure
from localpmc import LocalIdempotent, LocalStrandAlgebra, PMCSplitting
from utility import subset
from utility import F2
import itertools

class CobordismDA(ExtendedDAStructure):
    """Responsible for producing a type DA structure for the cobordism, using
    local actions.

    """
    def __init__(self, cob):
        """Specifies the cobordism to use. cobrodism should be of type
        Cobordism.

        """
        self.cob = cob
        self.genus, self.c_pair = self.cob.genus, self.cob.c_pair
        self.n, self.side = self.cob.n, self.cob.side
        self.start_pmc, self.end_pmc = self.cob.start_pmc, self.cob.end_pmc
        self.large_pmc, self.small_pmc = self.cob.large_pmc, self.cob.small_pmc
        self.c1, self.c2 = self.cob.c1, self.cob.c2

        # Four possible local cases
        self.MIDDLE, self.NEXT_TOP, self.TOP, self.BOTTOM = 0, 1, 2, 3

        if self.c2 == self.c1 + 3:
            self.d, self.u = self.cob.d, self.cob.u
            self.to_s, self.pair_to_s = self.cob.to_s, self.cob.pair_to_s
            self.d_pair, self.u_pair = self.cob.d_pair, self.cob.u_pair
            self.du_pair = self.cob.du_pair  # refers to the smaller PMC.
            self.sd, self.su = self.small_pmc.pairs[self.du_pair]

            u1, u2 = self.large_pmc.pairs[self.u_pair]
            assert u1 == self.u
            if u2 == u1 + 3:
                # First case: c-pair at least two pairs from top
                self.case = self.MIDDLE
                self.splitting_small = PMCSplitting(
                    self.small_pmc, [(self.su-1, self.su)])
                self.splitting_large = PMCSplitting(
                    self.large_pmc, [(self.c1, u2)])
            else:
                # Second case: c-pair is the next-to-topmost pair
                assert u2 == u1 + 2
                self.case = self.NEXT_TOP
                self.splitting_small = PMCSplitting(
                    self.small_pmc, [(self.su, self.su)])
                self.splitting_large = PMCSplitting(
                    self.large_pmc, [(self.c1, u2)])
        else:
            assert self.c2 == self.c1 + 2
            if self.c2 == self.n-1:
                # Third case: degenerate at top
                self.case = self.TOP
                self.splitting_small = PMCSplitting(
                    self.small_pmc, [(self.n-5, self.n-5)])  # topmost point
                self.splitting_large = PMCSplitting(
                    self.large_pmc, [(self.c1-2, self.c2)])
            else:
                # Fourth case: degenerate at bottom
                assert self.c1 == 0
                self.case = self.BOTTOM
                self.splitting_small = PMCSplitting(self.small_pmc, [(0, 0)])
                self.splitting_large = PMCSplitting(self.large_pmc, [(0, 4)])

        assert self.side == LEFT
        self.splitting1 = self.splitting_large
        self.splitting2 = self.splitting_small
        self.local_pmc1 = self.splitting1.local_pmc
        self.local_pmc2 = self.splitting2.local_pmc

        # Required so the left to right transition on the outside can proceed.
        assert self.splitting1.outer_pmc == self.splitting2.outer_pmc

        self.local_da = self.getLocalDAStructure()

        if self.case == self.MIDDLE:
            autoCompleteDA(self.local_da, [2, 1, 0, 5, 6])
        elif self.case == self.NEXT_TOP:
            autoCompleteDA(self.local_da, [2, 1, 0])
        elif self.case == self.TOP:
            autoCompleteDA(self.local_da, [3, 1])
        else:
            autoCompleteDA(self.local_da, [0, 3])

        # Initiate the ExtendedDAStructure
        ExtendedDAStructure.__init__(
            self, self.local_da, self.splitting1, self.splitting2)

    def getLocalDAStructure(self):
        """Returns the local type DA structure associated to this cobordism."""
        # Compute the set of arrow patterns
        if self.case == self.MIDDLE:
            patterns_raw = self._patterns_left_middle()
        elif self.case == self.NEXT_TOP:
            patterns_raw = self._patterns_left_next_top()
        elif self.case == self.TOP:
            patterns_raw = self._patterns_left_top()
        else:
            assert self.case == self.BOTTOM
            patterns_raw = self._patterns_left_bottom()

        arrow_patterns = {}
        for pattern in patterns_raw:
            coeffs_a = []
            for i in range(len(pattern)-1):
                coeffs_a.append(self.local_pmc2.sd(pattern[i]))
            coeffs_a = tuple(coeffs_a)
            if coeffs_a not in arrow_patterns:
                arrow_patterns[coeffs_a] = []
            arrow_patterns[coeffs_a].append(self.local_pmc1.sd(pattern[-1]))

        # Start construction of the local DA structure.
        alg1 = LocalStrandAlgebra(F2, self.local_pmc1)
        alg2 = LocalStrandAlgebra(F2, self.local_pmc2)
        if self.case == self.MIDDLE:
            local_dastr = LocalDAStructure(
                F2, alg1, alg2, single_idems1 = [1, 3], single_idems2 = [1, 0])
        else:
            local_dastr = LocalDAStructure(F2, alg1, alg2)

        # Compute the set of local generators.
        if self.case == self.MIDDLE:
            da_idems = [([0], []), ([0, 2], [1]), ([0, 1], [1]),
                        ([0, 3], [0]), ([0, 2, 3], [0, 1]), ([0, 1, 3], [0, 1])]
        elif self.case == self.NEXT_TOP:
            da_idems = [([0], []), ([0, 2], [0]), ([0, 1], [0])]
        elif self.case == self.TOP:
            da_idems = [([2], []), ([1, 2], [0])]
        else:
            da_idems = [([0], []), ([0, 2], [0])]

        for i in range(len(da_idems)):
            l_idem, r_idem = da_idems[i]
            local_dastr.addGenerator(SimpleDAGenerator(
                local_dastr, LocalIdempotent(self.local_pmc1, l_idem),
                LocalIdempotent(self.local_pmc2, r_idem), "%d" % i))
        mod_gens = local_dastr.getGenerators()

        # After having added all generators, create u_map:
        local_dastr.auto_u_map()

        # Add arrows according to arrow_pattern.
        for coeffs_a in arrow_patterns.keys():
            if len(coeffs_a) == 1 and coeffs_a[0].isIdempotent():
                continue
            for coeff_d in arrow_patterns[coeffs_a]:
                for x, y in itertools.product(mod_gens, mod_gens):
                    if DAStructure.idemMatchDA(x, y, coeff_d, coeffs_a):
                        local_dastr.addDelta(x, y, coeff_d, coeffs_a, 1)

        return local_dastr

    @staticmethod
    def _patterns_left_middle():
        """Patterns for self.side == LEFT and self.case == self.MIDDLE."""
        # - Local PMC at left (D-side) is 0*-1-2-3-4-5-6-7*, with pairs
        #   0:(1, 4) and 2:(3, 6), and single points 1:(2,) and 3:(5,).
        # - Local PMC at right (A-side) is 0*-1-2-3*, with two single points.
        # - Single point 1:(2,) at left corresponds to single point 1:(2,) at
        #   right, and single point 3:(5,) at left corresponds to single point
        #   0:(1,) at right.
        patterns_raw = [
            ### Seeds
            ([1, (2, 3)],),
            ([(1, 2), (3, 4)],),
            ([(0, 1)], [1, (0, 5)]),
            ([(1, 2)], [1, (5, 6)]),
            ([(2, 3)], [1, (6, 7)]),
        ]
        return patterns_raw

    @staticmethod
    def _patterns_left_next_top():
        """Patterns for self.side == LEFT and self.case == self.NEXT_TOP."""
        # - Local PMC at left (D-side) is 0*-1-2-3-4-5, with pairs 0:(1, 4) and
        #   2:(3, 5), and single point 1:(2,).
        # - Local PMC at right (A-side) is 0*-1, with single point 1.
        patterns_raw = [
            ### Seeds
            ([1, (2, 3)],),
            ([(1, 2), (3, 4)],),
            ([(0, 1)], [1, (0, 5)]),
        ]
        return patterns_raw

    @staticmethod
    def _patterns_left_top():
        """Patterns for self.side == LEFT and self.case == self.TOP."""
        # - Local PMC at left (D-side) is 0*-1-2-3-4-5, with pairs 0:(1, 4) and
        #   2:(3, 5), and single point 1:(2,).
        # - Local PMC at right (A-side) is 0*-1, with single point 1.
        patterns_raw = [
            ### Seeds
            ([(3, 5)],),
            ([1, (3, 5)],),
            ([(0, 1)], [3, (0, 2)]),
        ]
        return patterns_raw

    @staticmethod
    def _patterns_left_bottom():
        """Patterns for self.side == LEFT and self.case == self.BOTTOM."""
        # - Local PMC at left (D-side) is 0-1-2-3-4-5*, with pairs 0:(0, 2) and
        #   1:(1, 4), and single point 2:(3,).
        # - Local PMC at right (A-side) is 0-1*, with single point 0.
        patterns_raw = [
            ### Seeds
            ([(0, 2)],),
            ([1, (0, 2)],),
            ([(0, 1)], [0, (3, 5)]),
        ]
        return patterns_raw
