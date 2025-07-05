"""Producing type DA structures for cobordisms, using local actions."""

from .arcslide import Arcslide
from .arcslideda import ArcslideDA
from .autocompleteda import autoCompleteDA
from .cobordism import Cobordism
from .cobordism import LEFT, RIGHT
from .dastructure import ComposedDAStructure, DAStructure, SimpleDAGenerator
from .extendbyid import ExtendedDAStructure, LocalDAStructure
from .hdiagram import getCobordismDiagramLeft, getSimpleCobordismDiagram
from .localpmc import LocalIdempotent, LocalStrandAlgebra, PMCSplitting
from .pmc import PMC
from .pmc import linearPMC
from .utility import subset
from .utility import F2
from .data import data_file_contents
import ast
import itertools

class CobordismDALeft(ExtendedDAStructure):
    """Responsible for producing a type DA structure for the cobordism with
    larger PMC on the left, using local actions.

    """
    def __init__(self, cob):
        """Specifies the cobordism to use. cob should be of type Cobordism,
        with cob.side == LEFT.

        """
        assert cob.side == LEFT

        self.cob = cob
        self.genus, self.c_pair, self.n = cob.genus, cob.c_pair, cob.n
        self.start_pmc, self.end_pmc = cob.start_pmc, cob.end_pmc
        self.large_pmc, self.small_pmc = cob.large_pmc, cob.small_pmc
        self.is_degenerate = cob.is_degenerate
        self.c1, self.c2 = cob.c1, cob.c2

        # Four possible local cases
        self.MIDDLE, self.NEXT_TOP, self.TOP, self.BOTTOM = 0, 1, 2, 3

        if not self.is_degenerate:
            self.d, self.u = self.cob.d, self.cob.u
            self.to_s, self.pair_to_s = self.cob.to_s, self.cob.pair_to_s
            self.d_pair, self.u_pair = self.cob.d_pair, self.cob.u_pair
            self.du_pair = self.cob.du_pair  # refers to the smaller PMC.
            self.sd, self.su = self.small_pmc.pairs[self.du_pair]

            u1, u2 = self.large_pmc.pairs[self.u_pair]
            assert u1 == self.u
            if self.c_pair < 2*self.genus-2:
                # First case: c-pair at least two pairs from top
                self.case = self.MIDDLE
                self.splitting_small = PMCSplitting(
                    self.small_pmc, [(self.su-1, self.su)])
                self.splitting_large = PMCSplitting(
                    self.large_pmc, [(self.c1, u2)])
            else:
                # Second case: c-pair is the next-to-topmost pair
                assert self.c_pair == 2*self.genus-2
                assert u2 == u1 + 2
                self.case = self.NEXT_TOP
                self.splitting_small = PMCSplitting(
                    self.small_pmc, [(self.su, self.su)])
                self.splitting_large = PMCSplitting(
                    self.large_pmc, [(self.c1, u2)])
        else:
            assert self.c2 == self.c1 + 2
            if self.c_pair == 2*self.genus-1:
                # Third case: degenerate at top
                self.case = self.TOP
                self.splitting_small = PMCSplitting(
                    self.small_pmc, [(self.n-5, self.n-5)])  # topmost point
                self.splitting_large = PMCSplitting(
                    self.large_pmc, [(self.c1-2, self.c2)])
            else:
                # Fourth case: degenerate at bottom
                assert self.c_pair == 0
                self.case = self.BOTTOM
                self.splitting_small = PMCSplitting(self.small_pmc, [(0, 0)])
                self.splitting_large = PMCSplitting(self.large_pmc, [(0, 4)])

        # Assumes the LEFT case
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

        # With generators set, add grading. Any generator can serve as base_gen
        for gen in self.generators:
            base_gen = gen
            break
        self.registerHDiagram(getCobordismDiagramLeft(self.cob), base_gen)

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

        arrow_patterns = dict()
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
            # 0 is the c-pair. u-d pairs are 1 and 2 at left and 1 at
            # right. 3 at left corresponds to 0 at right is free.
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
        for coeffs_a in list(arrow_patterns.keys()):
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

class SimpleCobordismDA(ExtendedDAStructure):
    """Cobordisms where a genus-1 split PMC is added somewhere in the starting
    PMC (on the left).

    """
    def __init__(self, start_pmc, insert_pos):
        """Specifies the starting PMC, and the point of inserting the genus-1
        split PMC.
        - insert_pos = 0: insert at bottom
        - insert_pos = start_pmc.n: insert at top

        """
        self.start_pmc = start_pmc
        self.insert_pos = insert_pos

        # Prepare end_pmc
        translate = dict()
        for p in range(self.start_pmc.n):
            if p < insert_pos:
                translate[p] = p
            else:  # p >= insert_pos
                translate[p] = p + 4
        self.end_pmc = PMC(
            [(translate[p], translate[q]) for p, q in self.start_pmc.pairs] +
            [(insert_pos, insert_pos+2), (insert_pos+1, insert_pos+3)])

        assert insert_pos >= 1 and insert_pos < start_pmc.n

        # Possible cases
        self.MIDDLE, self.NEXT_BOTTOM = 0, 1
        if insert_pos >= 2:
            self.case = self.MIDDLE
        else:
            self.case = self.NEXT_BOTTOM

        self.splitting1 = PMCSplitting(
            self.start_pmc, [(insert_pos-1, insert_pos-1)])
        self.splitting2 = PMCSplitting(
            self.end_pmc, [(insert_pos-1, insert_pos+3)])

        self.local_pmc1 = self.splitting1.local_pmc
        self.local_pmc2 = self.splitting2.local_pmc

        # Required so the left to right transition on the outside can proceed.
        assert self.splitting1.outer_pmc == self.splitting2.outer_pmc

        self.local_da = self.getLocalDAStructure()

        # Uncomment to use autoComplete. Note limit on len(coeffs_a) is 5.
        # if self.case == self.MIDDLE:
        #     autoCompleteDA(self.local_da, [0, 1])

        # Initiate the ExtendedDAStructure
        ExtendedDAStructure.__init__(
            self, self.local_da, self.splitting1, self.splitting2)

        # With generators set, add grading. Choose a generator of class zero as
        # base generator.
        for gen in self.generators:
            if gen.name[0] == "0":
                base_gen = gen
                break
        self.registerHDiagram(
            getSimpleCobordismDiagram(start_pmc, insert_pos), base_gen)

    def getLocalDAStructure(self):
        """Returns the local type DA structure associated to this simple
        cobordism.

        """
        # Construct arrow_patterns
        if self.case == self.MIDDLE:
            patterns_raw = self._patterns_middle()
        else:
            patterns_raw = self._patterns_next_bottom()

        arrow_patterns = dict()
        for pattern in patterns_raw:
            start_class, end_class = pattern[0], pattern[1]
            coeffs_a = []
            for i in range(2, len(pattern)-1):
                coeffs_a.append(self.local_pmc2.sd(pattern[i]))
            key = (start_class, end_class, tuple(coeffs_a))
            if key not in arrow_patterns:
                arrow_patterns[key] = []
            arrow_patterns[key].append(self.local_pmc1.sd(pattern[-1]))

        alg1 = LocalStrandAlgebra(F2, self.local_pmc1)
        alg2 = LocalStrandAlgebra(F2, self.local_pmc2)
        local_da = LocalDAStructure(F2, alg1, alg2)

        # The original part
        da_idems = [([], [2]), ([0], [0, 2])]
        for i in range(len(da_idems)):
            l_idem, r_idem = da_idems[i]
            local_da.addGenerator(SimpleDAGenerator(
                local_da, LocalIdempotent(self.local_pmc1, l_idem),
                LocalIdempotent(self.local_pmc2, r_idem), "0_%d" % i))

        # Part added due to finger-push
        da_idems_id = [([], [1]), ([0], [0, 1])]
        for i in range(len(da_idems_id)):
            l_idem, r_idem = da_idems_id[i]
            for gen_type in [1, 2]:
                local_da.addGenerator(SimpleDAGenerator(
                    local_da,
                    LocalIdempotent(self.local_pmc1, l_idem),
                    LocalIdempotent(self.local_pmc2, r_idem),
                    "%d_%d" % (gen_type, i)))

        mod_gens = local_da.getGenerators()

        # Manually take care of u_maps
        single_idems1 = local_da.single_idems1
        single_idems2 = local_da.single_idems2
        for i in range(len(single_idems1)):
            i1, i2 = single_idems1[i], single_idems2[i]
            for local_gen in mod_gens:
                idem1, idem2 = local_gen.idem1, local_gen.idem2
                if i1 in idem1 and i2 in idem2:
                    # local_gen is eligible for u_maps[i]
                    target_idem1 = idem1.removeSingleHor([i1])
                    target_idem2 = idem2.removeSingleHor([i2])
                    target_gen = [target for target in mod_gens
                                  if target.idem1 == target_idem1 and
                                  target.idem2 == target_idem2 and
                                  target.name[0] == local_gen.name[0]]
                    assert len(target_gen) == 1
                    local_da.add_u_map(i, local_gen, target_gen[0])

        # Check all u_map have been filled
        local_da.auto_u_map()

        # Add arrows according to arrow_pattern
        for key in list(arrow_patterns.keys()):
            start_class, end_class, coeffs_a = key
            if len(coeffs_a) == 1 and coeffs_a[0].isIdempotent():
                continue
            for coeff_d in arrow_patterns[key]:
                used = False
                for x, y in itertools.product(mod_gens, mod_gens):
                    if x.name[0] == "%d" % start_class and \
                       y.name[0] == "%d" % end_class and \
                       DAStructure.idemMatchDA(x, y, coeff_d, coeffs_a):
                        local_da.addDelta(x, y, coeff_d, coeffs_a, 1)
                        used = True
                if not used:
                    print("Warning: unused arrow: %s %s" % \
                        (coeffs_a, coeff_d))

        return local_da

    def _patterns_next_bottom(self):
        """Patterns in the middle case."""
        # - Local PMC at left (D-side) is 0-1*.
        # - Local PMC at right (A-side) is 0-1-2-3-4-5*, with pairs (1, 3)
        #   and (2, 4).

        # Simply take _pattern_middle and subtract one.
        patterns_raw_middle = self._patterns_middle()
        patterns_raw = []

        def translate(lst):
            """Similar to translate in arcslideda.py."""
            result = []
            for entry in lst:
                if isinstance(entry, int):
                    result.append(entry-1)
                    if result[-1] == -1:
                        return None
                else:  # entry must be a pair
                    result.append((entry[0]-1, entry[1]-1))
                    if result[-1][0] == -1 or result[-1][1] == -1:
                        return None
            return result

        for pattern in patterns_raw_middle:
            new_pattern = [translate(coeff) for coeff in pattern[2:]]
            if all([entry != None for entry in new_pattern]):
                patterns_raw.append([pattern[0], pattern[1]] + new_pattern)

        return patterns_raw

    def _patterns_middle(self):
        """Patterns in the middle case."""
        # - Local PMC at left (D-side) is 0*-1-2*.
        # - Local PMC at right (A-side) is 0*-1-2-3-4-5-6*, with pairs (2, 4)
        #   and (3, 5).
        input_patterns = data_file_contents("simplecob_arrows.data")
        patterns_raw = ast.literal_eval(input_patterns)
        return patterns_raw

class CobordismDARight(ComposedDAStructure):
    """Produces the type DA structure for a cobordism with larger PMC on the
    right, as the boxed tensor product of a SimpleCobordism with several
    arcslides.

    """
    def __init__(self, cob):
        """Specifies the cobordism to use. cob should be of type Cobordism,
        with cob.side == RIGHT.

        """
        assert cob.side == RIGHT

        self.n, self.genus, self.c_pair = cob.n, cob.genus, cob.c_pair
        start_pmc = linearPMC(self.genus-1)

        # Divide into four cases like in CobordismDALeft
        if self.c_pair == 0:
            # Bottom case:
            self.start_da = SimpleCobordismDA(start_pmc, 1)
            self.slides = [(0, 1)]
        elif self.c_pair == 2*self.genus-2:
            # Next to top case:
            self.start_da = SimpleCobordismDA(start_pmc, self.n-5)
            self.slides = [(self.n-1, self.n-2)]
        elif self.c_pair == 2*self.genus-1:
            # Top case:
            self.start_da = SimpleCobordismDA(start_pmc, self.n-5)
            self.slides = [(self.n-2, self.n-3), (self.n-1, self.n-2)]
        else:
            # Middle case
            self.c1 = cob.c1
            self.start_da = SimpleCobordismDA(start_pmc, self.c1)
            self.slides = [(self.c1+4, self.c1+3), (self.c1+1, self.c1),
                           (self.c1+2, self.c1+1), (self.c1+5, self.c1+4)]

        self.da_list = [self.start_da]
        for b1, c1 in self.slides:
            self.da_list.append(
                ArcslideDA(Arcslide(self.da_list[-1].pmc2, b1, c1)))

        ComposedDAStructure.__init__(self, self.da_list)
