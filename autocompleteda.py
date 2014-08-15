"""Auto-completion of arrows in a type DA structure, by solving certain
equations in linear algebra.

This module is used to produce arrows in the local type DA structure for
arcslides, in arcslidedatest.py.

"""

from algebra import E0
from dastructure import DAStructure
from extendbyid import LocalDAStructure
from linalg import F2RowSystem
from localpmc import LocalStrandDiagram
from utility import memorizeHash
import itertools
from Queue import Queue

class _DAArrow(tuple):
    """Structure representing a type DA arrow."""
    def __new__(cls, coeff_d, coeffs_a, source, target):
        return tuple.__new__(cls, (coeff_d, coeffs_a, source, target))

class AutoCompleteDAStructure:
    """A routine to complete arrows in type DA structures by solving certain
    equations in linear algebra.

    """
    def __init__(self, raw_da, d_side_order):
        """Initialize the needed values."""
        assert isinstance(raw_da, LocalDAStructure)
        self.raw_da = raw_da
        self.d_side_order = d_side_order
        self.mod_gens = self.raw_da.getGenerators()

        # Prepare single idems in raw_da as a list of pairs (idem_d, idem_a).
        self.single_idems = []
        for i in range(self.raw_da.num_single_idems):
            self.single_idems.append(
                (self.raw_da.single_idems1[i], self.raw_da.single_idems2[i]))

    def _getDerivedTwoStepArrows(self, arrows_base, arrow_new):
        """Find all ways of deriving two-step arrows from arrow_new, including
        by combining it with elements of arrows_base.

        """
        result = []
        coeff_d, coeffs_a, x, y = arrow_new
        if len(coeffs_a) > 0 and coeffs_a[0].isIdempotent():
            return result
        # Take anti-differential of one of coeffs_a
        for i in range(len(coeffs_a)):
            for anti_da, coeff in coeffs_a[i].antiDiff().items():
                result.append(_DAArrow(
                    coeff_d, coeffs_a[:i] + (anti_da,) + coeffs_a[i+1:], x, y))
            for (a, b), coeff in coeffs_a[i].factor().items():
                if a.isIdempotent() or b.isIdempotent():
                    continue
                result.append(_DAArrow(
                    coeff_d, coeffs_a[:i] + (a, b) + coeffs_a[i+1:], x, y))
        # Take differential of coeff_d
        for dd, coeff in coeff_d.diff().items():
            result.append(_DAArrow(dd, coeffs_a, x, y))
        # Multiply two together
        for coeff_d2, coeffs_a2, x2, y2 in arrows_base:
            if len(coeffs_a2) > 0 and coeffs_a2[0].isIdempotent():
                continue
            # One direction
            if y == x2 and coeff_d * coeff_d2 != E0:
                result.append(_DAArrow(
                    (coeff_d * coeff_d2).getElt(), coeffs_a + coeffs_a2, x, y2))
            # Other direction
            if y2 == x and coeff_d2 * coeff_d != E0:
                result.append(_DAArrow(
                    (coeff_d2 * coeff_d).getElt(), coeffs_a2 + coeffs_a, x2, y))
        return result

    def _getAltFactorizations(self, arrows_base, arrow_new):
        result = []
        coeff_d, coeffs_a, x, y = arrow_new
        # Take differential of one of coeffs_a
        for i in range(len(coeffs_a)):
            for da, coeff in coeffs_a[i].diff().items():
                result.append(_DAArrow(
                    coeff_d, coeffs_a[:i] + (da,) + coeffs_a[i+1:], x, y))
            if i > 0 and coeffs_a[i-1] * coeffs_a[i] != E0:
                result.append(_DAArrow(
                    coeff_d, coeffs_a[:i-1] + \
                    ((coeffs_a[i-1] * coeffs_a[i]).getElt(),) +
                    coeffs_a[i+1:], x, y))
        # Take anti-differential of coeff_d
        for anti_dd, coeff in coeff_d.antiDiff().items():
            result.append(_DAArrow(anti_dd, coeffs_a, x, y))
        # Split into two sequences
        for (a, b), coeff in coeff_d.factor().items():
            # Number of A-inputs in the first sequence
            for i in range(len(coeffs_a)):
                for z in self.mod_gens:  # intermediate generator
                    if (a, coeffs_a[:i], x, z) in arrows_base:
                        result.append(_DAArrow(b, coeffs_a[i:], z, y))
                    if (b, coeffs_a[i:], z, y) in arrows_base:
                        result.append(_DAArrow(a, coeffs_a[:i], x, z))
        return result

    def _getAltIdempotents(self, arrows):
        """Arrows is a list of tuples (coeff_d, coeffs_a).
        single_idems is a list of tuples (idem_d, idem_a), specifying the ID of
        single idempotents on the D-side and A-side.
        Returns a list of arrows that are different from those in the input by
        only the single idempotents.

        """
        def uses_idempotent(alg_gen, idem):
            return idem in alg_gen.left_idem or idem in alg_gen.right_idem

        def has_singlehor(alg_gen, idem):
            return idem in alg_gen.single_hor

        def add_singlehor(alg_gen, idem):
            return LocalStrandDiagram(
                alg_gen.parent, [idem] + list(alg_gen.left_idem),
                alg_gen.strands)

        def remove_singlehor(alg_gen, idem):
            return LocalStrandDiagram(
                alg_gen.parent, [i for i in alg_gen.left_idem if i != idem],
                alg_gen.strands)

        results = []
        for i in range(len(self.single_idems)):
            idem_d, idem_a = self.single_idems[i]
            for coeff_d, coeffs_a, x, y in results + arrows:
                if has_singlehor(coeff_d, idem_d) and \
                   all([has_singlehor(coeff, idem_a) for coeff in coeffs_a]):
                    new_arrow = _DAArrow(
                        remove_singlehor(coeff_d, idem_d),
                        tuple([remove_singlehor(coeff, idem_a)
                               for coeff in coeffs_a]),
                        self.raw_da.u_maps[i][x], self.raw_da.u_maps[i][y])
                    if not new_arrow in results + arrows:
                        results.append(new_arrow)
                if (not uses_idempotent(coeff_d, idem_d)) and \
                   all([not uses_idempotent(coeff, idem_a)
                        for coeff in coeffs_a]):
                    new_arrow = _DAArrow(
                        add_singlehor(coeff_d, idem_d),
                        tuple([add_singlehor(coeff, idem_a)
                               for coeff in coeffs_a]),
                        self.raw_da.uinv_maps[i][x],
                        self.raw_da.uinv_maps[i][y])
                    if not new_arrow in results + arrows:
                        results.append(new_arrow)
        return results

    def _autoCompleteByLinAlg(self, arrows_base, arrows_new):
        """Auto-complete arrows by solving a system of linear equations mod 2.
        Use only when it is clear that no two added arrows can be composed
        together.

        Returns the list of suggested arrows.

        """
        # First step: find all possible arrows, and all possible two-step
        # arrows.
        one_step_arrows = set()
        two_step_arrows = set()
        one_step_arrows_queue = Queue()
        two_step_arrows_queue = Queue()
        for arrow in arrows_new:
            one_step_arrows_queue.put(arrow)
            one_step_arrows.add(arrow)

        while not one_step_arrows_queue.empty() or \
              not two_step_arrows_queue.empty():
            if not one_step_arrows_queue.empty():
                cur_arrow = one_step_arrows_queue.get()
                for arrow in self._getDerivedTwoStepArrows(
                        arrows_base, cur_arrow):
                    if arrow not in two_step_arrows:
                        two_step_arrows.add(arrow)
                        two_step_arrows_queue.put(arrow)
                for arrow in self._getAltIdempotents([cur_arrow]):
                    if arrow not in one_step_arrows:
                        one_step_arrows.add(arrow)
                        one_step_arrows_queue.put(arrow)
            else:
                cur_arrow = two_step_arrows_queue.get()
                for arrow in self._getAltFactorizations(arrows_base, cur_arrow):
                    if arrow not in one_step_arrows:
                        one_step_arrows.add(arrow)
                        one_step_arrows_queue.put(arrow)

        # Combine some of the one-step arrows
        combined_one_step_arrows = []  # list of lists of arrows
        while len(one_step_arrows) != 0:
            arrow = one_step_arrows.pop()
            coeff_d, coeffs_a, gen_from, gen_to = arrow
            valid = True
            for idem_d, idem_a in self.single_idems:
                if idem_a in coeff_d.single_hor and \
                   not all(
                       idem_d in coeff_a.single_hor for coeff_a in coeffs_a):
                    valid = False
                    break
            if not valid:
                # These arrows should not be considered.
                continue

            alt_idems = self._getAltIdempotents([arrow])
            for alt_idem in alt_idems:
                one_step_arrows.remove(alt_idem)
            combined_one_step_arrows.append([arrow] + alt_idems)

        # Generate the matrix mapping from one-step arrows to two-step arrows
        two_step_arrows = list(two_step_arrows)
        matrix_map = [[0] * len(two_step_arrows)
                      for i in range(len(combined_one_step_arrows))]
        target_vec = [0] * len(two_step_arrows)
        for i in range(len(combined_one_step_arrows)):
            for one_step_arrow in combined_one_step_arrows[i]:
                derived_two_steps = self._getDerivedTwoStepArrows(
                    arrows_base, one_step_arrow)
                for j in range(len(two_step_arrows)):
                    if two_step_arrows[j] in derived_two_steps:
                        if one_step_arrow in arrows_new:
                            target_vec[j] += 1
                            target_vec[j] %= 2
                        else:
                            matrix_map[i][j] += 1
                            matrix_map[i][j] %= 2

        lin_sys = F2RowSystem(matrix_map)
        comb = lin_sys.getComb(target_vec)
        assert comb is not None

        result = []
        for i in range(len(combined_one_step_arrows)):
            if comb[i] != 0 and \
               combined_one_step_arrows[i][0] not in arrows_new:
                result.extend(combined_one_step_arrows[i])
        return result

    def complete(self):
        """Input raw_da is a local type DA structure with all generators and
        some arrows filled in. This function adds arrows to raw_da so that it
        satisfies the type DA structure equation, as well as the smearability
        condition.

        """
        # Complete in stages, each adding one interval on the D side, according
        # to d_side_order.
        for i in range(len(self.d_side_order)):
            # Lists of _DAArrow objects.
            arrows_base, arrows_seed = [], []
            # Figure out the base arrows and new arrows in this stage.
            for (gen_from, coeffs_a), target in self.raw_da.da_action.items():
                for (coeff_d, gen_to), ring_coeff in target.items():
                    arrow = _DAArrow(coeff_d, coeffs_a, gen_from, gen_to)
                    if all([coeff_d.multiplicity[p] == 0
                            for p in self.d_side_order[i+1:]]):
                        mult_at_i = coeff_d.multiplicity[self.d_side_order[i]]
                        if mult_at_i == 0:
                            arrows_base.append(arrow)
                        else:
                            assert mult_at_i == 1
                            arrows_seed.append(arrow)
            arrows_new = self._autoCompleteByLinAlg(arrows_base, arrows_seed)

            for coeff_d, coeffs_a, gen_from, gen_to in arrows_new:
                self.raw_da.addDelta(gen_from, gen_to, coeff_d, coeffs_a, 1)

        assert self.raw_da.testDelta()
        return self.raw_da
