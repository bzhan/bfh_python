"""Auto-completion of arrows in a type DA structure, by solving certain
equations in linear algebra.

This module is used to produce arrows in the local type DA structure for
arcslides, in arcslidedatest.py.

"""

from algebra import E0
from dastructure import DAStructure
from linalg import F2RowSystem
from localpmc import LocalStrandDiagram
import itertools
from Queue import Queue

class AutoCompleteDAStructure:
    """A routine to complete arrows in type DA structures by solving certain
    equations in linear algebra.

    """
    @staticmethod
    def getDerivedTwoStepArrows(arrows_base, arrow_new, mod_gens):
        """Find all ways of deriving two-step arrows from arrow_new, including
        by combining it with elements of arrows_base.

        """
        result = []
        coeff_d, coeffs_a = arrow_new
        if len(coeffs_a) > 0 and coeffs_a[0].isIdempotent():
            return result
        # Take anti-differential of one of coeffs_a
        for i in range(len(coeffs_a)):
            for anti_da, coeff in coeffs_a[i].antiDiff().items():
                result.append(
                    (coeff_d, coeffs_a[:i] + (anti_da,) + coeffs_a[i+1:]))
            for (a, b), coeff in coeffs_a[i].factor().items():
                if a.isIdempotent() or b.isIdempotent():
                    continue
                result.append((coeff_d, coeffs_a[:i] + (a, b) + coeffs_a[i+1:]))
        # Take differential of coeff_d
        for dd, coeff in coeff_d.diff().items():
            result.append((dd, coeffs_a))
        # Multiply two together
        for coeff_d2, coeffs_a2 in arrows_base:
            if len(coeffs_a2) > 0 and coeffs_a2[0].isIdempotent():
                continue
            # One direction
            if coeff_d * coeff_d2 != E0:
                if any([DAStructure.idemMatchDA(x, y, coeff_d, coeffs_a) and
                        DAStructure.idemMatchDA(y, z, coeff_d2, coeffs_a2)
                        for x, y, z in itertools.product(
                                mod_gens, mod_gens, mod_gens)]):
                    result.append(((coeff_d * coeff_d2).getElt(),
                                   coeffs_a + coeffs_a2))
            # Other direction
            if coeff_d2 * coeff_d != E0:
                if any([DAStructure.idemMatchDA(x, y, coeff_d2, coeffs_a2) and
                        DAStructure.idemMatchDA(y, z, coeff_d, coeffs_a)
                        for x, y, z in itertools.product(
                                mod_gens, mod_gens, mod_gens)]):
                    result.append(((coeff_d2 * coeff_d).getElt(),
                                   coeffs_a2 + coeffs_a))
        return result

    @staticmethod
    def getAltFactorizations(arrows_base, arrow_new, mod_gens):
        result = []
        coeff_d, coeffs_a = arrow_new
        # Take differential of one of coeffs_a
        for i in range(len(coeffs_a)):
            for da, coeff in coeffs_a[i].diff().items():
                result.append((coeff_d, coeffs_a[:i] + (da,) + coeffs_a[i+1:]))
            if i > 0 and coeffs_a[i-1] * coeffs_a[i] != E0:
                result.append((coeff_d, coeffs_a[:i-1] +
                               ((coeffs_a[i-1] * coeffs_a[i]).getElt(),) +
                               coeffs_a[i+1:]))
        # Take anti-differential of coeff_d
        for anti_dd, coeff in coeff_d.antiDiff().items():
            result.append((anti_dd, coeffs_a))
        # Split into two sequences
        for (a, b), coeff in coeff_d.factor().items():
            # Number of A-inputs in the first sequence
            for i in range(len(coeffs_a)):
                if (a, coeffs_a[:i]) in arrows_base:
                    arrow_used = False
                    for x, y in itertools.product(mod_gens, mod_gens):
                        if DAStructure.idemMatchDA(x, y, a, coeffs_a[:i]):
                            arrow_used = True
                    if arrow_used:
                        result.append((b, coeffs_a[i:]))
                if (b, coeffs_a[i:]) in arrows_base:
                    result.append((a, coeffs_a[:i]))
        return result

    @staticmethod
    def getAltIdempotents(arrows, single_idems):
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
        for idem_d, idem_a in single_idems:
            for coeff_d, coeffs_a in results + arrows:
                if has_singlehor(coeff_d, idem_d) and \
                   all([has_singlehor(coeff, idem_a) for coeff in coeffs_a]):
                    new_arrow = (remove_singlehor(coeff_d, idem_d),
                                 tuple([remove_singlehor(coeff, idem_a)
                                        for coeff in coeffs_a]))
                    if not new_arrow in results + arrows:
                        results.append(new_arrow)
                if (not uses_idempotent(coeff_d, idem_d)) and \
                   all([not uses_idempotent(coeff, idem_a)
                        for coeff in coeffs_a]):
                    new_arrow = (add_singlehor(coeff_d, idem_d),
                                 tuple([add_singlehor(coeff, idem_a)
                                        for coeff in coeffs_a]))
                    if not new_arrow in results + arrows:
                        results.append(new_arrow)
        return results

    @staticmethod
    def autoCompleteByLinAlg(arrows_base, arrows_new, single_idems, mod_gens):
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
                for arrow in AutoCompleteDAStructure.getDerivedTwoStepArrows(
                        arrows_base, cur_arrow, mod_gens):
                    if arrow not in two_step_arrows:
                        two_step_arrows.add(arrow)
                        two_step_arrows_queue.put(arrow)
                for arrow in AutoCompleteDAStructure.getAltIdempotents(
                        [cur_arrow], single_idems):
                    if arrow not in one_step_arrows:
                        one_step_arrows.add(arrow)
                        one_step_arrows_queue.put(arrow)
            else:
                cur_arrow = two_step_arrows_queue.get()
                for arrow in AutoCompleteDAStructure.getAltFactorizations(
                        arrows_base, cur_arrow, mod_gens):
                    if arrow not in one_step_arrows:
                        one_step_arrows.add(arrow)
                        one_step_arrows_queue.put(arrow)

        # Combine some of the one-step arrows
        combined_one_step_arrows = []  # list of lists of arrows
        while len(one_step_arrows) != 0:
            arrow = one_step_arrows.pop()
            coeff_d, coeffs_a = arrow
            valid = True
            for idem_d, idem_a in single_idems:
                if idem_a in coeff_d.single_hor and \
                   not all(
                       idem_d in coeff_a.single_hor for coeff_a in coeffs_a):
                    valid = False
                    break
            if not valid:
                # These arrows should not be considered.
                continue

            alt_idems = AutoCompleteDAStructure.getAltIdempotents(
                [arrow], single_idems)
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
                derived_two_steps = AutoCompleteDAStructure\
                    .getDerivedTwoStepArrows(
                        arrows_base, one_step_arrow, mod_gens)
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
