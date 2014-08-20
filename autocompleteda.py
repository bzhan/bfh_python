"""Auto-completion of arrows in a type DA structure, by solving certain
equations in linear algebra.

This module is used to produce arrows in the local type DA structure for
arcslides, in arcslidedatest.py.

"""

from algebra import solveOverF2
from algebra import E0
from dastructure import DAStructure, MorDAtoDAGenerator
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

class _AutoCompleteDAStructure:
    """A routine to complete arrows in type DA structures by solving certain
    equations in linear algebra.

    """
    def _getDerivedTwoStepArrows(self, arrows_base_left, arrows_base_right,
                                 arrow_new):
        """Find all ways of deriving two-step arrows from arrow_new, including
        by combining it with elements of arrows_base_left and arrows_base_right.

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
        # Multiply two together. One direction.
        for coeff_d2, coeffs_a2, x2, y2 in arrows_base_left:
            if len(coeffs_a2) > 0 and coeffs_a2[0].isIdempotent():
                continue
            if y2 == x and coeff_d2 * coeff_d != E0:
                result.append(_DAArrow(
                    (coeff_d2 * coeff_d).getElt(), coeffs_a2 + coeffs_a, x2, y))
        # The other direction.
        for coeff_d2, coeffs_a2, x2, y2 in arrows_base_right:
            if len(coeffs_a2) > 0 and coeffs_a2[0].isIdempotent():
                continue
            if y == x2 and coeff_d * coeff_d2 != E0:
                result.append(_DAArrow(
                    (coeff_d * coeff_d2).getElt(), coeffs_a + coeffs_a2, x, y2))
        return result

    def _getAltFactorizations(self, arrows_base_left_map, arrows_base_right_map,
                              arrow_new):
        """Find all one-step-arrows that can produce the two-step-arrow
        arrow_new. Here arrows_base_left and arrows_base_right are each passed
        as two maps. Each arrow (coeff_d, coeffs_a, x, z) in arrows_base_left is
        indexed in arrows_base_left_map as (coeff_d, coeffs_a, x) -> z, and each
        arrow (coeff_d, coeffs_a, z, y) in arrows_base_right is indexed in
        arrows_base_right_map as (coeff_d, coeffs_a, y) -> z.

        """
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
                if (a, coeffs_a[:i], x) in arrows_base_left_map:
                    for z in arrows_base_left_map[(a, coeffs_a[:i], x)]:
                        result.append(_DAArrow(b, coeffs_a[i:], z, y))
                if (b, coeffs_a[i:], y) in arrows_base_right_map:
                    for z in arrows_base_right_map[(b, coeffs_a[i:], y)]:
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
                if x in self.da_left.u_maps[i] and \
                   y in self.da_right.u_maps[i] and \
                   has_singlehor(coeff_d, idem_d) and \
                   all([has_singlehor(coeff, idem_a)
                        for coeff in coeffs_a]):
                    new_arrow = _DAArrow(
                        remove_singlehor(coeff_d, idem_d),
                        tuple([remove_singlehor(coeff, idem_a)
                               for coeff in coeffs_a]),
                        self.da_left.u_maps[i][x],
                        self.da_right.u_maps[i][y])
                    if not new_arrow in results + arrows:
                        results.append(new_arrow)
                if (x in self.da_left.uinv_maps[i] and \
                    y in self.da_right.uinv_maps[i] and \
                    (not uses_idempotent(coeff_d, idem_d)) and \
                    all([not uses_idempotent(coeff, idem_a)
                         for coeff in coeffs_a])):
                    new_arrow = _DAArrow(
                        add_singlehor(coeff_d, idem_d),
                        tuple([add_singlehor(coeff, idem_a)
                               for coeff in coeffs_a]),
                        self.da_left.uinv_maps[i][x],
                        self.da_right.uinv_maps[i][y])
                    if not new_arrow in results + arrows:
                        results.append(new_arrow)
        return results

    def _autoCompleteByLinAlg(self, arrows_base_left, arrows_base_right,
                              arrows_new):
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
        # Form the arrow base maps for faster computation in
        # _getAltFactorizations.
        arrows_base_left_map = dict()
        arrows_base_right_map = dict()
        for coeff_d, coeffs_a, source, target in arrows_base_left:
            key = (coeff_d, coeffs_a, source)
            if key not in arrows_base_left_map:
                arrows_base_left_map[key] = []
            arrows_base_left_map[key].append(target)
        for coeff_d, coeffs_a, source, target in arrows_base_right:
            key = (coeff_d, coeffs_a, target)
            if key not in arrows_base_right_map:
                arrows_base_right_map[key] = []
            arrows_base_right_map[key].append(source)

        while not one_step_arrows_queue.empty() or \
              not two_step_arrows_queue.empty():
            if not one_step_arrows_queue.empty():
                cur_arrow = one_step_arrows_queue.get()
                for arrow in self._getDerivedTwoStepArrows(
                        arrows_base_left, arrows_base_right, cur_arrow):
                    if arrow not in two_step_arrows:
                        two_step_arrows.add(arrow)
                        two_step_arrows_queue.put(arrow)
                for arrow in self._getAltIdempotents([cur_arrow]):
                    if arrow not in one_step_arrows:
                        one_step_arrows.add(arrow)
                        one_step_arrows_queue.put(arrow)
            else:
                cur_arrow = two_step_arrows_queue.get()
                for arrow in self._getAltFactorizations(
                        arrows_base_left_map, arrows_base_right_map, cur_arrow):
                    coeff_d, coeffs_a, x, y = arrow
                    # HACK - it appears considering one_step_arrows with at most
                    # four algebra inputs is sufficient. In the autocompletion
                    # for the anti-braid case, there are infinitely many
                    # reachable one_step_arrows, so we place this limit.
                    # Can change to 2 or 3 to see if simpler arrows are
                    # possible.
                    if len(coeffs_a) > 4:
                        continue
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
        num_row, num_col = len(combined_one_step_arrows), len(two_step_arrows)
        matrix_entries = set()
        target_vec = set()

        two_step_arrows_dict = dict()  # index the two step arrows.
        two_step_arrows = list(two_step_arrows)
        for i in range(len(two_step_arrows)):
            two_step_arrows_dict[two_step_arrows[i]] = i

        for i in range(len(combined_one_step_arrows)):
            for one_step_arrow in combined_one_step_arrows[i]:
                derived_two_steps = self._getDerivedTwoStepArrows(
                    arrows_base_left, arrows_base_right, one_step_arrow)
                for two_step_arrow in derived_two_steps:
                    j = two_step_arrows_dict[two_step_arrow]
                    if one_step_arrow in arrows_new:
                        if j in target_vec:
                            target_vec.remove(j)
                        else:
                            target_vec.add(j)
                    else:
                        if (i, j) in matrix_entries:
                            matrix_entries.remove((i, j))
                        else:
                            matrix_entries.add((i, j))

        comb = solveOverF2(num_row, num_col, list(matrix_entries),
                           list(target_vec))
        assert comb is not None

        result = []
        for term in comb:
            if combined_one_step_arrows[term][0] not in arrows_new:
                result.extend(combined_one_step_arrows[term])
        return result

    def _arrows_to_string(self, arrow_set):
        output_strs = set()  # remove duplicates
        for coeff_d, coeffs_a, gen_from, gen_to in arrow_set:
            has_class_from = isinstance(gen_from.name, str) and \
                             '_' in gen_from.name
            has_class_to = isinstance(gen_to.name, str) and '_' in gen_to.name
            if has_class_from or has_class_to:
                n1 = gen_from.name[0:1] if has_class_from else "-1"
                n2 = gen_to.name[0:1] if has_class_to else "-1"
                output_strs.add("(%s, %s, %s)," % (n1, n2, ", ".join(
                    coeff.inputForm()
                    for coeff in list(coeffs_a) + [coeff_d])))
            else:
                output_strs.add("(%s)," % ", ".join(
                    coeff.inputForm()
                    for coeff in list(coeffs_a) + [coeff_d]))
        return "\n".join(output_strs)

    def complete(self, raw_da, d_side_order):
        """Input raw_da is a local type DA structure with all generators and
        some arrows filled in. This function adds arrows to raw_da so that it
        satisfies the type DA structure equation, as well as the smearability
        condition.

        """
        # Initialize the needed values.
        assert isinstance(raw_da, LocalDAStructure)
        self.da_left = self.da_right = raw_da

        # Prepare single idems in raw_da as a list of pairs (idem_d, idem_a).
        self.single_idems = []
        for i in range(raw_da.num_single_idems):
            self.single_idems.append(
                (raw_da.single_idems1[i], raw_da.single_idems2[i]))

        # Complete in stages, each adding one interval on the D side, according
        # to d_side_order.
        for i in range(len(d_side_order)):
            # Lists of _DAArrow objects. In this case, arrows_base_left and
            # arrows_base_right are the same.
            arrows_base, arrows_seed = [], []
            # Figure out the base arrows and new arrows in this stage.
            for (gen_from, coeffs_a), target in raw_da.da_action.items():
                for (coeff_d, gen_to), ring_coeff in target.items():
                    arrow = _DAArrow(coeff_d, coeffs_a, gen_from, gen_to)
                    if all([coeff_d.multiplicity[p] == 0
                            for p in d_side_order[i+1:]]):
                        mult_at_i = coeff_d.multiplicity[d_side_order[i]]
                        if mult_at_i == 0:
                            arrows_base.append(arrow)
                        else:
                            assert mult_at_i == 1
                            arrows_seed.append(arrow)
            arrows_new = self._autoCompleteByLinAlg(
                arrows_base, arrows_base, arrows_seed)

            for coeff_d, coeffs_a, gen_from, gen_to in arrows_new:
                raw_da.addDelta(gen_from, gen_to, coeff_d, coeffs_a, 1)

            ### Uncomment to see the added arrows.
            # if i == 0:
            #     print "# Initial patterns: %s\n" % \
            #         self._arrows_to_string(arrows_base)
            # print "# Step %d, D side position %d:" % \
            #     (i+1, self.d_side_order[i])
            # print "# Seed arrows: %s\n# New arrows %s\n" % \
            #     (self._arrows_to_string(arrows_seed),
            #      self._arrows_to_string(arrows_new))

        # Final check
        assert raw_da.testDelta()
        return raw_da

    def completeMorphism(self, dastr1, dastr2, raw_morphism):
        """Input two local DA structures dastr1 and dastr2, and an incomplete
        morphism between the two. Add arrows to the morphism so that it
        satisfies the structure equations, as well as the smearability
        condition.

        """
        # Initialize the needed values.
        assert isinstance(dastr1, LocalDAStructure)
        assert isinstance(dastr2, LocalDAStructure)
        assert dastr1.algebra1 == dastr2.algebra1
        assert dastr1.algebra2 == dastr2.algebra2
        assert dastr1.single_idems1 == dastr2.single_idems1
        assert dastr1.single_idems2 == dastr2.single_idems2

        self.da_left, self.da_right = dastr1, dastr2
        self.raw_morphism = raw_morphism

        # Prepare single idems in dastr1 (or dastr2) as a list of pairs
        # (idem_d, idem_a).
        self.single_idems = []
        for i in range(dastr1.num_single_idems):
            self.single_idems.append(
                (dastr1.single_idems1[i], dastr1.single_idems2[i]))

        arrows_base_left, arrows_base_right, arrows_seed = [], [], []
        for (gen_from, coeffs_a), target in dastr1.da_action.items():
            for (coeff_d, gen_to), ring_coeff in target.items():
                arrows_base_left.append(
                    _DAArrow(coeff_d, coeffs_a, gen_from, gen_to))
        for (gen_from, coeffs_a), target in dastr2.da_action.items():
            for (coeff_d, gen_to), ring_coeff in target.items():
                arrows_base_right.append(
                    _DAArrow(coeff_d, coeffs_a, gen_from, gen_to))
        for gen in raw_morphism:
            x, (coeff_d, coeffs_a), y = gen.source, gen.coeff, gen.target
            arrows_seed.append(
                _DAArrow(coeff_d, tuple(coeffs_a), x, y))

        arrows_new = self._autoCompleteByLinAlg(
            arrows_base_left, arrows_base_right, arrows_seed)

        mor_parent = raw_morphism.getElt().parent
        for coeff_d, coeffs_a, gen_from, gen_to in arrows_new:
            raw_morphism += 1 * MorDAtoDAGenerator(
                mor_parent, coeff_d, coeffs_a, gen_from, gen_to)

        ### Uncomment to see the added arrows
        # print "New arrows:\n%s\n" % self._arrows_to_string(arrows_new)

        # Final check
        assert raw_morphism.diff() == 0
        return raw_morphism

def autoCompleteDA(raw_da, d_side_order):
    auto = _AutoCompleteDAStructure()
    auto.complete(raw_da, d_side_order)

def autoCompleteMorphism(dastr1, dastr2, raw_morphism):
    auto = _AutoCompleteDAStructure()
    auto.completeMorphism(dastr1, dastr2, raw_morphism)
