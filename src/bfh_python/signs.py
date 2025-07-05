"""Sign conventions."""

from fractions import Fraction

from .algebra import findRankOverF2
from .algebra import DGAlgebra, Element, Generator
from .algebra import E0
from .grading import standardRefinement, standardRefinementForIdem
from .grading import DEFAULT_REFINEMENT
from .linalg import F2RowSystem
from .pmc import StrandAlgebra, StrandDiagram
from .utility import memorize
from .utility import F2

class AbsZ2Grading(object):
    """Computes absolute Z/2Z grading for the strand algebra."""
    def __init__(self, algebra):
        """Creates data needed for computation of Z/2Z grading for the given
        strand algebra.

        """
        self.algebra = algebra
        gen_list = algebra.getGenerators()
        self.init_idem = gen_list[0].left_idem

        # Compute the set of idempotents that have offset 1/2 with init_idem.
        # This is stored as self.adjusted_idem.
        self.adjusted_idem = set()
        self.tested_idem = set([self.init_idem])

        def testAdd(gen, tested_idem, to_test_idem):
            # Compute absolute grading on gen to see if its two idempotents
            # should have 1/2 offsets.
            if to_test_idem in self.tested_idem:
                return
            self.tested_idem.add(to_test_idem)
            maslov = self._getAbsGradingRaw(gen)
            assert maslov.denominator <= 2

            is_different = (maslov.denominator != 1)
            prev_adjust = (tested_idem in self.adjusted_idem)
            if (is_different and not prev_adjust) or \
                    (not is_different and prev_adjust):
                self.adjusted_idem.add(to_test_idem)

        need_repeat = True
        while need_repeat:
            need_repeat = False
            for gen in gen_list:
                if gen.left_idem in self.tested_idem:
                    testAdd(gen, gen.left_idem, gen.right_idem)
                elif gen.right_idem in self.tested_idem:
                    testAdd(gen, gen.right_idem, gen.left_idem)
                else:
                    need_repeat = True
            
    def _getAbsGradingRaw(self, sd):
        """Returns the maslov grading from which the absolute grading is
        derived (no 1/2 offsets applied).

        """
        assert sd.parent == self.algebra
        small_gr = sd.getSmallGrading(refinement = standardRefinement);
        gr_group = small_gr.parent
        for i in range(len(small_gr.spinc)):
            mult = -small_gr.spinc[i]
            small_gr = small_gr * gr_group.basis(i).power(mult)
            small_gr.maslov += Fraction(-1, 2) * mult
        assert all([n == 0 for n in small_gr.spinc])
        return small_gr.maslov

    def getAbsGrading(self, sd):
        """Returns the absolute Z/2Z grading (returned value is either 0 or 1).

        """
        maslov = self._getAbsGradingRaw(sd)
        if sd.left_idem in self.adjusted_idem:
            maslov -= Fraction(1, 2)
        if sd.right_idem in self.adjusted_idem:
            maslov += Fraction(1, 2)
        assert maslov.denominator == 1
        return maslov.numerator % 2

class PreStrandDiagram(Generator):
    """A pre-strand diagram on n points.

    Difference with StrandDiagram in pmc.py is that there are no matching on
    points. There are only single horizontals, and there may be strands
    beginning or ending at any pair of points.

    """
    def __init__(self, parent, strands):
        """Specifies the parent algebra (of type PreStrandAlgebra), and a list
        of strands. Each element of strands is a pair specifying the starting
        and ending points.

        """
        Generator.__init__(self, parent)
        self.pmc = parent.pmc
        self.strands = tuple(sorted(strands))
        self.left_pt_idem = tuple(sorted([s for s, t in self.strands]))
        self.right_pt_idem = tuple(sorted([t for s, t in self.strands]))
        self.left_idem = tuple(sorted([self.pmc.pairid[s]
                                       for s in self.left_pt_idem]))
        self.right_idem = tuple(sorted([self.pmc.pairid[s]
                                        for s in self.right_pt_idem]))

    def numCrossing(self):
        """Returns the number of crossings between moving strands."""
        # Note this counts crossing between single horizontals
        return sum(1 for (s1, t1) in self.strands for (s2, t2) in self.strands
                   if s1 < s2 and t1 > t2)

    def getBigGrading(self):
        """Returns the big grading."""
        multiplicity = [0] * (self.pmc.n - 1)
        for s, t in self.strands:
            for pos in range(s, t):
                multiplicity[pos] += 1
        maslov = 0
        for s, t in self.strands:
            if s != self.pmc.n - 1:
                maslov -= Fraction(multiplicity[s], 2)
            if s != 0:
                maslov -= Fraction(multiplicity[s-1], 2)
        maslov += self.numCrossing()
        return self.pmc.big_gr(maslov, multiplicity)

    @memorize
    def getSmallGrading(self, refinement = DEFAULT_REFINEMENT):
        """Returns the small grading."""
        # standardRefinement ensures integrality of spin-c components of the
        # grading.
        assert refinement == standardRefinement
        refine_data = refinement(self.pmc, len(self.left_idem))
        p_l, p_r = [standardRefinementForIdem(self.pmc, idem)
                    for idem in [self.left_idem, self.right_idem]]
        return (p_l * self.getBigGrading() * p_r.inverse()).toSmallGrading()

    def isIdempotent(self):
        """Tests whether this generator is an idempotent."""
        return all([s == t for s, t in self.strands])

    def __str__(self):
        return "[%s]" % ",".join(["%s->%s" % (p, q) for p, q in self.strands])

    def __repr__(self):
        return str(self)

    def __eq__(self, other):
        if other is None:
            return False
        return self.parent == other.parent and self.strands == other.strands

    def __ne__(self, other):
        return not (self == other)

    def __hash__(self):
        return hash((self.parent, self.strands))

# Don't need anything beyond Element at this time.
PreStrandDiagram.ELT_CLASS = Element

class PreStrandAlgebra(DGAlgebra):
    """Represents the strand algebra of a local PMC."""

    def __init__(self, ring, pmc, idem_size):
        """Specifies the number of points n and the number of strands k."""
        DGAlgebra.__init__(self, ring)
        assert idem_size <= pmc.genus * 2, "idem_size too large"
        self.pmc = pmc
        self.idem_size = idem_size
        self.abs_gr = AbsZ2Grading(self)

    def __str__(self):
        return "Pre strand algebra over %s with idem_size = %d" % \
            (str(self.pmc), self.idem_size)

    def __eq__(self, other):
        return self.pmc == other.pmc and self.idem_size == other.idem_size

    def __ne__(self, other):
        return not (self == other)

    def __hash__(self):
        return hash(("PreStrandAlgebra", self.pmc, self.idem_size))

    @memorize
    def getGenerators(self):
        result = []
        def helper(strands):
            # strands is the current list of strands.
            if len(strands) == self.idem_size:
                result.append(PreStrandDiagram(self, strands))
                return
            start = 0
            if len(strands) > 0:
                start = strands[-1][0] + 1
            for p in range(start, self.pmc.n):
                for q in range(p, self.pmc.n):
                    if q not in [end for start, end in strands]:
                        helper(strands + [(p, q)])
        helper([])
        return result

    @memorize
    def getGeneratorsForPtIdem(self, l_pt_idem = None, r_pt_idem = None):
        """Returns the list of idempotents with the specified left_pt_idem and
        right_pt_idem. Giving None as input means no constraints there.

        """
        return [gen for gen in self.getGenerators() if
                (l_pt_idem is None or gen.left_pt_idem == l_pt_idem) and
                (r_pt_idem is None or gen.right_pt_idem == r_pt_idem)]

    @memorize
    def _diffRaw(self, gen):
        """Returns a list of elements of the form ((n1, n2), diff_term), where
        n1, n2 are indices of strands in gen that crosses, and diff_term is a
        generator in gen.diff() obtained by uncrossing these two strands.
        Together they specify all terms in gen.diff(). Elements in the list are
        sorted by (n1, n2).

        """
        target_num_crossing = gen.numCrossing() - 1
        result = []
        for n1 in range(len(gen.strands)):
            for n2 in range(len(gen.strands)):
                s1, t1 = gen.strands[n1]
                s2, t2 = gen.strands[n2]
                if s1 < s2 and t1 > t2:
                    new_strands = list(gen.strands)
                    new_strands.remove((s1, t1))
                    new_strands.remove((s2, t2))
                    new_strands.extend([(s1, t2), (s2, t1)])
                    diff_term = PreStrandDiagram(self, new_strands)
                    if diff_term.numCrossing() == target_num_crossing:
                        result.append(((n1, n2), diff_term))
        return result

    @memorize
    def diff(self, gen):
        result = E0
        for (n1, n2), diff_term in self._diffRaw(gen):
            if self.ring is F2:
                result += diff_term.elt()
            else:
                result += self.diffSign(gen, n1, n2) * diff_term
        return result

    @memorize
    def _multiplyRaw(self, gen1, gen2):
        """If gen1 and gen2 can be multiplied, return the generator that is
        their product. Otherwise, return None.

        """
        if gen1.right_pt_idem != gen2.left_pt_idem:
            return None
        new_strands = []
        for s1, t1 in gen1.strands:
            t2_list = [t2 for s2, t2 in gen2.strands if s2 == t1]
            assert len(t2_list) == 1
            new_strands.append((s1, t2_list[0]))
        product = PreStrandDiagram(self, new_strands)
        if product.numCrossing() == gen1.numCrossing() + gen2.numCrossing():
            return product
        else:
            return None

    @memorize
    def multiply(self, gen1, gen2):
        if not isinstance(gen1, PreStrandDiagram):
            return NotImplemented
        if not isinstance(gen2, PreStrandDiagram):
            return NotImplemented
        product = self._multiplyRaw(gen1, gen2)
        if product is None:
            return E0

        if self.ring is F2:
            return product.elt()
        else:
            return self.multiplySign(gen1, gen2) * product

    @memorize
    def diffSign(self, gen, n1, n2):
        """Returns the sign (+/-1) of the differential of gen when resolving the
        crossing between strands indexed n1 and n2.

        The cancellation in d^2=0 we use is:
        gen -- (n1, n2) --> n_term -- (nk1, nk2) --> nk_term
        gen -- (k1, k2) --> k_term -- (kn1, kn2) --> kn_term

        """
        dgen_raw = self._diffRaw(gen)
        n_term = None
        for (m1, m2), diff_term in dgen_raw:
            if (m1, m2) == (n1, n2):
                n_term = diff_term
                break
        assert n_term is not None

        # Key pair (k1, k2) satisfy the property that there are no k' such that
        # (k1, k') is resolvable (by ordering) or (k', k2) is resolvable.
        (k1, k2), k_term = dgen_raw[0]
        for (m1, m2), m_term in dgen_raw:
            if k1 < m1 and k2 == m2:
                (k1, k2), k_term = (m1, m2), m_term

        if (n1, n2) == (k1, k2):
            return 1  # This arrow is defined to be positive

        # Find a cancellation term
        dn_raw = self._diffRaw(n_term)
        dk_raw = self._diffRaw(k_term)
        for (nk1, nk2), nk_term in dn_raw:
            for (kn1, kn2), kn_term in dk_raw:
                if nk_term == kn_term:
                    return self.diffSign(n_term, nk1, nk2) * \
                        self.diffSign(k_term, kn1, kn2) * -1

        # Should not get to this point
        assert False

    @memorize
    def multiplySign(self, a, b):
        """Returns the sign (+/-1) of the multiplication between a and b.
        
        """
        prod = self._multiplyRaw(a, b)
        assert prod is not None

        # All multiplication involving idempotents are positive.
        if a.isIdempotent() or b.isIdempotent():
            return 1

        # prod is differentiable. Define using product on generators with a
        # smaller number of crossings.
        dprod_raw = self._diffRaw(prod)
        if len(dprod_raw) > 0:
            da_raw = self._diffRaw(a)
            db_raw = self._diffRaw(b)
            # Pick a term in d(prod), and find the term in either da*b or a*db
            # that cancels it.
            (p1, p2), dprod_term = dprod_raw[0]
            for (a1, a2), da_term in da_raw:
                if self._multiplyRaw(da_term, b) == dprod_term:
                    return self.diffSign(prod, p1, p2) * \
                        self.diffSign(a, a1, a2) * \
                        self.multiplySign(da_term, b)
            for (b1, b2), db_term in db_raw:
                if self._multiplyRaw(a, db_term) == dprod_term:
                    return self.diffSign(prod, p1, p2) * \
                        self.diffSign(b, b1, b2) * \
                        self.multiplySign(a, db_term) * self.grSign(a)
            # Should not get to this point
            assert False

        # Now prod (and therefore a and b) has no crossings. (k_s, k_t) is the
        # initial piece of the last moving strand in prod (this ensures that the
        # chord can always be factored from the left).
        prod_s, prod_t = next((s, t) for s, t in reversed(prod.strands)
                              if s < t)
        k_s, k_t = (prod_s, prod_s + 1)
        other_chords_prod = [(s, t) for s, t in prod.strands if s != prod_s]
        k_a = PreStrandDiagram(
            self, [(s, s) for s, t in other_chords_prod] + [(k_s, k_t)])
        k_b = PreStrandDiagram(self, other_chords_prod + [(k_t, prod_t)])
        assert self._multiplyRaw(k_a, k_b) == prod

        if a == k_a:
            return 1  # This product is defined to be positive

        # Two cases for when prod has no crossings
        ka_s, ka_t = next((s, t) for s, t in a.strands if s == k_s)
        if ka_t > ka_s:  # Case 1: (k_s, k_t) is part of a
            other_chords_a = [(s, t) for s, t in a.strands if s != k_s]
            l_a = PreStrandDiagram(self, other_chords_a + [(k_t, ka_t)])
            assert self._multiplyRaw(k_a, l_a) == a
            assert self._multiplyRaw(k_a, self._multiplyRaw(l_a, b)) == prod
            return self.multiplySign(k_a, l_a) * self.multiplySign(l_a, b)
        else:  # Case 2: (k_s, k_t) is not part of a, ka_s == k_s == ka_t
            other_chords_b = [(s, t) for s, t in b.strands if s != k_s]
            kb_s, kb_t = prod_s, prod_t
            assert (kb_s, kb_t) in b.strands
            l1_b = PreStrandDiagram(
                self, [(s, s) for s, t in other_chords_b] + [(k_s, k_t)])
            l2_b = PreStrandDiagram(self, other_chords_b + [(k_t, prod_t)])
            if l2_b.isIdempotent():
                # Will cause an infinite loop. This reflects the assumption that
                # a * l1_b = k_a * k_b with positive sign.
                return 1
            assert self._multiplyRaw(l1_b, l2_b) == b
            assert self._multiplyRaw(self._multiplyRaw(a, l1_b), l2_b) == prod
            return self.multiplySign(l1_b, l2_b) * \
                self.multiplySign(a, l1_b) * \
                self.multiplySign(self._multiplyRaw(a, l1_b), l2_b)

class SignLinAlg(object):
    """Try to obtain a sign convention by solving a linear algebra problem. """
    def __init__(self, algebra):
        self.algebra = algebra
        self.abs_gr = AbsZ2Grading(algebra)

    def grSign(self, gen):
        """Returns (-1)^gr(gen), where gr is the absolute Z/2Z grading."""
        return 1 - 2*self.abs_gr.getAbsGrading(gen)

    def createRowSystem(self):
        """Create row system. Each row represents either a multiplication or an
        arrow in the differential. Each column represents one of the constraints
        that must be satisfied.

        """
        all_gens = [gen for gen in self.algebra.getGenerators()
                    if not gen.isIdempotent()]
        print("Number of generators:", len(all_gens))
        # Maps multiplication / differential to the column index
        self.index = dict()
        for gen1 in all_gens:
            for gen2 in self.algebra.getGeneratorsForIdem(
                    left_idem = gen1.right_idem):
                if (not gen2.isIdempotent()) and gen1 * gen2 != 0:
                    self.index[("M", gen1, gen2)] = len(self.index)
        for gen in all_gens:
            for dgen_term in gen.diff():
                self.index[("D", gen, dgen_term)] = len(self.index)
        num_row = len(self.index)
        print("Number of operations:", num_row)

        # Linear combination of rows to look for.
        expected_sums = []
        # Pairs (i, j), zero-based indices indicating a_ij is one.
        entries = []
        # Current number of columns
        num_col = 0

        def addColumn(ops, expected_sum):
            for op in ops:
                entries.append((self.index[op], num_col))
            expected_sums.append(expected_sum)

        # Now create row for each relation
        for gen in all_gens:
            # First part: d^2 = 0
            # Create a map from terms in ddgen to list of terms in dgen
            dd_to_d_map = dict()
            for dgen in gen.diff():
                for ddgen in dgen.diff():
                    if ddgen not in dd_to_d_map:
                        dd_to_d_map[ddgen] = []
                    dd_to_d_map[ddgen].append(dgen)
            # Now use that map to produce the relations in d^2=0
            for ddgen, dgen_list in list(dd_to_d_map.items()):
                assert len(dgen_list) == 2
                dgen1, dgen2 = dgen_list
                addColumn([("D", gen, dgen1), ("D", gen, dgen2),
                           ("D", dgen1, ddgen), ("D", dgen2, ddgen)], 1)
                num_col += 1

        for gen1 in all_gens:
            for gen2 in self.algebra.getGeneratorsForIdem(
                    left_idem = gen1.right_idem):
                if gen2.isIdempotent() or gen1 * gen2 == 0:
                    continue
                # Second part: d(ab) = da * b + (-1)^gr(a) a*db
                if (gen1 * gen2).diff() != 0:
                    for dgen12 in (gen1 * gen2).diff():
                        for dgen1 in gen1.diff():
                            if dgen1 * gen2 == dgen12.elt():
                                addColumn([("M", gen1, gen2),
                                           ("D", (gen1*gen2).getElt(), dgen12),
                                           ("D", gen1, dgen1),
                                           ("M", dgen1, gen2)], 0)
                                num_col += 1
                        for dgen2 in gen2.diff():
                            if gen1 * dgen2 == dgen12.elt():
                                addColumn([("M", gen1, gen2),
                                           ("D", (gen1*gen2).getElt(), dgen12),
                                           ("D", gen2, dgen2),
                                           ("M", gen1, dgen2)],
                                          self.abs_gr.getAbsGrading(gen1))
                                num_col += 1
                # Third part: (ab)c = a(bc)
                for gen3 in self.algebra.getGeneratorsForIdem(
                        left_idem = gen2.right_idem):
                    if gen3.isIdempotent() or gen2 * gen3 == 0:
                        continue
                    if (gen1*gen2) * gen3 == 0:
                        continue
                    addColumn(
                        [("M", gen1, gen2), ("M", (gen1*gen2).getElt(), gen3),
                         ("M", gen2, gen3), ("M", gen1, (gen2*gen3).getElt())],
                        0)
                    num_col += 1

        print("Number of constraints:", num_col)

        # Use the following to find a solution
        # matrix = [[0] * num_col for i in range(num_row)]
        # for i, j in entries:
        #     matrix[i][j] = 1
        # row_sys = F2RowSystem(matrix)
        # comb = row_sys.getComb(expected_sums)
        # assert comb is not None, "Cannot be solved"
        # for op, index in self.index.items():
        #     print op, comb[index]

        # Use simplify method
        row_rank = findRankOverF2(num_row, num_col, entries)
        print("Rank:", row_rank)

        # Form row system of gauge equivalences
        gen_index = dict()
        for i in range(len(all_gens)):
            gen_index[all_gens[i]] = i
        # Pairs (i, j), zero-based indices indicating a_ij is one in gauge
        # matrix.
        gauge_entries = []

        for op, index in list(self.index.items()):
            assert op[1] != op[2]
            gauge_entries.append((index, gen_index[op[1]]))
            gauge_entries.append((index, gen_index[op[2]]))
            if op[0] == "M":
                gauge_entries.append((index, gen_index[(op[1]*op[2]).getElt()]))

        gauge_rank = findRankOverF2(num_row, len(all_gens), gauge_entries)
        print("Rank of gauge equivalences:", gauge_rank)
        print("Free choice:", num_row - row_rank - gauge_rank)
