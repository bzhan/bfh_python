"""Arcslides and their type DD structures."""

from ddstructure import DDStrFromChords
from grading import SimpleDbGradingSet, SimpleDbGradingSetElement, \
    SmallGradingGroup
from hdiagram import getArcslideDiagram
from pmc import Idempotent, PMC, Strands
from pmc import connectSumPMC, splitPMC
from utility import memorize
from utility import ACTION_LEFT

# Two types of arcslides:
UNDER_SLIDE, OVER_SLIDE = 0, 1

class Arcslide:
    """Represents an arcslide."""
    def __init__(self, start_pmc, b1, c1):
        """Specifies the starting pmc, the sliding point (b1), and the point it
        slides over (c1).

        """
        assert b1 == c1-1 or b1 == c1+1
        self.start_pmc, self.b1, self.c1 = start_pmc, b1, c1
        self.n = self.start_pmc.n
        self.b_pair = self.start_pmc.pairid[b1]
        self.c_pair = self.start_pmc.pairid[c1]
        self.b2 = start_pmc.otherp[self.b1]
        self.c2 = start_pmc.otherp[self.c1]

        if self.c1 < self.b1 < self.c2 or self.c2 < self.b1 < self.c1:
            self.slide_type = UNDER_SLIDE
        else:
            self.slide_type = OVER_SLIDE

        if self.b1 == self.c1-1:
            if self.c2 > self.c1:
                self.to_r = self._getShiftMap(self.b1, self.c1, self.c2)
            else: # self.c2 < self.c1
                self.to_r = self._getShiftMap(self.b1, self.c2+1, self.c1-2)
        else: # self.b1 == self.c1+1
            if self.c2 > self.c1:
                self.to_r = self._getShiftMap(self.b1, self.c1+2, self.c2-1)
            else: # self.c2 < self.c1
                self.to_r = self._getShiftMap(self.b1, self.c2, self.c1)

        self.end_pmc = PMC([(self.to_r[p], self.to_r[q])
                            for p, q in self.start_pmc.pairs])
        self.pair_to_r = dict()
        for i in range(self.n/2):
            p, q = self.start_pmc.pairs[i]
            self.pair_to_r[i] = self.end_pmc.pairid[self.to_r[p]]

    def _getShiftMap(self, p, q1, q2):
        """Move p to the other side of the interval [q1,q2] (inclusive)."""
        assert p == q1-1 or p == q2+1
        n = self.start_pmc.n
        result = dict()
        for i in range(n):
            if i == p:
                if p == q1-1: result[i] = q2
                else: result[i] = q1
            elif q1 <= i <= q2:
                if p == q1-1: result[i] = i-1
                else: result[i] = i+1
            else:
                result[i] = i
        return result

    def __eq__(self, other):
        return self.start_pmc == other.start_pmc and self.b1 == other.b1 \
            and self.c1 == other.c1

    def __ne__(self, other):
        return not (self == other)

    def __hash__(self):
        return hash((self.start_pmc, self.b1, self.c1, "Arcslide"))

    def __str__(self):
        if self.slide_type == UNDER_SLIDE:
            result = "Underslide of "
        else:
            result = "Overslide of "
        result += "point %d over %d starting at %s" % \
                  (self.b1, self.c1, self.start_pmc)
        return result

    def __repr__(self):
        return str(self)

    def inverse(self):
        """Returns the inverse arcslide."""
        return Arcslide(self.end_pmc, self.to_r[self.b1], self.to_r[self.c2])

    @memorize
    def getDDStructure(self, abs_gr_info = None):
        """Returns the type DD structure corresponding to this arcslide."""
        self.all_idems = self.getIdems()
        self.all_chords = []
        if self.slide_type == UNDER_SLIDE:
            for chord_type in self._UChords:
                chord_type(self)
        else:
            for chord_type in self._OChords:
                chord_type(self)
        alg1 = self.start_pmc.getAlgebra()
        alg2 = self.end_pmc.opp().getAlgebra()
        ddstr = DDStrFromChords(alg1, alg2, self.all_idems, self.all_chords)
        if abs_gr_info is not None:
            self._getAbsGrading(ddstr, abs_gr_info)
        else:
            for gen in ddstr.getGenerators():
                base_gen = gen
                break
            ddstr.registerHDiagram(getArcslideDiagram(self), base_gen)
        return ddstr

    def _getAbsGrading(self, ddstr, abs_gr_info, expand_after = True):
        """Find absolute grading for type DD structure (to replace finding a
        relative grading). It can be shown that value of expand_after does not
        matter.

        """
        # Find grading of elements:
        # Create Heegaard diagram for the expanded arcslide.
        genus = self.start_pmc.genus
        if expand_after:
            ex_start_pmc = connectSumPMC(self.start_pmc, splitPMC(genus))
            ex_slide = Arcslide(ex_start_pmc, self.b1, self.c1)
        else:
            ex_start_pmc = connectSumPMC(splitPMC(genus), self.start_pmc)
            ex_slide = Arcslide(ex_start_pmc, 4*genus+self.b1, 4*genus+self.c1)
        ex_end_pmc = ex_slide.end_pmc
        ex_hdiagram = getArcslideDiagram(ex_slide)
        hgens = ex_hdiagram.getHFGenerators()

        # First step, find grading of two extreme generators, using one of
        # which as the base generator.
        base_idem = [Idempotent(ex_start_pmc, range(2*genus)),
                     Idempotent(ex_end_pmc.opp(), range(2*genus))]
        base_idem2 = [Idempotent(ex_start_pmc, range(2*genus, 4*genus)),
                      Idempotent(ex_end_pmc.opp(), range(2*genus, 4*genus))]
        if abs_gr_info == 0:
            base_gen = ex_hdiagram.getGeneratorByIdem(base_idem, True)
        else:
            base_gen = ex_hdiagram.getGeneratorByIdem(base_idem2, True)
        ex_gr_set, ex_grs = ex_hdiagram.computeDDGrading(base_gen)

        # Second step
        # # print ex_gr_set
        # base_gr1 = ex_grs[base_gen]
        # base_gr2 = ex_grs[base_gen2]
        # # print base_gr1, base_gr2
        # spinc11, spinc12 = base_gr1.data[0].spinc, base_gr1.data[1].spinc
        # spinc21, spinc22 = base_gr2.data[0].spinc, base_gr2.data[1].spinc
        # spincmavg1 = [-(a+b)/Fraction(2) for a,b in zip(spinc11, spinc21)]
        # spincmavg2 = [-(a+b)/Fraction(2) for a,b in zip(spinc12, spinc22)]
        # maslovavg = base_gr1.data[0].maslov + base_gr2.data[0].maslov + \
        #     base_gr1.data[1].maslov + base_gr2.data[1].maslov
        # mavg1 = SmallGradingElement(
        #     base_gr1.data[0].parent, -maslovavg, spincmavg1)
        # mavg2 = SmallGradingElement(
        #     base_gr2.data[0].parent, -Fraction(1,16), spincmavg2)
        # ex_gr_set, ex_grs = ex_hdiagram.computeDDGrading(base_gen,
        #                                                  (mavg1, mavg2))
        # # print ex_gr_set
        # # print ex_grs[base_gen], ex_grs[base_gen2]

        # Form the grading set for the original arcslide from the extended one
        periodic_domains = []
        for ex_gr1, ex_gr2 in ex_gr_set.periodic_domains:
            ex_spinc1, ex_spinc2 = ex_gr1.spinc, ex_gr2.spinc
            if expand_after:
                spinc1, spinc2 = ex_spinc1[0:2*genus], ex_spinc2[2*genus:]
            else:
                spinc1, spinc2 = ex_spinc1[2*genus:], ex_spinc2[0:2*genus]
            if not (all([n == 0 for n in spinc1]) and
                    all([n == 0 for n in spinc2])):
                gr1 = self.start_pmc.small_gr(ex_gr1.maslov, spinc1)
                gr2 = self.end_pmc.opp().small_gr(ex_gr2.maslov, spinc2)
                periodic_domains.append([gr1, gr2])
        ddstr.gr_set = SimpleDbGradingSet(
            SmallGradingGroup(self.start_pmc), ACTION_LEFT,
            SmallGradingGroup(self.end_pmc.opp()), ACTION_LEFT,
            periodic_domains)

        # Now obtain the grading of generators
        ddstr.grading = dict()
        for hgen, ex_gr in ex_grs.items():
            ex_idem1, ex_idem2 = hgen.getDIdem()
            ex_gr1, ex_gr2 = ex_gr.data
            if expand_after:
                pairs1 = [n for n in ex_idem1 if n < 2*genus]
                pairs2 = [n-2*genus for n in ex_idem2 if n >= 2*genus]
            else:
                pairs1 = [n-2*genus for n in ex_idem1 if n >= 2*genus]
                pairs2 = [n for n in ex_idem2 if n < 2*genus]
            if len(pairs1) == genus:
                # Get the corresponding generator in ddstr
                cur_idem = [Idempotent(self.start_pmc, pairs1),
                            Idempotent(self.end_pmc.opp(), pairs2)]
                cur_gen = [gen for gen in ddstr.getGenerators()
                           if [gen.idem1, gen.idem2] == cur_idem]
                assert len(cur_gen) == 1
                cur_gen = cur_gen[0]
                # Finally get the grading
                ex_spinc1, ex_spinc2 = ex_gr1.spinc, ex_gr2.spinc
                if expand_after:
                    for i in range(2*genus):
                        assert ex_spinc2[i] == ex_spinc1[4*genus-i-1]
                    spinc1 = ex_spinc1[0:2*genus]
                    spinc2 = ex_spinc2[2*genus:]
                else:
                    for i in range(2*genus, 4*genus):
                        assert ex_spinc2[i] == ex_spinc1[4*genus-i-1]
                    spinc1 = ex_spinc1[2*genus:]
                    spinc2 = ex_spinc2[0:2*genus]
                gr = SimpleDbGradingSetElement(
                    ddstr.gr_set,
                    [self.start_pmc.small_gr(ex_gr1.maslov, spinc1),
                     self.end_pmc.opp().small_gr(ex_gr2.maslov, spinc2)])
                if cur_gen in ddstr.grading:
                    assert ddstr.grading[cur_gen] == gr
                else:
                    ddstr.grading[cur_gen] = gr
        ddstr.checkGrading()
        return ddstr

    def getIdems(self):
        """Returns the set of possible idempotent-pairs for generators."""
        all_idems = []

        def shift_idem(idem):
            """Find the corresponding idempotent at right."""
            return Idempotent(self.end_pmc, [self.pair_to_r[i] for i in idem])

        left_idems = self.start_pmc.getIdempotents()
        # Generators of type X (complementary)
        for idem in left_idems:
            all_idems.append((idem, shift_idem(idem).opp().comp()))
        # Generators of type Y (sub-complementary)
        for idem in left_idems:
            if self.c_pair in idem and not self.b_pair in idem:
                idem_data = list(shift_idem(idem).comp())
                idem_data.remove(self.pair_to_r[self.b_pair])
                idem_data.append(self.pair_to_r[self.c_pair])
                right_idem = Idempotent(self.end_pmc, idem_data).opp()
                all_idems.append((idem, right_idem))
        return all_idems

    def _addPair(self, data1, data2):
        """Add this pair of chord data."""
        # For left chords, don't need to change anything
        chord_left = Strands(self.start_pmc, data1)
        # For right chords, find corresponding position at right and then take
        # opposite.
        data2 = [(self.to_r[p], self.to_r[q]) for p, q in data2]
        data2 = [(self.n-1-q, self.n-1-p) for p, q in data2]
        chord_right = Strands(self.end_pmc.opp(), data2)
        self.all_chords.append((chord_left, chord_right))

    def _U1Chords(self):
        for x in range(self.n):
            for y in range(x+1, self.n):
                if x != self.b1 and y != self.b1 and \
                        self.start_pmc.otherp[x] != y:
                    self._addPair([(x, y)], [(x, y)])

    def _U2Chords(self):
        b1, c1, c2 = self.b1, self.c1, self.c2
        if b1 < c1: # so to_r[c2] < to_r[b1]
            self._addPair([(b1, c1)], [])
            self._addPair([], [(c2, b1)])
        if b1 > c1: # so to_r[b1] < to_r[c2]
            self._addPair([(c1, b1)], [])
            self._addPair([], [(b1, c2)])

    def _U3Chords(self):
        b1, c1, c2 = self.b1, self.c1, self.c2
        if b1 < c1:
            for x in range(c1+1, self.n):
                self._addPair([(b1, x)], [(c1, x)])
            for x in range(0, c2):
                self._addPair([(x, c2)], [(x, b1)])
        if b1 > c1:
            for x in range(c2+1, self.n):
                self._addPair([(c2, x)], [(b1, x)])
            for x in range(0, c1):
                self._addPair([(x, b1)], [(x, c1)])

    def _U4Chords(self):
        b1, c1, c2 = self.b1, self.c1, self.c2
        # Two connected chords
        if b1 < c1:
            for x in range(0, b1):
                self._addPair([(x, b1)], [(x, c1)])
            for x in range(c2+1, self.n):
                if x != b1:
                    self._addPair([(c2, x)], [(b1, x)])
        if b1 > c1:
            for x in range(b1+1, self.n):
                self._addPair([(b1, x)], [(c1, x)])
            for x in range(0, c2):
                if x != b1:
                    self._addPair([(x, c2)], [(x, b1)])
        # Three connected chords
        if b1 < c1:
            for x in range(0, b1):
                for y in range(c1+1, self.n):
                    self._addPair([(x, b1), (c1, y)], [(x, y)])
            for x in range(0, c2):
                for y in range(c2+1, self.n):
                    if y != b1:
                        self._addPair([(x, y)], [(x, c2), (b1, y)])
        if b1 > c1:
            for x in range(0, c1):
                for y in range(b1+1, self.n):
                    self._addPair([(x, c1), (b1, y)], [(x, y)])
            for x in range(0, c2):
                for y in range(c2+1, self.n):
                    if x != b1:
                        self._addPair([(x, y)], [(x, b1), (c2, y)])

    def _U5Chords(self):
        b1, c1, c2 = self.b1, self.c1, self.c2
        sc, bc = min(c1, c2), max(c1, c2)
        for x in range(0, sc):
            for y in range(bc+1, self.n):
                self._addPair([(x, sc), (bc, y)], [(x, sc), (bc, y)])

    def _U6Chords(self):
        b1, c1, c2 = self.b1, self.c1, self.c2
        if b1 < c1:
            for x in range(c2+1, self.n):
                if x != b1 and x != c1:
                    self._addPair([(c2, x), (b1, c1)], [(b1, x)])
            for x in range(0, b1):
                if x != c2:
                    self._addPair([(x, b1)], [(x, c1), (c2, b1)])
        if b1 > c1:
            for x in range(0, c2):
                if x != b1 and x != c1:
                    self._addPair([(x, c2), (c1, b1)], [(x, b1)])
            for x in range(b1+1, self.n):
                if x != c2:
                    self._addPair([(b1, x)], [(c1, x), (b1, c2)])

    _UChords = [_U1Chords, _U2Chords, _U3Chords, _U4Chords, _U5Chords, \
                    _U6Chords]

    _O1Chords = _U1Chords
    _O2Chords = _U2Chords
    _O6Chords = _U6Chords

    def _O3Chords(self):
        b1, c1, c2 = self.b1, self.c1, self.c2
        if b1 < c1:
            for x in range(c1+1, self.n):
                if x != c2:
                    self._addPair([(b1, x)], [(c1, x)])
            for x in range(0, c2):
                if x != b1 and x != c1:
                    self._addPair([(x, c2)], [(x, b1)])
        if b1 > c1:
            for x in range(0, c1):
                if x != c2:
                    self._addPair([(x, b1)], [(x, c1)])
            for x in range(c2+1, self.n):
                if x != b1 and x != c1:
                    self._addPair([(c2, x)], [(b1, x)])

    def _O4Chords(self):
        b1, c1, c2 = self.b1, self.c1, self.c2
        # Two connected chords
        if b1 < c1:
            for x in range(0, b1):
                self._addPair([(x, b1)], [(x, c1)])
            for x in range(c2+1, self.n):
                if x != b1:
                    self._addPair([(c2, x)], [(b1, x)])
        if b1 > c1:
            for x in range(b1+1, self.n):
                self._addPair([(b1, x)], [(c1, x)])
            for x in range(0, c2):
                if x != b1:
                    self._addPair([(x, c2)], [(x, b1)])
        # Three connected chords
        if b1 < c1:
            for x in range(0, b1):
                for y in range(c1+1, self.n):
                    if y != c2:
                        self._addPair([(x, b1), (c1, y)], [(x, y)])
            for x in range(0, c2):
                for y in range(c2+1, self.n):
                    if x != b1 and x != c1:
                        self._addPair([(x, y)], [(x, c2), (b1, y)])
        if b1 > c1:
            for x in range(0, c1):
                for y in range(b1+1, self.n):
                    if x != c2:
                        self._addPair([(x, c1), (b1, y)], [(x, y)])
            for x in range(0, c2):
                for y in range(c2+1, self.n):
                    if y != b1 and y != c1:
                        self._addPair([(x, y)], [(x, b1), (c2, y)])

    def _O5Chords(self):
        b1, b2, c1, c2 = self.b1, self.b2, self.c1, self.c2
        sc, bc = min(c1, c2), max(c1, c2)
        # Chords can be disjoint
        for x in range(sc+1, bc):
            for y in range(x+1, bc):
                if x not in (b1, b2) and y not in (b1, b2):
                    self._addPair([(sc, x), (y, bc)], [(sc, x), (y, bc)])
        # Or one can be contained in the other
        for x in range(bc+1, self.n):
            for y in range(sc+1, bc):
                if x not in (b1, b2) and y not in (b1, b2):
                    self._addPair([(sc, x), (y, bc)], [(sc, x), (y, bc)])
        for x in range(sc+1, bc):
            for y in range(0, sc):
                if x not in (b1, b2) and y not in (b1, b2):
                    self._addPair([(sc, x), (y, bc)], [(sc, x), (y, bc)])

    def _O_BasicChoice(self):
        # Uses one standard way of forming basic choice
        b1, c1, c2 = self.b1, self.c1, self.c2
        # Type O3: those using sigma', not those using sigma
        if b1 < c1:
            self._addPair([(c1, c2)], [(c1, b1)])
        if b1 > c1:
            self._addPair([(c2, c1)], [(b1, c1)])
        # Type O4: three connected chords, two on left
        if b1 < c1:
            for x in range(0, b1):
                self._addPair([(x, b1), (c1, c2)], [(x, c2)])
        if b1 > c1:
            for y in range(b1+1, self.n):
                self._addPair([(c2, c1), (b1, y)], [(c2, y)])
        # Type O7: those with a break on the left
        if c1 < c2:
            for x in range(c1+1, c2):
                self._addPair([(c1, x), (x, c2)], [(c1, c2)])
        if c2 < c1:
            for x in range(c2+1, c1):
                self._addPair([(c2, x), (x, c1)], [(c2, c1)])
        # Type O8: None chosen

    _OChords = [_O1Chords, _O2Chords, _O3Chords, _O4Chords, _O5Chords, \
                    _O6Chords, _O_BasicChoice]
