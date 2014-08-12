"""Dehn twists starting at linear PMC and their type DD structures."""

from algebra import TensorDGAlgebra, TensorGenerator
from algebra import E0
from ddstructure import MorDDtoDDComplex, MorDDtoDDGenerator, \
    SimpleDDGenerator, SimpleDDStructure
from ddstructure import DDStrFromChords, identityDD
from pmc import Idempotent, Strands, StrandDiagram
from pmc import linearPMC
from utility import memorize
from utility import F2, NEG, POS
import itertools

class DehnTwist:
    """Represents a Dehn twist starting at linear PMC."""
    def __init__(self, genus, c_pair, orientation):
        """Specifies genus of the starting pmc, id of the pair of Dehn twist,
        and orientation of the twist (POS or NEG).

        """
        self.genus = genus
        self.orientation = orientation
        self.start_pmc = linearPMC(genus)
        self.end_pmc = self.start_pmc
        self.n = 4 * genus
        self.c1, self.c2 = self.start_pmc.pairs[c_pair]
        self.c_pair = c_pair
        assert self.c2 == self.c1 + 3
        # Two positions between c1 and c2, for (d)own and (u)p
        self.d = self.c1 + 1
        self.u = self.c1 + 2
        self.d_pair = self.start_pmc.pairid[self.d]
        self.u_pair = self.start_pmc.pairid[self.u]

    @memorize
    def getDDStructure(self):
        """Returns the type DD structure corresponding to this dehn twist."""
        self.all_idems = self._getIdems()
        self.all_chords = []
        for chord_type in self._PChords:
            chord_type(self)
        if self.orientation == NEG:
            self.all_chords = [(r.opp(), l.opp()) for l, r in self.all_chords]
        # mult_one = False case is more complicated
        alg1 = self.start_pmc.getAlgebra(mult_one = True)
        alg2 = alg1
        ddstr = DDStrFromChords(alg1, alg2, self.all_idems, self.all_chords)
        assert ddstr.testDelta()
        return ddstr

    def _getIdems(self):
        """Returns the set of possible idempotent-pairs for generators."""
        all_idems = []

        left_idems = self.start_pmc.getIdempotents()
        # Generators of type X (complementary)
        for idem in left_idems:
            all_idems.append((idem, idem.opp().comp()))
        # Generators of type Y (near-complementary)
        for idem in left_idems:
            for rem_pair in (self.d_pair, self.u_pair):
                if self.c_pair in idem and not rem_pair in idem:
                    idem_data = list(idem.comp())
                    idem_data.remove(rem_pair)
                    idem_data.append(self.c_pair)
                    right_idem = Idempotent(self.end_pmc, idem_data).opp()
                    all_idems.append((idem, right_idem))
        return all_idems

    def _addPair(self, data1, data2):
        """Add this pair of chord data."""
        chord_left = Strands(self.start_pmc, data1)
        data2 = [(self.n-1-q, self.n-1-p) for p, q in data2]
        chord_right = Strands(self.end_pmc.opp(), data2)
        self.all_chords.append((chord_left, chord_right))

    def _P1Chords(self):
        for x in range(self.n):
            for y in range(x+1, self.n):
                if any([(x, y) == (self.d, self.u),
                        x < self.c1 and self.c2 <= y,
                        x <= self.c1 and self.c2 < y,
                        y <= self.c1,
                        x >= self.c2]):
                    self._addPair([(x, y)], [(x, y)])

    def _P2Chords(self):
        for x in range(self.c2, self.n):
            self._addPair([(self.d, x)], [(self.c1, x)])
        for x in range(0, self.c1):
            for y in range(self.c2+1, self.n):
                self._addPair([(x, self.c1), (self.d, y)], [(x, y)])
        for x in range(0, self.c1+1):
            self._addPair([(x, self.c2)], [(x, self.u)])
        for x in range(0, self.c1):
            for y in range(self.c2+1, self.n):
                self._addPair([(x, y)], [(x, self.u), (self.c2, y)])

    def _P3Chords(self):
        for x in range(0, self.c1+1):
            self._addPair([(x, self.c2)], [(x, self.d)])
        for x in range(0, self.c1):
            for y in range(self.c2+1, self.n):
                self._addPair([(x, y)], [(x, self.d), (self.c2, y)])
        for x in range(self.c2, self.n):
            self._addPair([(self.u, x)], [(self.c1, x)])
        for x in range(0, self.c1):
            for y in range(self.c2+1, self.n):
                self._addPair([(x, self.c1), (self.u, y)], [(x, y)])

    def _P4Chords(self):
        self._addPair([], [(self.u, self.c2)])
        for x in range(self.c2+1, self.n):
            self._addPair([(self.c2, x)], [(self.u, x)])
        self._addPair([(self.c1, self.d)], [])
        for x in range(0, self.c1):
            self._addPair([(x, self.d)], [(x, self.c1)])

    def _P5Chords(self):
        self._addPair([(self.c1, self.u)], [])
        for x in range(0, self.c1):
            self._addPair([(x, self.u)], [(x, self.c1)])
        self._addPair([], [(self.d, self.c2)])
        for x in range(self.c2+1, self.n):
            self._addPair([(self.c2, x)], [(self.d, x)])

    def _P6Chords(self):
        self._addPair([(self.d, self.u)], [])
        self._addPair([], [(self.d, self.u)])

    def _P7Chords(self):
        for x in range(0, self.c1):
            for y in range(self.c2+1, self.n):
                self._addPair([(x, self.c1), (self.c2, y)],
                              [(x, self.c1), (self.c2, y)])

    def _P8Chords(self):
        for x in range(self.c2, self.n):
            self._addPair([(self.u, x)], [(self.c1, self.d), (self.u, x)])
        for x in range(0, self.c1):
            for y in range(self.c2, self.n):
                self._addPair([(x, self.c1), (self.u, y)],
                              [(x, self.d), (self.u, y)])
        for x in range(0, self.c1+1):
            self._addPair([(x, self.d), (self.u, self.c2)], [(x, self.d)])
        for x in range(0, self.c1):
            for y in range(self.c2+1, self.n):
                self._addPair([(x, self.d), (self.u, y)],
                              [(x, self.d), (self.c2, y)])

    def _P9Chords(self):
        for x in range(0, self.c1):
            for y in range(self.c2+1, self.n):
                self._addPair([(x, self.d), (self.u, y)], [(x, y)])
                self._addPair([(x, y)], [(x, self.d), (self.u, y)])

    _PChords = [_P1Chords, _P2Chords, _P3Chords, _P4Chords, _P5Chords,
                _P6Chords, _P7Chords, _P8Chords, _P9Chords]

class AntiBraid:
    """Represents the anti-braid resolution."""
    def __init__(self, genus, c_pair):
        """Specifies genus of the starting pmc and the id of the pair of
        anti-braid resolution.

        """
        self.genus = genus
        self.c_pair = c_pair
        self.start_pmc = linearPMC(genus)
        self.end_pmc = self.start_pmc

        self.n = 4 * genus
        self.c1, self.c2 = self.start_pmc.pairs[c_pair]

        if self.c2 == self.c1 + 3:
            self.is_degenerate = False
        else:
            assert self.c2 == self.c1 + 2
            assert self.c1 == 0 or self.c2 == self.n - 1
            self.is_degenerate = True

        if self.is_degenerate:
            # One position between c1 and c2, called p
            self.p = self.c1 + 1
            self.p_pair = self.start_pmc.pairid[self.p]
        else:
            # Two positions between c1 and c2, for (d)own and (u)p
            self.d = self.c1 + 1
            self.u = self.c1 + 2
            self.d_pair = self.start_pmc.pairid[self.d]
            self.u_pair = self.start_pmc.pairid[self.u]

    @memorize
    def getDDStructure(self):
        """Returns the type DD structure corresponding to this dehn twist."""
        all_idems = self._getIdems()
        all_chords = []
        for chord_type in self._getChordsList():
            all_chords.extend(chord_type())
        all_chords = [self._StrandsFromChords(chord1, chord2)
                      for chord1, chord2 in all_chords]

        alg1 = self.start_pmc.getAlgebra(mult_one = True)
        alg2 = alg1
        ddstr = DDStrFromChords(alg1, alg2, all_idems, all_chords)
        assert ddstr.testDelta()
        return ddstr

    @memorize
    def getAdmissibleDDStructure(self):
        """Returns the type DD structure corresponding to the Heegaard diagram
        created by a finger move of the beta circle to the right.

        """
        alg1 = self.start_pmc.getAlgebra(mult_one = True)
        alg2 = alg1
        ddstr = SimpleDDStructure(F2, alg1, alg2)

        # Add generators for the non-admissible case - that is, those generators
        # that do not contain the two intersections created by the finger move.
        original_idems = self._getIdems()
        for i in range(len(original_idems)):
            left_idem, right_idem = original_idems[i]
            ddstr.addGenerator(
                SimpleDDGenerator(ddstr, left_idem, right_idem, "0_%d" % i))

        # Now add the new generators. These just correspond to the complementary
        # idempotents with c_pair on the left, repeated twice.
        left_idems = [idem for idem in self.start_pmc.getIdempotents()
                      if self.c_pair in idem]
        for i in range(len(left_idems)):
            left_idem = left_idems[i]
            right_idem = left_idem.opp().comp()
            ddstr.addGenerator(
                SimpleDDGenerator(ddstr, left_idem, right_idem, "1_%d" % i))
            ddstr.addGenerator(
                SimpleDDGenerator(ddstr, left_idem, right_idem, "2_%d" % i))

        gen_set = []
        for i in range(3):
            gen_set.append([gen for gen in ddstr.getGenerators()
                            if gen.name[:1] == "%d" % i])

        # Enumerate the non-special chords (those that do not dependent on the
        # idempotent. See the functions themselves for the format of all_chords.
        if self.is_degenerate:
            all_chords = self._getAdmissibleNonSpecialChordsDegenerate()
        else:
            all_chords = self._getAdmissibleNonSpecialChords()
        for i, j in itertools.product(range(3), range(3)):
            all_chords[i][j] = [self._StrandsFromChords(chord1, chord2)
                                for chord1, chord2 in all_chords[i][j]]

        # Now we emulate the logic in ddstructure.DDStrFromChords, except we
        # distinguish between ''classes'' of generators, by the first character
        # of the name of the generator.
        for i, j in itertools.product(range(3), range(3)):
            for x, y in itertools.product(gen_set[i], gen_set[j]):
                for l_chord, r_chord in all_chords[i][j]:
                    if l_chord.idemCompatible(x.idem1, y.idem1) and \
                            r_chord.idemCompatible(x.idem2, y.idem2):
                        ddstr.addDelta(x, y,
                                       StrandDiagram(alg1, x.idem1, l_chord),
                                       StrandDiagram(alg2, x.idem2, r_chord), 1)

        # Special handling for these. From class 2 to class 1, add only if the
        # c-pair is occupied on the left side (and not on the right).
        # Non-degenerate cases only.
        sp_chords = []
        if not self.is_degenerate:
            for x in range(0, self.c1):
                for y in range(self.c2+1, self.n):
                    sp_chords.append(([(x, y)], [(x, self.u), (self.u, y)]))
                    sp_chords.append(([(x, y)], [(x, self.d), (self.d, y)]))
                    sp_chords.append(([(x, self.d), (self.u, y)],
                                      [(x, self.d), (self.u, y)]))

        sp_chords = [self._StrandsFromChords(chord1, chord2)
                     for chord1, chord2 in sp_chords]
        for x, y in itertools.product(gen_set[2], gen_set[1]):
            for l_chord, r_chord in sp_chords:
                if self.c_pair in x.idem1 and \
                        l_chord.idemCompatible(x.idem1, y.idem1) and \
                        r_chord.idemCompatible(x.idem2, y.idem2):
                    assert self.c_pair not in x.idem2.opp() and \
                        self.c_pair in y.idem1 and \
                        self.c_pair not in y.idem2.opp()
                    ddstr.addDelta(x, y,
                                   StrandDiagram(alg1, x.idem1, l_chord),
                                   StrandDiagram(alg2, x.idem2, r_chord), 1)

        assert ddstr.testDelta()
        return ddstr

    def _getAdmissibleNonSpecialChordsDegenerate(self):
        """Returns the non-special chords (those that do not depend on the
        idempotent in the degenerate admissible case.

        """
        # Initialize all_chords to be a 3*3 matrix of lists, with each entry
        # (i, j) containing the chord pairs from generators of class i to
        # generators of class j.
        all_chords = []
        for i in range(3):
            all_chords.append([])
            for j in range(3):
                all_chords[i].append([])

        # Idempotent actions from 1 to 2.
        all_chords[1][2].append(([], []))

        # Basic chords.
        all_chords[0][2].append(([], [(self.p, self.c2)]))
        all_chords[1][0].append(([], [(self.c1, self.p)]))

        # Incorporate the left side.
        all_chords[0][0].append(([(self.c1, self.c2)], []))
        all_chords[2][0].append(([(self.c1, self.c2)], [(self.c1, self.p)]))
        all_chords[0][1].append(([(self.c1, self.c2)], [(self.p, self.c2)]))

        # Identity away from the anti-braid.
        for x in range(0, self.n):
            for y in range(x+1, self.n):
                if y < self.c1 or x > self.c2:
                    for i in range(3):
                        all_chords[i][i].append(([(x, y)], [(x, y)]))
        return all_chords

    def _getAdmissibleNonSpecialChords(self):
        """Returns the non-special chords (those that do not depend on the
        idempotent in the non-degenerate admissible case.

        """
        # Initialize all_chords to be a 3*3 matrix of lists, with each entry
        # (i, j) containing the chord pairs from generators of class i to
        # generators of class j.
        all_chords = []
        for i in range(3):
            all_chords.append([])
            for j in range(3):
                all_chords[i].append([])

        # Idempotent actions from 1 to 2.
        all_chords[1][2].append(([], []))

        # Basic chords.
        all_chords[0][0].append(([], [(self.d, self.u)]))
        all_chords[0][2].append(([], [(self.u, self.c2)]))
        all_chords[1][0].append(([], [(self.c1, self.d)]))
        all_chords[0][2].append(([], [(self.d, self.c2)]))
        all_chords[1][0].append(([], [(self.c1, self.u)]))

        # Incorporate the intervals (c1-1, c1) and (c2, c2+1).
        for x in range(0, self.c1):
            for y in range(self.c2+1, self.n):
                all_chords[0][0].append(
                    ([(x, self.c1), (self.c2, y)],
                     [(x, self.c1), (self.c2, y)]))
                all_chords[2][0].append(
                    ([(x, self.c1), (self.c2, y)], [(x, self.u), (self.c2, y)]))
                all_chords[2][0].append(
                    ([(x, self.c1), (self.c2, y)], [(x, self.d), (self.c2, y)]))
                all_chords[0][1].append(
                    ([(x, self.c1), (self.c2, y)], [(x, self.c1), (self.u, y)]))
                all_chords[0][1].append(
                    ([(x, self.c1), (self.c2, y)], [(x, self.c1), (self.d, y)]))
                for i in range(3):
                    all_chords[i][i].append(
                        ([(x, self.c1), (self.c2, y)], [(x, y)]))

        # Incorporate the left side.
        all_chords[0][0].extend(
            [([(self.d, self.u)], []),
             ([(self.c1, self.d), (self.u, self.c2)], []),
             ([(self.c1, self.c2)], [])])
        for i in (1, 2):
            all_chords[i][i].append(([(self.d, self.u)], [(self.d, self.u)]))

        all_chords[2][0].extend(
            [([(self.c1, self.d), (self.u, self.c2)], [(self.c1, self.d)]),
             ([(self.c1, self.c2)], [(self.c1, self.d)]),
             ([(self.c1, self.d), (self.u, self.c2)], [(self.c1, self.u)]),
             ([(self.c1, self.c2)], [(self.c1, self.u)])])

        all_chords[0][1].extend(
            [([(self.c1, self.d), (self.u, self.c2)], [(self.u, self.c2)]),
             ([(self.c1, self.c2)], [(self.u, self.c2)]),
             ([(self.c1, self.d), (self.u, self.c2)], [(self.d, self.c2)]),
             ([(self.c1, self.c2)], [(self.d, self.c2)])])

        # Again add intervals (c1-1, c1) and (c2, c2+1).
        for x in range(0, self.c1):
            for y in range(self.c2+1, self.n):
                all_chords[0][0].extend(
                    [([(x, y)], [(x, y)]),
                     ([(x, y)], [(x, self.d), (self.u, y)]),
                     ([(x, y)], [(x, self.c1), (self.c2, y)]),
                     ([(x, self.c1), (self.c2, y)], [(x, self.d), (self.u, y)]),
                     ([(x, self.d), (self.u, y)], [(x, self.c1), (self.c2, y)]),
                     ([(x, self.d), (self.u, y)], [(x, y)])])

                all_chords[2][0].extend(
                    [([(x, y)], [(x, self.d), (self.c2, y)]),
                     ([(x, y)], [(x, self.u), (self.c2, y)]),
                     ([(x, self.d), (self.u, y)], [(x, self.d), (self.c2, y)]),
                     ([(x, self.d), (self.u, y)], [(x, self.u), (self.c2, y)])])

                all_chords[0][1].extend(
                    [([(x, y)], [(x, self.c1), (self.u, y)]),
                     ([(x, y)], [(x, self.c1), (self.d, y)]),
                     ([(x, self.d), (self.u, y)], [(x, self.c1), (self.u, y)]),
                     ([(x, self.d), (self.u, y)], [(x, self.c1), (self.d, y)])])

        # Identity away from the anti-braid.
        for x in range(0, self.n):
            for y in range(x+1, self.n):
                if y < self.c1 or x > self.c2:
                    for i in range(3):
                        all_chords[i][i].append(([(x, y)], [(x, y)]))
        return all_chords

    def _getIdems(self):
        """Returns the set of possible idempotent-pairs for generators."""
        all_idems = []

        left_idems = self.start_pmc.getIdempotents()
        # Near-complementary generators
        if self.is_degenerate:
            rem_pairs = (self.p_pair,)
        else:
            rem_pairs = (self.d_pair, self.u_pair)

        for idem in left_idems:
            for rem_pair in rem_pairs:
                if self.c_pair in idem and not rem_pair in idem:
                    idem_data = list(idem.comp())
                    idem_data.remove(rem_pair)
                    idem_data.append(self.c_pair)
                    right_idem = Idempotent(self.end_pmc, idem_data).opp()
                    all_idems.append((idem, right_idem))
        return all_idems

    def _StrandsFromChords(self, chord1, chord2):
        """Create strand objects from lists of chords. Points in chord2 are
        reversed (refer to the opposite pmc).

        """
        chord_left = Strands(self.start_pmc, chord1)
        chord2 = [(self.n-1-q, self.n-1-p) for p, q in chord2]
        chord_right = Strands(self.end_pmc.opp(), chord2)
        return (chord_left, chord_right)

    def _B1Chords(self):
        return [([(x, y)], [(x, y)])
                for x in range(self.n) for y in range(x+1, self.n)
                if (x < self.c1 or x > self.c2) and \
                    (y < self.c1 or y > self.c2)]

    def _B2Chords(self):
        return [([(self.d, self.u)], []), ([], [(self.d, self.u)])]

    def _B3Chords(self):
        return [([(self.c1, self.d), (self.u, self.c2)], []),
                ([], [(self.c1, self.d), (self.u, self.c2)])]

    def _B4Chords(self):
        return [([(self.c1, self.c2)], []), ([], [(self.c1, self.c2)])]

    def _B5Chords(self):
        result = []
        for x in range(0, self.c1):
            for y in range(self.c2+1, self.n):
                result.append(([(x, self.c1), (self.c2, y)], [(x, y)]))
                result.append(([(x, y)], [(x, self.c1), (self.c2, y)]))
        return result

    def _B6Chords(self):
        result = []
        for x in range(0, self.c1):
            for y in range(self.c2+1, self.n):
                result.append(([(x, self.d), (self.u, y)], [(x, y)]))
                result.append(([(x, y)], [(x, self.d), (self.u, y)]))
        return result

    def _B7Chords(self):
        return [([(x, self.c1), (self.c2, y)], [(x, self.c1), (self.c2, y)])
                for x in range(0, self.c1) for y in range(self.c2+1, self.n)]

    def _B8Chords(self):
        result = []
        for x in range(0, self.c1):
            for y in range(self.c2+1, self.n):
                result.append(([(x, self.d), (self.u, y)],
                               [(x, self.c1), (self.c2, y)]))
                result.append(([(x, self.c1), (self.c2, y)],
                               [(x, self.d), (self.u, y)]))
        return result

    def _getChordsList(self):
        if self.is_degenerate:
            return [self._B1Chords, self._B4Chords]
        else:
            return [self._B1Chords, self._B2Chords, self._B3Chords,
                    self._B4Chords, self._B5Chords, self._B6Chords,
                    self._B7Chords, self._B8Chords]

class DehnSurgery:
    """Represents a Dehn twist starting at linear PMC."""
    def __init__(self, genus, c_pair, orientation):
        """Specifies genus of the starting pmc, id of the pair of Dehn twist,
        and orientation of the twist (POS or NEG).

        """
        self.genus = genus
        self.orientation = orientation
        self.start_pmc = linearPMC(genus)
        self.end_pmc = self.start_pmc
        self.n = 4 * genus
        self.c1, self.c2 = self.start_pmc.pairs[c_pair]
        self.c_pair = c_pair

        if self.c2 == self.c1 + 3:
            self.is_degenerate = False
        else:
            assert self.c2 == self.c1 + 2
            assert self.c1 == 0 or self.c2 == self.n - 1
            self.is_degenerate = True

        if not self.is_degenerate:
            # Two positions between c1 and c2, for (d)own and (u)p
            self.d = self.c1 + 1
            self.u = self.c1 + 2

    def _addPair(self, data1, data2):
        """Add this pair of chord data."""
        chord_left = Strands(self.start_pmc, data1)
        data2 = [(self.n-1-q, self.n-1-p) for p, q in data2]
        chord_right = Strands(self.end_pmc.opp(), data2)
        self.all_chords.append((chord_left, chord_right))

    @memorize
    def getMorphism(self):
        id_dd = identityDD(self.start_pmc)
        anti_braid = AntiBraid(self.genus, self.c_pair)
        ab_dd = anti_braid.getDDStructure()

        if self.orientation == NEG:
            source = id_dd
            target = ab_dd
        else:
            source = ab_dd
            target = id_dd

        morphism_cx = MorDDtoDDComplex(F2, source, target)

        self.all_chords = []
        for chord_type in self._getChordsList():
            chord_type()

        morphism = E0
        alg = self.start_pmc.getAlgebra()
        assert alg.mult_one is True  # not prepared for the other case
        tensor_alg = TensorDGAlgebra((alg, alg))
        # Similar to method in DDStrFromChords
        for x in source.getGenerators():
            for y in target.getGenerators():
                for l_chord, r_chord in self.all_chords:
                    if l_chord.idemCompatible(x.idem1, y.idem1) and \
                            r_chord.idemCompatible(x.idem2, y.idem2):
                        a1 = StrandDiagram(alg, x.idem1, l_chord)
                        a2 = StrandDiagram(alg, x.idem2, r_chord)
                        coeff = TensorGenerator((a1, a2), tensor_alg)
                        morphism += 1*MorDDtoDDGenerator(
                            morphism_cx, x, coeff, y)
        return morphism

    def _N1Chords(self):
        self._addPair([(self.c2-1, self.c2)], [])
        self._addPair([], [(self.c1, self.c1+1)])

    def _N2Chords(self):
        for x in range(self.c2+1, self.n):
            self._addPair([(self.c2-1, x)], [(self.c2, x)])
        for x in range(0, self.c1):
            self._addPair([(x, self.c1)], [(x, self.c1+1)])

    def _N3Chords(self):
        self._addPair([(self.d, self.c2)], [])
        self._addPair([], [(self.c1, self.u)])

    def _N4Chords(self):
        for x in range(self.c2+1, self.n):
            self._addPair([(self.d, x)], [(self.c2, x)])
        for x in range(0, self.c1):
            self._addPair([(x, self.c1)], [(x, self.u)])

    def _P1Chords(self):
        self._addPair([(self.c1, self.c1+1)], [])
        self._addPair([], [(self.c2-1, self.c2)])

    def _P2Chords(self):
        for x in range(0, self.c1):
            self._addPair([(x, self.c1+1)], [(x, self.c1)])
        for x in range(self.c2+1, self.n):
            self._addPair([(self.c2, x)], [(self.c2-1, x)])

    def _P3Chords(self):
        self._addPair([(self.c1, self.u)], [])
        self._addPair([], [(self.d, self.c2)])

    def _P4Chords(self):
        for x in range(0, self.c1):
            self._addPair([(x, self.u)], [(x, self.c1)])
        for x in range(self.c2+1, self.n):
            self._addPair([(self.c2, x)], [(self.d, x)])

    def _getChordsList(self):
        if self.orientation == NEG:
            if self.is_degenerate:
                return [self._N1Chords, self._N2Chords]
            else:
                return [self._N1Chords, self._N2Chords, self._N3Chords,
                        self._N4Chords]
        else:
            if self.is_degenerate:
                return [self._P1Chords, self._P2Chords]
            else:
                return [self._P1Chords, self._P2Chords, self._P3Chords,
                        self._P4Chords]
