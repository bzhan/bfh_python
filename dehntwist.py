"""Dehn twists starting at linear PMC and their type DD structures."""

from algebra import TensorDGAlgebra, TensorGenerator
from algebra import E0
from ddstructure import MorDDtoDDComplex, MorDDtoDDGenerator
from ddstructure import DDStrFromChords, identityDD
from pmc import Idempotent, Strands, StrandDiagram
from pmc import linearPMC
from utility import memorize
from utility import F2, NEG, POS

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
        self.all_idems = self._getIdems()
        self.all_chords = []
        for chord_type in self._getChordsList():
            chord_type()

        alg1 = self.start_pmc.getAlgebra(mult_one = True)
        alg2 = alg1
        ddstr = DDStrFromChords(alg1, alg2, self.all_idems, self.all_chords)
        assert ddstr.testDelta()
        return ddstr

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

    def _addPair(self, data1, data2):
        """Add this pair of chord data."""
        chord_left = Strands(self.start_pmc, data1)
        data2 = [(self.n-1-q, self.n-1-p) for p, q in data2]
        chord_right = Strands(self.end_pmc.opp(), data2)
        self.all_chords.append((chord_left, chord_right))

    def _B1Chords(self):
        for x in range(self.n):
            for y in range(x+1, self.n):
                if (x < self.c1 or x > self.c2) and \
                        (y < self.c1 or y > self.c2):
                    self._addPair([(x, y)], [(x, y)])

    def _B2Chords(self):
        self._addPair([(self.d, self.u)], [])
        self._addPair([], [(self.d, self.u)])

    def _B3Chords(self):
        self._addPair([(self.c1, self.d), (self.u, self.c2)], [])
        self._addPair([], [(self.c1, self.d), (self.u, self.c2)])

    def _B4Chords(self):
        self._addPair([(self.c1, self.c2)], [])
        self._addPair([], [(self.c1, self.c2)])

    def _B5Chords(self):
        for x in range(0, self.c1):
            for y in range(self.c2+1, self.n):
                self._addPair([(x, self.c1), (self.c2, y)], [(x, y)])
                self._addPair([(x, y)], [(x, self.c1), (self.c2, y)])

    def _B6Chords(self):
        for x in range(0, self.c1):
            for y in range(self.c2+1, self.n):
                self._addPair([(x, self.d), (self.u, y)], [(x, y)])
                self._addPair([(x, y)], [(x, self.d), (self.u, y)])

    def _B7Chords(self):
        for x in range(0, self.c1):
            for y in range(self.c2+1, self.n):
                self._addPair([(x, self.c1), (self.c2, y)],
                              [(x, self.c1), (self.c2, y)])

    def _B8Chords(self):
        for x in range(0, self.c1):
            for y in range(self.c2+1, self.n):
                self._addPair([(x, self.d), (self.u, y)],
                              [(x, self.c1), (self.c2, y)])
                self._addPair([(x, self.c1), (self.c2, y)],
                              [(x, self.d), (self.u, y)])

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
