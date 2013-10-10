"""Dehn twists starting at linear PMC and their type DD structures."""

from ddstructure import DDStrFromChords
from pmc import Idempotent, Strands
from pmc import linearPMC
from utility import memorize
from utility import NEG, POS

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
