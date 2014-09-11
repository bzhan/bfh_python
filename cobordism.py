"""Type DD structures for cobordisms between linear pointed matched circles."""

from ddstructure import DDStrFromChords
from pmc import linearPMC
from pmc import Idempotent, Strands
from utility import memorize

# Two sides for the larger PMC:
LEFT, RIGHT = 0, 1

class Cobordism:
    """Represents a cobordism."""
    def __init__(self, genus, c_pair, side):
        """Specifies the genus of the starting (linear PMC), the c-pair at which
        the cobordism occurred, and the side of the larger PMC (LEFT or RIGHT).

        """
        self.genus = genus  # genus of larger PMC
        self.n = 4 * self.genus  # number of points in the larger PMC
        self.c_pair = c_pair
        self.side = side  # LEFT or RIGHT

        if side == LEFT:
            # Just keep the dual cobordism and use it for everything.
            self.dual_cobordism = Cobordism(genus, c_pair, RIGHT)
            self.start_pmc = linearPMC(self.genus)
            self.end_pmc = linearPMC(self.genus-1)
            return

        assert side == RIGHT
        if c_pair == 0 or c_pair == 2*self.genus-1:
            self.is_degenerate = True
        else:
            self.is_degenerate = False

        self.start_pmc = linearPMC(self.genus-1)
        self.end_pmc = linearPMC(self.genus)

        # Some special points and pairs
        self.c1, self.c2 = self.end_pmc.pairs[self.c_pair]
        if self.is_degenerate:
            assert self.c2 == self.c1 + 2
            self.p = self.c1 + 1
            self.p_pair = self.end_pmc.pairid[self.p]
        else:
            assert self.c2 == self.c1 + 3
            self.d, self.u = self.c1 + 1, self.c1 + 2
            self.d_pair = self.end_pmc.pairid[self.d]
            self.u_pair = self.end_pmc.pairid[self.u]

        # Construct the to_l dictionary. Keys are points on the right that match
        # points on the left. Value is the point that it matches.
        self.to_l = dict()
        cur_pt = 0
        for i in range(self.n):
            if self.is_degenerate:
                pair_i = self.end_pmc.pairid[i]
                if pair_i == self.c_pair or pair_i == self.p_pair:
                    continue
            else:
                if self.c1 <= i <= self.c2:
                    continue
            self.to_l[i] = cur_pt
            cur_pt += 1

        # The pair_to_l dictionary is similar, but for pairs. The c-pair does
        # not match to anything. In the non-degenerate case, the u and d pairs
        # both match the (u',d') pair on the left. In the degenerate case, the
        # p-pair also does not match anything.
        self.pair_to_l = dict()
        for i in range(self.n/2):
            for p in self.end_pmc.pairs[i]:  # the larger PMC
                if p in self.to_l:
                    self.pair_to_l[i] = self.start_pmc.pairid[self.to_l[p]]

        if not self.is_degenerate:
            self.du_pair = self.pair_to_l[self.d_pair]  # special pair at left
            assert self.du_pair == self.pair_to_l[self.u_pair]

    def __eq__(self, other):
        return self.genus == other.genus and self.c_pair == other.c_pair and \
            self.side == other.side

    def __ne__(self, other):
        return not (self == other)

    def __hash__(self):
        return hash((self.genus, self.c_pair, self.side, "Cobordism"))

    @memorize
    def getDDStructure(self):
        """Returns the type DD structure corresponding to this cobordism."""
        if self.side == LEFT:
            return self.dual_cobordism.getDDStructure().dual()

        assert self.side == RIGHT
        all_idems = self._getIdems()
        all_chords = self._getChords()
        all_chords = [self._StrandsFromChords(chord1, chord2)
                      for chord1, chord2 in all_chords]

        alg1 = self.start_pmc.getAlgebra(mult_one = True)
        alg2 = self.end_pmc.getAlgebra(mult_one = True)
        ddstr = DDStrFromChords(alg1, alg2, all_idems, all_chords)
        assert ddstr.testDelta()
        return ddstr

    def _StrandsFromChords(self, chord1, chord2):
        """Create strand objects from lists of chords. Points in chord2 are
        reversed (refer to the opposite pmc).

        """
        assert self.side == RIGHT
        chord1 = [(self.to_l[p], self.to_l[q]) for p, q in chord1]
        chord_left = Strands(self.start_pmc, chord1)
        chord2 = [(self.n-1-q, self.n-1-p) for p, q in chord2]
        chord_right = Strands(self.end_pmc, chord2)
        return (chord_left, chord_right)

    @memorize
    def _getIdems(self):
        """Returns the set of possible idempotent-pairs for generators.

        In the non-degenerate case: the c-pair must be on the right, and at most
        one of u and d-pairs are on the right. The left idempotent is the
        complement of the right idempotent, under the mapping given by
        self.pair_to_l.

        In the degenerate case: the c-pair must be on the right, and the p-pair
        must not be on the right. The left idempotent is again the complement
        of the right idempotent.

        """
        assert self.side == RIGHT
        all_idems = []

        right_idems = self.end_pmc.getIdempotents()
        for right_idem in right_idems:
            if self.c_pair not in right_idem:
                continue
            if self.is_degenerate:
                if self.p_pair in right_idem:
                    continue
            else:
                if self.u_pair in right_idem and self.d_pair in right_idem:
                    continue
            left_idem_comp = Idempotent(
                self.start_pmc,
                [self.pair_to_l[i] for i in right_idem if i != self.c_pair])
            left_idem = left_idem_comp.comp()
            all_idems.append((left_idem, right_idem.opp()))

        return all_idems

    def _getChords(self):
        """Returns the chords in the RIGHT case."""
        assert self.side == RIGHT
        all_chords = []
        for x in range(self.n):
            for y in range(x+1, self.n):
                if x in self.to_l and y in self.to_l:
                    all_chords.append(([(x, y)], [(x, y)]))

        all_chords.append(([], [(self.c1, self.c2)]))

        if not self.is_degenerate:
            all_chords.append(([], [(self.d, self.u)]))
            all_chords.append(([], [(self.c1, self.d), (self.u, self.c2)]))

            for x in range(0, self.c1):
                for y in range(self.c2+1, self.n):
                    all_chords.append(([(x, y)], [(x, self.c1), (self.c2, y)]))
                    all_chords.append(([(x, y)], [(x, self.d), (self.u, y)]))

        return all_chords
