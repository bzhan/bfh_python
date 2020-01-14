"""Type DD structures for cobordisms between linear pointed matched circles."""

from ddstructure import DDStrFromChords
from pmc import linearPMC
from pmc import Idempotent, Strands
from utility import memorize

# Two sides for the larger PMC:
LEFT, RIGHT = 0, 1

class Cobordism(object):
    """Represents a cobordism."""
    def __init__(self, genus, c_pair, side):
        """Specifies the genus of the larger (linear PMC), the c-pair at which
        the cobordism occurred, and the side of the larger PMC (LEFT or RIGHT).

        """
        self.genus = genus  # genus of larger PMC
        self.n = 4 * self.genus  # number of points in the larger PMC
        self.c_pair = c_pair
        self.side = side  # LEFT or RIGHT

        if c_pair == 0 or c_pair == 2*self.genus-1:
            self.is_degenerate = True
        else:
            self.is_degenerate = False

        self.large_pmc = linearPMC(self.genus)
        self.small_pmc = linearPMC(self.genus-1)

        # Some special points and pairs
        self.c1, self.c2 = self.large_pmc.pairs[self.c_pair]
        if self.is_degenerate:
            assert self.c2 == self.c1 + 2
            self.p = self.c1 + 1
            self.p_pair = self.large_pmc.pairid[self.p]
        else:
            assert self.c2 == self.c1 + 3
            self.d, self.u = self.c1 + 1, self.c1 + 2
            self.d_pair = self.large_pmc.pairid[self.d]
            self.u_pair = self.large_pmc.pairid[self.u]

        # Construct the to_s dictionary. Keys are points on the large PMC that
        # match points on the small PMC. Value is the point that it matches.
        self.to_s = dict()
        cur_pt = 0
        for i in range(self.n):
            if self.is_degenerate:
                pair_i = self.large_pmc.pairid[i]
                if pair_i == self.c_pair or pair_i == self.p_pair:
                    continue
            else:
                if self.c1 <= i <= self.c2:
                    continue
            self.to_s[i] = cur_pt
            cur_pt += 1

        # The pair_to_s dictionary is similar, but for pairs. The c-pair does
        # not match to anything. In the non-degenerate case, the u and d pairs
        # both match the (u',d') pair on the left. In the degenerate case, the
        # p-pair also does not match anything.
        self.pair_to_s = dict()
        for i in range(self.n//2):
            for p in self.large_pmc.pairs[i]:
                if p in self.to_s:
                    self.pair_to_s[i] = self.small_pmc.pairid[self.to_s[p]]

        if not self.is_degenerate:
            # Special pair on the small PMC
            self.du_pair = self.pair_to_s[self.d_pair]
            assert self.du_pair == self.pair_to_s[self.u_pair]

        if self.side == LEFT:
            self.start_pmc, self.end_pmc = self.large_pmc, self.small_pmc
        else:
            self.start_pmc, self.end_pmc = self.small_pmc, self.large_pmc

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
        chord1 = [(self.to_s[p], self.to_s[q]) for p, q in chord1]
        chord_small = Strands(self.small_pmc, chord1)
        chord2 = [(p, q) for p, q in chord2]
        chord_large = Strands(self.large_pmc, chord2)
        if self.side == LEFT:
            return (chord_large, chord_small.opp())
        else:
            return (chord_small, chord_large.opp())

    @memorize
    def _getIdems(self):
        """Returns the set of possible idempotent-pairs for generators.

        In the non-degenerate case: the c-pair must be on the right, and at most
        one of u and d-pairs are on the right. The left idempotent is the
        complement of the right idempotent, under the mapping given by
        self.pair_to_s.

        In the degenerate case: the c-pair must be on the right, and the p-pair
        must not be on the right. The left idempotent is again the complement
        of the right idempotent.

        """
        all_idems = []

        large_idems = self.large_pmc.getIdempotents()
        for large_idem in large_idems:
            if self.c_pair not in large_idem:
                continue
            if self.is_degenerate:
                if self.p_pair in large_idem:
                    continue
            else:
                if self.u_pair in large_idem and self.d_pair in large_idem:
                    continue
            small_idem_comp = Idempotent(
                self.small_pmc,
                [self.pair_to_s[i] for i in large_idem if i != self.c_pair])
            small_idem = small_idem_comp.comp()
            if self.side == LEFT:
                all_idems.append((large_idem, small_idem.opp()))
            else:
                all_idems.append((small_idem, large_idem.opp()))

        return all_idems

    def _getChords(self):
        """Returns the chords in the RIGHT case. The chords in the LEFT case are
        formed by switching the components of the pair.

        """
        all_chords = []
        for x in range(self.n):
            for y in range(x+1, self.n):
                if x in self.to_s and y in self.to_s:
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
