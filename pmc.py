"""Pointed matched circle and its algebras."""

from algebra import *
from grading import *
from fractions import Fraction
import itertools

class PMC:
    """Represents a pointed matched circle."""
    def __init__(self, matching):
        """Creates a pointed matched circle from a list of matched pairs. The
        indices start at 0.

        """
        self.n = len(matching) * 2
        self.num_pair = self.n / 2
        self.genus = self.n / 4
        # otherp[i] is the point paired to i (0 <= i < n)
        self.otherp = [0] * self.n
        for p, q in matching:
            self.otherp[p] = q
            self.otherp[q] = p
        # pairid[i] is the ID of the pair containing i (0 <= pairid[i] < n/2)
        self.pairid = [-1] * self.n
        # pairs[i] is the pair of points with ID i (0 <= i < n/2)
        self.pairs = []
        pairCount = 0
        for pos in range(self.n):
            if self.pairid[pos] == -1:
                self.pairid[pos] = self.pairid[self.otherp[pos]] = pairCount
                self.pairs.append((pos, self.otherp[pos]))
                pairCount += 1

    def __eq__(self, other):
        return self.otherp == other.otherp

    def __ne__(self, other):
        return not (self == other)

    @memorizeHash
    def __hash__(self):
        return hash(tuple(self.otherp))

    def __str__(self):
        return str(self.pairs)

    def __repr__(self):
        return "PMC(%s)" % str(self)

    @memorize
    def opp(self):
        """Returns a new PMC that represents the opposite PMC."""
        return PMC([(self.n-1-p, self.n-1-q) for p, q in self.pairs])

    def sd(self, data, mult_one = MULT_ONE):
        """Simple way to obtain a strand diagram for this PMC. Each element of
        data is either an integer or a pair. An integer specifies a double
        horizontal at this position (and its paired position). A pair (p, q)
        specifies a strand from p to q.

        """
        idem_size = len(data)
        parent = StrandAlgebra(F2, self, idem_size, mult_one)
        left_idem = []
        strands = []
        for d in data:
            if isinstance(d, int):
                left_idem.append(self.pairid[d])
            else:
                left_idem.append(self.pairid[d[0]])
                strands.append(d)
        return StrandDiagram(parent, left_idem, strands)

    def idem(self, data):
        """Simple way to obtain an idempotent for this PMC. Each element of
        data specifies a double horizontal at this position (and its paired
        position).

        """
        idem_data = [self.pairid[d] for d in data]
        return Idempotent(self, idem_data)

    def big_gr(self, maslov, spinc):
        """Simple way to obtain an element of the big grading group for this
        PMC.

        """
        grading_group = BigGradingGroup(self)
        return BigGradingElement(grading_group, maslov, spinc)

    def small_gr(self, maslov, spinc):
        """Simple way to obtain an element of the small grading group for this
        PMC.

        """
        grading_group = SmallGradingGroup(self)
        return SmallGradingElement(grading_group, maslov, spinc)

    def getAlgebra(self, idem_size = None, mult_one = MULT_ONE):
        """Returns the algebra with a given size of idempotent (the default
        value, with size half the number of pairs, is most used).

        """
        if idem_size == None: idem_size = self.genus
        return StrandAlgebra(F2, self, idem_size, mult_one)

    def getIdempotents(self, idem_size = None):
        """Get the list of all idempotents."""
        if idem_size == None: idem_size = self.genus
        return [Idempotent(self, data) for data in
                itertools.combinations(range(self.num_pair), idem_size)]

    def getStrandDiagrams(self, idem_size = None):
        """Get the list of generators of the (full) strand algebra."""
        if idem_size == None: idem_size = self.genus
        algebra = self.getAlgebra(idem_size, mult_one = False)
        result = []
        def helper(l_idem, r_idem, strands, pos):
            # Both l_idem and r_idem are lists of pair ID's. The first
            # 'pos' of them are already uesd to generate strands or double
            # horizontals. 'strands' keep track of strands generated.
            if pos == idem_size:
                result.append(StrandDiagram(algebra, l_idem, strands))
                return
            for i in range(pos, idem_size):
                r_idem[i], r_idem[pos] = r_idem[pos], r_idem[i]
                if l_idem[pos] == r_idem[pos]:
                    helper(l_idem, r_idem, strands, pos+1)
                for p in self.pairs[l_idem[pos]]:
                    for q in self.pairs[r_idem[pos]]:
                        if p < q:
                            helper(l_idem, r_idem, strands + [(p, q)], pos+1)
                r_idem[i], r_idem[pos] = r_idem[pos], r_idem[i]
        idems = self.getIdempotents(idem_size)
        for l_idem in idems:
            for r_idem in idems:
                helper(list(l_idem), list(r_idem), [], 0)
        return result

    def getMultOneStrandDiagrams(self, idem_size = None):
        """Get the list of generators of the multiplicity one strand algebra
        (exclude those with multiplicity greater than one at some position).

        """
        allSDs = self.getStrandDiagrams(idem_size)
        multOneSDs = []
        algebra = self.getAlgebra(idem_size, mult_one = True)
        for sd in allSDs:
            if all([x <= 1 for x in sd.multiplicity]):
                multOneSDs.append(
                    StrandDiagram(algebra, sd.left_idem, sd.strands))
        return multOneSDs

def splitPMC(genus):
    """Returns the split pmc with a given genus."""
    return PMC(sum([[(4*i, 4*i+2), (4*i+1, 4*i+3)] for i in range(0,genus)],[]))

def linearPMC(genus):
    """Returns the linear pmc with a given genus."""
    matching = [(0, 2),(4*genus-3, 4*genus-1)]
    matching += [(2*i-1, 2*i+2) for i in range(1, 2*genus-1)]
    return PMC(matching)

def antipodalPMC(genus):
    """Returns the antipodal pmc with a given genus."""
    return PMC([(i, 2*genus+i) for i in range(2*genus)])

def connectSumPMC(pmc1, pmc2):
    """Return the connect sum of two PMC's."""
    pairs2 = [(p+pmc1.n, q+pmc1.n) for p, q in pmc2.pairs]
    return PMC(pmc1.pairs + pairs2)

@memorize
def unconnectSumPMC(pmc, genus1):
    """Returns a pair (pmc1, pmc2) such that pmc1 has genus1 and
    pmc1 # pmc2 = pmc.

    """
    cut_point = 4 * genus1
    for p, q in pmc.pairs:
        assert (p < cut_point and q < cut_point) or \
            (p >= cut_point and q >= cut_point)
    pmc1 = PMC([(p, q) for p, q in pmc.pairs if p < cut_point])
    pmc2 = PMC([(p-cut_point, q-cut_point)
                for p, q in pmc.pairs if p >= cut_point])
    return (pmc1, pmc2)

class Idempotent(tuple):
    """Represents an idempotent in a certain PMC. Stored as a tuple of pairid
    of occupied pairs.

    """
    def __new__(cls, pmc, data):
        return tuple.__new__(cls, tuple(sorted(data)))

    def __init__(self, pmc, data):
        self.pmc = pmc

    def __eq__(self, other):
        if isinstance(other, Idempotent):
            return self.pmc == other.pmc and tuple.__eq__(self, other)
        else:
            return False

    def __ne__(self, other):
        return not (self == other)

    def __hash__(self):
        return hash((self.pmc, tuple(self), "Idempotent"))

    def __str__(self):
        return repr(self)

    def __repr__(self):
        return "(%s)" % ",".join(str(self.pmc.pairs[i]) for i in self)

    @memorize
    def opp(self):
        """Get the same idempotent in the opposite PMC."""
        pmc, pmcopp = self.pmc, self.pmc.opp()
        return Idempotent(pmcopp, [pmcopp.pairid[pmc.n-1-pmc.pairs[i][0]]
                                   for i in self])

    def comp(self):
        """Get the complementary idempotent in the same PMC."""
        return Idempotent(self.pmc,
                          set(range(self.pmc.num_pair))-set(self))

    def toAlgElt(self, parent):
        """Get the strand diagram corresponding to this idempotent, in the
        specified strand algebra.

        """
        return StrandDiagram(parent, self, [])

def unconnectSumIdem(idem, genus1):
    """Returns the pair of idempotents (idem1, idem2) in (pmc1, pmc2), where
    (pmc1, pmc2) is the pair returned by unconnectSumPMC(idem.pmc, genus1).

    """
    cut_pair = 2 * genus1
    pmc1, pmc2 = unconnectSumPMC(idem.pmc, genus1)
    return (Idempotent(pmc1, [pair for pair in idem if pair < cut_pair]),
            Idempotent(pmc2,
                       [pair-cut_pair for pair in idem if pair >= cut_pair]))

class Strands(tuple):
    """Represents a (fixed) list of strands in a certain PMC. Stored as a tuple
    of pairs.

    """
    def __new__(cls, pmc, data):
        return tuple.__new__(cls, tuple(sorted(data)))

    def __init__(self, pmc, data):
        self.pmc = pmc
        # Compute multiplicity at each interval
        self.multiplicity = [0] * (self.pmc.n - 1)
        for st in self:
            assert len(st) == 2 and st[0] < st[1]
            for pos in range(st[0], st[1]):
                self.multiplicity[pos] += 1

    def __eq__(self, other):
        if isinstance(other, Strands):
            return self.pmc == other.pmc and tuple.__eq__(self, other)
        else:
            return False

    def __ne__(self, other):
        return not (self == other)

    def __hash__(self):
        return hash((self.pmc, tuple(self), "Strands"))

    def __str__(self):
        return "(%s)" % ",".join("%d->%d" % (s, t) for s, t in self)

    def __repr__(self):
        return str(self)

    def opp(self):
        """Returns the same strands (with direction reversed) in the opposite
        PMC.

        """
        n = self.pmc.n
        return Strands(self.pmc.opp(), [(n-1-q, n-1-p) for p, q in self])

    def leftCompatible(self, idem):
        """Tests whether this set of strands is compatible with a given left
        idempotent.

        """
        return self.propagateRight(idem) is not None

    def rightCompatible(self, idem):
        """Tests whether this set of strands is compatible with a given right
        idempotent.

        """
        return self.propagateLeft(idem) is not None

    def idemCompatible(self, left_idem, right_idem):
        """Tests whether this set of strands is compatible with the given
        idempotents on the two sides. Note this is not the same as left and
        right compatible.

        """
        return self.leftCompatible(left_idem) and \
            self.propagateRight(left_idem) == right_idem

    def propagateRight(self, left_idem):
        """Find the right_idem given left_idem and strand info. If not
        compatible, return None.

        """
        idemCount = [0] * self.pmc.num_pair
        for pair in left_idem:
            idemCount[pair] += 1
        for st in self:
            if idemCount[self.pmc.pairid[st[0]]] == 0: return None
            idemCount[self.pmc.pairid[st[0]]] -= 1
        for st in self:
            if idemCount[self.pmc.pairid[st[1]]] == 1: return None
            idemCount[self.pmc.pairid[st[1]]] += 1
        right_idem = [i for i in range(self.pmc.n/2) if idemCount[i] == 1]
        return Idempotent(self.pmc, right_idem)

    def propagateLeft(self, right_idem):
        """Find the left_idem given right_idem and strand info. If not
        compatible, return None.

        """
        idemCount = [0] * self.pmc.num_pair
        for pair in right_idem:
            idemCount[pair] += 1
        for st in self:
            if idemCount[self.pmc.pairid[st[1]]] == 0: return None
            idemCount[self.pmc.pairid[st[1]]] -= 1
        for st in self:
            if idemCount[self.pmc.pairid[st[0]]] == 1: return None
            idemCount[self.pmc.pairid[st[0]]] += 1
        left_idem = [i for i in range(self.pmc.n/2) if idemCount[i] == 1]
        return Idempotent(self.pmc, left_idem)

    def isMultOne(self):
        """Tests whether this set of strands have total multiplicity <= 1
        everywhere.

        """
        return all([n <= 1 for n in self.multiplicity])

def unconnectSumStrands(strands, genus1):
    """Returns pairs of strands (strand1, strand2) in (pmc1, pmc2), where
    (pmc1, pmc2) is the pair returned by unconnectSumPMC(strands.pmc, genus1).

    """
    cut_pos = 4 * genus1
    pmc1, pmc2 = unconnectSumPMC(strands.pmc, genus1)
    for p, q in strands:
        assert q < cut_pos or p >= cut_pos
    return (Strands(pmc1, [(p, q) for p, q in strands if q < cut_pos]),
            Strands(pmc2, [(p-cut_pos, q-cut_pos)
                           for p, q in strands if p >= cut_pos]))

class StrandDiagram(Generator):
    """Represents a strand diagram, or a generator of the strand algebra."""

    def __init__(self, parent, left_idem, strands, right_idem = None):
        """Specifies PMC, left idempotent and right idempotent as list of pair
        ID's, and strands as a list of pairs (start, end).
        For example, in the split PMC of genus 2, the strand diagram with
        double horizontal at (1,3) and strand from 2 to 5 would be encoded as:
        left_idem = [1,2], right_idem = [1,3], strands = [(2,5)], since pair
        (1,3) has index 1, pair (2,4) has index 2, and pair (5,7) has index 3.

        """
        Generator.__init__(self, parent)
        self.pmc = parent.pmc
        self.mult_one = parent.mult_one

        self.strands = strands
        if not isinstance(self.strands, Strands):
            self.strands = Strands(self.pmc, self.strands)

        # Calculate left idempotent if necessary
        if left_idem is None:
            assert right_idem is not None
            left_idem = self.strands.propagateLeft(right_idem)
        self.left_idem = left_idem
        if not isinstance(self.left_idem, Idempotent):
            self.left_idem = Idempotent(self.pmc, self.left_idem)

        # Calculate right idempotent if necessary
        if right_idem is None:
            right_idem = self.strands.propagateRight(self.left_idem)
        self.right_idem = right_idem
        if not isinstance(self.right_idem, Idempotent):
            self.right_idem = Idempotent(self.pmc, self.right_idem)

        # Enumerate double horizontals
        self.double_hor = list(self.left_idem)
        for st in self.strands:
            self.double_hor.remove(self.pmc.pairid[st[0]])
        self.double_hor = tuple(self.double_hor)

        # Get multiplicity from strands
        self.multiplicity = self.strands.multiplicity

    def getBigGrading(self):
        if not hasattr(self, "big_gr"):
            self.big_gr = self.pmc.big_gr(self.maslov(), self.multiplicity)
        return self.big_gr

    def getSmallGrading(self):
        if not hasattr(self, "small_gr"):
            refine_data = DEFAULT_REFINEMENT(self.pmc, len(self.left_idem))
            p_l, p_r = [refine_data[i]
                        for i in (self.left_idem, self.right_idem)]
            self.small_gr = (p_l * self.getBigGrading() * p_r.inverse()). \
                toSmallGrading()
        return self.small_gr

    def getGrading(self):
        if DEFAULT_GRADING == BIG_GRADING:
            return self.getBigGrading()
        else: # DEFAULT_GRADING == SMALL_GRADING
            return self.getSmallGrading()

    def __eq__(self, other):
        return self.parent == other.parent \
            and self.left_idem == other.left_idem \
            and self.strands == other.strands

    def __ne__(self, other):
        return not (self == other)

    @memorizeHash
    def __hash__(self):
        return hash((self.parent, self.left_idem, self.strands))

    def __str__(self):
        return "[%s]" % \
            ",".join([str(self.pmc.pairs[i]) for i in self.double_hor] +
                     ["%s->%s" % (p, q) for (p, q) in self.strands])

    def __repr__(self):
        return str(self)

    def isIdempotent(self):
        """Tests whether this generator is an idempotent."""
        return len(self.strands) == 0

    @memorize
    def opp(self):
        """Returns the opposite strand diagram in the opposite strand
        algebra.

        """
        return StrandDiagram(self.parent.opp(), self.right_idem.opp(),
                             self.strands.opp(), self.left_idem.opp())

    def numCrossing(self):
        """Returns the number of crossings between moving strands."""
        return sum(1 for (s1, t1) in self.strands for (s2, t2) in self.strands
                   if s1 < s2 and t1 > t2)

    def maslov(self):
        """Returns the Maslov index, defined as i(a) = inv(a) - m([a],S)."""
        maslov = Fraction()
        for s, t in self.strands:
            maslov -= Fraction(self.multiplicity[s], 2)
            if s != 0:
                maslov -= Fraction(self.multiplicity[s-1], 2)
        maslov += self.numCrossing()
        return maslov

    def getLeftIdem(self):
        """Return the left idempotent."""
        return self.left_idem

    def getRightIdem(self):
        """Return the right idempotent."""
        return self.right_idem

    def antiDiff(self):
        """Computes the dual of differential. See corresponding function in
        parent for details.

        """
        return self.parent.antiDiff(self)

    def factor(self):
        """Find all ways to factor this generator into a product of two
        generators. See corresponding function in parent for details.

        """
        return self.parent.factor(self)

def unconnectSumStrandDiagram(sd, genus1):
    """Returns a pair of strand diagrams (sd1, sd2) in the algebra of
    (pmc1, pmc2), where (pmc1, pmc2) is the pair returned by
    unconnectSumPMC(sd.pmc, genus1).

    """
    pmc1, pmc2 = unconnectSumPMC(sd.pmc, genus1)
    left_idem1, left_idem2 = unconnectSumIdem(sd.left_idem, genus1)
    strands1, strands2 = unconnectSumStrands(sd.strands, genus1)
    alg1 = StrandAlgebra(F2, pmc1, len(left_idem1), sd.mult_one)
    alg2 = StrandAlgebra(F2, pmc2, len(left_idem2), sd.mult_one)
    return (StrandDiagram(alg1, left_idem1, strands1),
            StrandDiagram(alg2, left_idem2, strands2))

class StrandAlgebra(DGAlgebra):
    """Represents the strand algebra of a PMC."""

    def __init__(self, ring, pmc, idem_size, mult_one = False):
        """Specifies the PMC, size of idempotent, and whether this is a
        multiplicity one algebra.

        """
        DGAlgebra.__init__(self, ring)
        self.pmc = pmc
        self.idem_size = idem_size
        self.mult_one = mult_one

    def __str__(self):
        return "Strand algebra over %s with idem_size = %d and mult_one = %r" \
            % (str(self.pmc), self.idem_size, self.mult_one)

    def __eq__(self, other):
        if not isinstance(other, StrandAlgebra):
            return False
        return self.pmc == other.pmc and self.idem_size == other.idem_size \
            and self.mult_one == other.mult_one

    def __ne__(self, other):
        return not (self == other)

    def __hash__(self):
        return hash((self.pmc, self.idem_size, self.mult_one))

    def opp(self):
        """Returns the opposite algebra, as the strand algebra associated to
        the opposite PMC.

        """
        return StrandAlgebra(self.ring, self.pmc.opp(), self.idem_size,
                             self.mult_one)

    @memorize
    def diff(self, gen):
        target_maslov = gen.maslov() - 1
        cur_strands = gen.strands
        result = [E0]
        def appendCandidate(new_strands):
            # Same info except strands, then check grading
            diff_term = StrandDiagram(self, gen.left_idem, new_strands,
                                      gen.right_idem)
            if diff_term.maslov() == target_maslov:
                result[0] += diff_term.elt()

        # Uncross two moving strands
        for st1 in cur_strands:
            for st2 in cur_strands:
                if st1[0] < st2[0] and st1[1] > st2[1]:
                    new_strands = list(cur_strands)
                    new_strands.remove(st1)
                    new_strands.remove(st2)
                    new_strands.extend([(st1[0], st2[1]), (st2[0], st1[1])])
                    appendCandidate(new_strands)

        # Uncross a moving strand with a double horizontal
        for st in cur_strands:
            for i in gen.double_hor:
                for p in gen.pmc.pairs[i]:
                    if st[0] <= p and p <= st[1]:
                        new_strands = list(cur_strands)
                        new_strands.remove(st)
                        new_strands.extend([(st[0], p), (p, st[1])])
                        appendCandidate(new_strands)
        return result[0]

    @memorize
    def getGenerators(self):
        if not self.mult_one:
            return self.pmc.getStrandDiagrams()
        else:
            return self.pmc.getMultOneStrandDiagrams()

    @memorize
    def getIdempotents(self):
        """Returns the set of idempotents. Use corresponding function in PMC.

        """
        return self.pmc.getIdempotents()

    @memorize
    def multiply(self, gen1, gen2):
        if not isinstance(gen1, StrandDiagram):
            return NotImplemented
        if not isinstance(gen2, StrandDiagram):
            return NotImplemented
        assert gen1.parent == self and gen2.parent == self, \
            "Algebra not compatible."
        if gen1.right_idem != gen2.left_idem:
            return E0
        if self.mult_one:
            # Enforce the multiplicity one condition
            total_mult = [m1+m2 for m1, m2 in zip(gen1.multiplicity,
                                                  gen2.multiplicity)]
            if not all([x <= 1 for x in total_mult]):
                return E0

        pmc = gen1.pmc
        new_strands = []

        # Keep track of which strands at right are not yet used.
        strands_right = list(gen2.strands)
        for sd in gen1.strands:
            mid_idem = pmc.pairid[sd[1]]
            possible_match = [sd2 for sd2 in strands_right
                              if pmc.pairid[sd2[0]] == mid_idem]
            if len(possible_match) == 0:
                new_strands.append(sd)
            else: # len(possible_match) == 1
                sd2 = possible_match[0]
                if sd2[0] != sd[1]:
                    return E0
                else:
                    new_strands.append((sd[0], sd2[1]))
                    strands_right.remove(sd2)

        new_strands.extend(strands_right)
        mult_term = StrandDiagram(self, gen1.left_idem, new_strands,
                                  gen2.right_idem)
        if mult_term.getBigGrading() == \
                gen1.getBigGrading() * gen2.getBigGrading():
            return mult_term.elt()
        else:
            return E0

    @memorize
    def _getAntiDiffMap(self):
        """Helper function generating tables of dual of differential, for calls
        to antiDiff.

        """
        antiDiffMap = {}
        for gen in self.getGenerators():
            antiDiffMap[gen] = E0
        for gen in self.getGenerators():
            for dgen, coeff in gen.diff().items():
                antiDiffMap[dgen] += coeff * gen
        return antiDiffMap

    def antiDiff(self, gen):
        """Returns the dual of the differential of gen, as an element of this
        algebra. The element is a sum of all terms c*y, for which gen appears
        in the dy with coefficient c.

        """
        return self._getAntiDiffMap()[gen]

    @memorize
    def _getFactorMap(self):
        """Helper function generating tables for factoring generators, for
        calls to factor.

        """
        factorMap = {}
        for gen in self.getGenerators():
            factorMap[gen] = E0
        # Resulting element lies in the tensor product of A with itself
        parent = Tensor((self, self))
        for gen1 in self.getGenerators():
            for gen2 in self.getGenerators():
                for prod, coeff in (gen1*gen2).items():
                    tensor_gen = TensorGenerator((gen1, gen2), parent)
                    factorMap[prod] += coeff * tensor_gen
        return factorMap

    def factor(self, gen):
        """Returns an element of (A tensor A), where A is the present algebra.
        The element is the sum of all terms c*(p tensor q), for which gen is a
        term in p*q with coefficient c.

        """
        return self._getFactorMap()[gen]

class StrandAlgebraElement(Element):
    """An element of strand algebra."""
    def isIdempotent(self):
        """Tests whether this element is an idempotent."""
        for sd, coeff in self.items():
            if not sd.isIdempotent():
                return False
        return True

    def invertible(self):
        """Tests whether this element is invertible."""
        return self != 0 and self.isIdempotent()

    def inverse(self):
        """Returns the inverse of this element, if invertible. Undefined
        behavior if the element is not invertible.

        """
        return self

StrandDiagram.ELT_CLASS = StrandAlgebraElement
