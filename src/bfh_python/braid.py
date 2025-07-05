"""Handles braids and their type DD structures."""

import sys
from .arcslide import Arcslide
from .arcslideda import ArcslideDA
from .cobordism import Cobordism
from .cobordism import LEFT, RIGHT
from .cobordismda import CobordismDALeft, CobordismDARight
from .dehntwistda import DehnSurgeryDA
from .digraph import computeATensorD, computeATensorDD, computeDATensorD
from .dstructure import infTypeD, platTypeD, zeroTypeD
from .pmc import linearPMC, splitPMC
from .utility import memorize
from .utility import NEG, POS, PRINT_PROGRESS

class Braid(object):
    """Represents a braid with a fix number of strands. Each braid generator is
    represented by an integer n. If 1 <= n <= num_strands-1, it represents
    moving strand n over strand n+1. If -(num_strand-1) <= n <= -1, it
    represents moving strand |n| under strand |n|+1.

    """
    def __init__(self, num_strands):
        """Specifies the number of strands in the braid."""
        self.num_strands = num_strands
        self.genus = (num_strands - 2)//2
        self.pmc = linearPMC(self.genus)

    def getArcslides(self, word):
        """Get the sequence of arcslides corresponding to a braid generator or
        a list of braid generators.

        """
        if isinstance(word, list):
            return sum([self.getArcslides(n) for n in word], [])

        abs_word = abs(word)
        assert 1 <= abs_word <= self.num_strands-1
        if abs_word == 1:
            slides_info = [(1, 0)]
        elif 1 < abs_word < self.num_strands-2:
            slides_info = [(2*abs_word-2, 2*abs_word-3)] * 2
        elif abs_word == self.num_strands-2:
            slides_info = [(2*abs_word-2, 2*abs_word-3)]
        else: # abs_word = self.num_strands-1
            slides_info = [(4*i-3, 4*i-4) for i in range(self.genus, 0, -1)]
            slides_info += [(2*i+2, 2*i+1) for i in range(2*self.genus-2)]
        slides = [Arcslide(self.pmc, slides_info[0][0], slides_info[0][1])]
        for i in range(1, len(slides_info)):
            slides.append(Arcslide(slides[i-1].end_pmc,
                                   slides_info[i][0], slides_info[i][1]))
        if word > 0:
            return slides
        else: # n < 0
            return list(reversed([slide.inverse() for slide in slides]))

def composeDD(dstr, ddstr_list, is_dual = False, method = "Tensor"):
    """Successively compute morphisms from DD structures in ``ddstr_list`` to
    the type D structure, simplifying at each step.

    """
    if PRINT_PROGRESS > 0:
        print("(compose %s %d)" % (method, len(ddstr_list)), end=' ')
    for ddstr in ddstr_list:
        if PRINT_PROGRESS > 0:
            sys.stdout.write("%d," % len(dstr))
            sys.stdout.flush()
        if method == "Mor":
            assert not is_dual
            dstr = ddstr.morToD(dstr)
        elif method == "Tensor":
            if is_dual:
                dstr = computeATensorDD(dstr, ddstr)
            else:
                dstr = computeDATensorD(ddstr, dstr)
        else:
            assert False, "Unknown method"
        dstr.simplify()
        dstr.reindex()
    if PRINT_PROGRESS > 0:
        print("%d\n" % len(dstr), end=' ')
    return dstr

@memorize
def platTypeD2(genus, is_dual = False):
    """Obtain linear handlebody from inf-framed handlebody by a sequence of
    arcslides. As the inf-framed handlebody has absolute grading, we can get
    absolute grading on this D structure.

    """
    start = infTypeD(genus, is_dual)
    slides = []
    for i in range(genus-1):
        slides += [(4*i+3,4*i+4), (4*i+6,4*i+7), (4*i+5,4*i+6)]
    cur_pmc = splitPMC(genus)
    for i in range(len(slides)):
        b1, c1 = slides[i]
        slides[i] = Arcslide(cur_pmc, b1, c1)
        cur_pmc = slides[i].end_pmc
    if not is_dual:
        slides = [slide.inverse() for slide in slides]
    slides_dd = [slide.getDDStructure() for slide in slides]
    return composeDD(start, slides_dd, is_dual)

class BraidCap(object):
    """Represents a capping of a braid."""
    def __init__(self, matching):
        """Specifies the matching of strands."""
        self.matching = tuple(matching)
        self.num_strands = len(self.matching)
        self.genus = (len(self.matching) - 2)//2

    def __eq__(self, other):
        return self.matching == other.matching

    def __ne__(self, other):
        return not (self == other)

    def __hash__(self):
        return hash(self.matching)

    def getLastCobordism(self):
        """Returns the position of the leftmost pair of adjacent points matched
        under this braid cap. This is a pair that can be considered as last
        added in a sequence of cobordisms forming this braid cap. Also returns
        the braid cap before this cobordism.

        """
        if self.genus == 0:
            return None
        # At each step, find the left-most ending point of a matching.
        for i in range(len(self.matching)):
            if self.matching[i] <= i+1:
                assert self.matching[i] == i
                cur_move = i-1
                new_matching = []
                for n in self.matching[:i-1] + self.matching[i+1:]:
                    assert n != i and n != i+1
                    if n < i:
                        new_matching.append(n)
                    else:
                        new_matching.append(n-2)
                new_cap = BraidCap(new_matching)
                return (cur_move, new_cap)

    @memorize
    def getCobordismSequence(self):
        """Returns a sequence of cobordisms that will close off this braid cap
        to two strands (genus-0). The cobordisms are labeled 0-based starting
        from the left. The right-most cobordism is never used.

        For example:
        (4, 3, 2, 1) --> [1]
        (2, 1, 4, 3) --> [0]
        (6, 5, 4, 3, 2, 1) --> [2, 1]
        (6, 3, 2, 5, 4, 1) --> [1, 1]
        (2, 1, 6, 5, 4, 3) --> [0, 1]

        """
        if self.genus == 0:
            return []
        else:
            last_move, prev_cap = self.getLastCobordism()
            return [last_move] + prev_cap.getCobordismSequence()

    @memorize
    def openCap(self):
        """Returns the type D structure corresponding to this handlebody by
        tensoring type DA bimodules with the genus-1 handlebody.

        """
        assert self.genus > 0
        last_move, prev_cap = self.getLastCobordism()
        if self.genus == 1:
            if last_move == 0:
                return zeroTypeD(1)
            else:
                return infTypeD(1)
        else:
            cur_da = CobordismDALeft(Cobordism(self.genus, last_move, LEFT))
            dstr = cur_da.tensorD(prev_cap.openCap())
            dstr.simplify()
            dstr.reindex()
        return dstr

    def closeCap(self, dstr, cancellation_constraint = None):
        """Computes the chain complex obtained by closing off this cap on dstr.
        That is, compute the box tensor product of the type A module
        corresponding to this cap with dstr.

        This is obtained by tensoring dstr with a sequence of right-side
        cobordisms, and finishing off by computing morToD with either
        zeroTypeD(1) or infTypeD(1).

        """
        assert self.genus > 0
        if self.genus <= 3:
            # morToD is efficient up to genus = 3.
            cx = dstr.morToD(self.openCap())
            cx.reindex()
            cx.simplify(cancellation_constraint = cancellation_constraint)
            return cx
        else:
            last_move, prev_cap = self.getLastCobordism()
            cur_da = CobordismDARight(Cobordism(self.genus, last_move, RIGHT))
            dstr = cur_da.tensorD(
                dstr, cancellation_constraint = cancellation_constraint)
            dstr.reindex()
            dstr.simplify(cancellation_constraint = cancellation_constraint)
            return prev_cap.closeCap(dstr, cancellation_constraint)

class BridgePresentation(object):
    """Represents a bridge presentation of a knot. Computes HF of branched
    double cover from bridge presentation.

    """
    def __init__(self, name, start, braid_word, end):
        """Specifies start (of type BraidCap), braid_word (of type list), and
        end (of type BraidCap).

        """
        self.name = name
        self.num_strands = len(start)
        self.start = start
        self.braid_word = braid_word
        self.end = end

    def getHF(self, method = "Mor"):
        """Computes HF of branched double cover."""
        assert method in ("Mor", "Tensor")
        start_d = BraidCap(self.start).openCap()
        slides = Braid(self.num_strands).getArcslides(self.braid_word)
        if method == "Tensor":
            slides = [slide.inverse() for slide in slides]
        slides_dd = [slide.getDDStructure() for slide in slides]
        start_d = composeDD(
            start_d, slides_dd, is_dual = False, method = method)
        if method == "Tensor":
            end_d = BraidCap(self.end).openCap().dual()
            cx = computeATensorD(end_d, start_d)
        else:
            end_d = BraidCap(self.end).openCap()
            cx = end_d.morToD(start_d)
        cx.simplify()
        cx.reindex()
        return cx

    def getHFByLocalDA(self):
        """Compute HF of branched double cover, using local type DA structures
        for arcslides.

        """
        start_d = BraidCap(self.start).openCap()
        slides = Braid(self.num_strands).getArcslides(self.braid_word)
        slides = [slide.inverse() for slide in slides]
        for slide in slides:
            print("%d" % len(start_d), end=' ')
            sys.stdout.flush(),
            start_d = ArcslideDA(slide).tensorD(start_d)
            start_d.reindex()
            start_d.simplify()
        # Close off using cobordisms
        cx = BraidCap(self.end).closeCap(start_d)
        return cx

    def addStrandsAtRight(self):
        """Return a bridge presentation for the same knot, with two more strands
        at right that are not involved in any crossings. The new bridge
        presentation will always be acceptable to getSpecSeq below.

        """
        n = self.num_strands
        new_start = self.start + [n+2, n+1]
        new_end = self.end + [-1, -1]
        for i in range(len(new_end)):
            if new_end[i] == n:
            # Strand i+1 is matched with n+2
                new_end[i] = n+2
                new_end[n+1] = i+1
        # Strand n matched with n+1 (note indices are 0-based).
        new_end[n-1] = n+1
        new_end[n] = n
        return BridgePresentation(
            self.name, new_start, self.braid_word, new_end)

    def getSpecSeq(self):
        """Compute the spectral sequence from Khovanov homology to HF of
        branched double cover, using local type DA structures for Dehn twists
        (as mapping cone between identity and anti-braid).

        Returns a list of lists. The first list contains counts of filtration
        gradings of generators in the Khovanov homology, starting from the one
        with lowest to the one with the highest filtration grading. The next
        list contains counts of filtration gradings in the E_3 page, etc, until
        the pages have stabilized.

        """
        if self.num_strands-2 in self.braid_word:
            return self.addStrandsAtRight().getSpecSeq()

        start_d = BraidCap(self.start).openCap()
        genus = self.num_strands//2 - 1
        for twist in self.braid_word:
            print("%d" % len(start_d), end=' ')
            sys.stdout.flush()
            abs_twist = abs(twist)
            assert 1 <= abs_twist <= self.num_strands-2
            # Choice of orientation for the knot
            if twist < 0:
                surgery = DehnSurgeryDA(genus, abs_twist-1, POS)
            else:
                surgery = DehnSurgeryDA(genus, abs_twist-1, NEG)
            surgery_da = surgery.getMappingCone()
            start_d = surgery_da.tensorD(start_d)
            start_d.reindex()
            # Must be done in two steps
            start_d.simplify(cancellation_constraint = lambda x, y: (
                sum(x.filtration) == sum(y.filtration)))
            start_d.simplify(cancellation_constraint = lambda x, y: (
                sum(x.filtration) + 1 >= sum(y.filtration)))
        # Must not simplify everything immediately.
        cx = BraidCap(self.end).closeCap(
            start_d, cancellation_constraint = lambda x, y: (
                sum(x.filtration) == sum(y.filtration)))
        cx.reindex()
        cx.checkDifferential()
        cx.simplify(cancellation_constraint = lambda x, y: (
            sum(x.filtration) == sum(y.filtration)))

        result = []
        filt_diff = 1
        while any(x.diff() != 0 for x in cx.getGenerators()):
            cx.simplify(cancellation_constraint = lambda x, y: (
                sum(x.filtration) + filt_diff >= sum(y.filtration)))
            filt_grs = [sum(gen.filtration) for gen in cx.getGenerators()]
            if filt_diff == 1:
                # Find minimum and maximum at the second page
                min_filt, max_filt = min(filt_grs), max(filt_grs)
            result.append([filt_grs.count(i)
                           for i in range(min_filt, max_filt+1)])
            filt_diff += 1
        return result

    def __str__(self):
        return str(self.name)

    def __repr__(self):
        result = "Bridge presentation for %s:\n" % str(self.name)
        result += "Start = %s\n" % str(self.start)
        result += "Braid word = %s\n" % str(self.braid_word)
        result += "End = %s\n" % str(self.end)
        return result

def readBridgePresentation(str_input):
    """Read bridge presentation from string input. The format is as follows:

    str_input is a line with space-separated tokens. The first token is the
    name of the knot (can be anything that doesn't contain a space). The
    remaining tokens are integers. The first integer is the integer k/2, where
    k is the number of strands in the braid.

    The next k integers 0 <= a_0, a_1, ... a_{k-1} < k specify the top closure
    of the braid. That is, the i'th strand is paired with the (a_i)'th strand
    (so we must have a_{a_i} = i).

    The next integer is the number of crossings n. This is followed by n pairs
    of integers, with each pair i, j meaning moving the i'th strand over the
    j'th strand (so |i - j| = 1).

    The next k integers specify the bottom closure of the braid, in the same
    format as that for the top closure.

    """
    tokens = str_input.split()
    name = tokens[0]
    rest = [int(token) for token in tokens[1:]]
    bridge_size = rest[0]
    start = [1+n for n in rest[1:1+bridge_size*2]]
    rest = rest[1+bridge_size*2:]
    num_cross = rest[0]
    braid_word = []
    for i in range(num_cross):
        p, q = rest[2*i+1], rest[2*i+2]
        if p == q-1:
            braid_word.append(q)
        elif p == q+1:
            braid_word.append(-p)
        else:
            assert False
    rest = rest[1+num_cross*2:]
    end = [1+n for n in rest]
    return BridgePresentation(name, start, braid_word, end)
