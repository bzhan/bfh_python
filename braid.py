"""Handles braids and their type DD structures."""

import sys
from arcslide import Arcslide
from arcslideda import ArcslideDA
from dehntwistda import DehnSurgeryDA
from digraph import computeATensorD, computeATensorDD, computeDATensorD
from dstructure import infTypeD, platTypeD
from pmc import linearPMC, splitPMC
from utility import memorize
from utility import NEG, POS, PRINT_PROGRESS

class Braid():
    """Represents a braid with a fix number of strands. Each braid generator is
    represented by an integer n. If 1 <= n <= num_strands-1, it represents
    moving strand n over strand n+1. If -(num_strand-1) <= n <= -1, it
    represents moving strand |n| under strand |n|+1.

    """
    def __init__(self, num_strands):
        """Specifies the number of strands in the braid."""
        self.num_strands = num_strands
        self.genus = (num_strands - 2)/2
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
        print "(compose %d)" % len(ddstr_list),
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
        print "%d\n" % len(dstr),
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

class BraidCap:
    """Represents a capping of a braid."""
    def __init__(self, matching):
        """Specifies the matching of strands."""
        self.matching = tuple(matching)
        self.num_strands = len(self.matching)
        self.genus = (len(self.matching) - 2)/2
        self.fix = self.ends[self.matching]

    def __eq__(self, other):
        return self.matching == other.matching

    def __ne__(self, other):
        return not (self == other)

    def __hash__(self):
        return hash(self.matching)

    @memorize
    def getHandlebody(self, is_dual = False):
        """Returns the handlebody corresponding to this matching."""
        if PRINT_PROGRESS:
            print "Get handlebody for %s" % str(self.matching)
        slides = Braid(self.num_strands).getArcslides(self.fix)
        if not is_dual:
            slides = [slide.inverse() for slide in slides]
        slides_dd = [slide.getDDStructure() for slide in slides]
        dstr = composeDD(platTypeD2(self.genus, is_dual), slides_dd, is_dual)
        dstr.checkGrading()
        return dstr

    @memorize
    def getHandlebodyByLocalDA(self):
        """Returns the handlebody corresponding to this matching. Use local type
        DA structures associated to arcslides. Note this corresponds to
        dual = False case of getHandlebody().

        """
        slides = Braid(self.num_strands).getArcslides(self.fix)
        slides = [slide.inverse() for slide in slides]
        dstr = platTypeD(self.genus)
        for slide in slides:
            dstr = ArcslideDA(slide).tensorD(dstr)
            dstr.reindex()
            dstr.simplify()
        return dstr

    ends = {
        # Genus 3
        (6,5,4,3,2,1) : [3,4],
        (6,3,2,5,4,1) : [],
        (4,3,2,1,6,5) : [1,2,3,4,2,3],
        (2,1,6,5,4,3) : [1,2],
        (2,1,4,3,6,5) : [1,2,3,4],
        # Genus 4
        (8,7,6,5,4,3,2,1) : [3,4,5,4,3,2],
        (8,7,4,3,6,5,2,1) : [3,4,5,6],
        (8,5,4,3,2,7,6,1) : [3,4],
        (8,3,2,7,6,5,4,1) : [5,6],
        (8,3,2,5,4,7,6,1) : [],
        (6,5,4,3,2,1,8,7) : [1,2,3,4,5,6,2,3,4,5,3,4],
        (6,3,2,5,4,1,8,7) : [1,2,3,4,5,6,2,3,4,5],
        (4,3,2,1,8,7,6,5) : [1,2,3,4,2,3],
        (4,3,2,1,6,5,8,7) : [1,2,3,4,5,6,2,3],
        (2,1,8,7,6,5,4,3) : [1,2,5,6],
        (2,1,8,5,4,7,6,3) : [1,2],
        (2,1,6,5,4,3,8,7) : [1,2,3,4,5,6,4,5],
        (2,1,4,3,8,7,6,5) : [1,2,3,4],
        (2,1,4,3,6,5,8,7) : [1,2,3,4,5,6],
        # Genus 5 (incomplete)
        (2,1,4,3,6,5,8,7,10,9) : [1,2,3,4,5,6,7,8],
        (2,1,10,7,6,5,4,9,8,3) : [1,2,5,6],
        (4,3,2,1,6,5,10,9,8,7) : [1,2,3,4,5,6,2,3],
        (4,3,2,1,10,9,8,7,6,5) : [1,2,3,4,2,3,7,8],
        (6,5,4,3,2,1,8,7,10,9) : [1,2,3,4,5,6,7,8,2,3,4,5,3,4],
        (6,5,4,3,2,1,10,9,8,7) : [1,2,3,4,5,6,2,3,4,5,3,4],
        (8,3,2,7,6,5,4,1,10,9) : [1,2,3,4,5,6,7,8,2,3,4,5,6,7,5,6],
        (8,5,4,3,2,7,6,1,10,9) : [1,2,3,4,5,6,7,8,2,3,4,5,6,7,3,4],
        (10,3,2,5,4,7,6,9,8,1) : [],
        (10,3,2,7,6,5,4,9,8,1) : [5,6],
        (10,9,8,7,6,5,4,3,2,1) : [3,4,5,6,7,8,4,5,6,7,5,6],
        # Genus 6 (incomplete)
        (12,3,2,5,4,7,6,9,8,11,10,1) : [],
        (12,11,10,9,8,7,6,5,4,3,2,1) :
            [3,4,5,6,7,8,9,10,4,5,6,7,8,9,5,6,7,8,6,7],
        }

    @staticmethod
    def verifyEnds():
        """Verify the corrections for ends below."""
        for matching, fix in BraidCap.ends.items():
            assert len(matching) % 2 == 0
            num_bridge = len(matching) / 2
            # Form starting matching
            cur_match = [num_bridge*2]
            for i in range(1, num_bridge):
                cur_match += [2*i+1, 2*i]
            cur_match += [1]
            # Now perform exchange for each step
            for i in range(len(fix)):
                # Swaps points fix[i] and fix[i]+1, which are stored at IDs
                # fix[i]-1 and fix[i] (a bit awkward).
                id1, id2 = fix[i]-1, fix[i]
                pt1, pt2 = fix[i], fix[i]+1
                assert cur_match[id1] != pt2 and cur_match[id2] != pt1
                cur_match[id1], cur_match[id2] = cur_match[id2], cur_match[id1]
                for j in range(len(cur_match)):
                    if cur_match[j] == pt1:
                        cur_match[j] = pt2
                    elif cur_match[j] == pt2:
                        cur_match[j] = pt1
            assert matching == tuple(cur_match)


class BridgePresentation:
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

    def getHF(self):
        """Computes HF of branched double cover."""
        start_d = BraidCap(self.start).getHandlebody(False)
        slides = Braid(self.num_strands).getArcslides(self.braid_word)
        # Use this line only for "Tensor"
        # slides = [slide.inverse() for slide in slides]
        slides_dd = [slide.getDDStructure() for slide in slides]
        end_d = BraidCap(self.end).getHandlebody(True)
        start_d = composeDD(start_d, slides_dd, is_dual = False, method = "Mor")
        cx = computeATensorD(end_d, start_d)
        cx.simplify()
        cx.reindex()
        return cx

    def getHFByLocalDA(self):
        """Compute HF of branched double cover, using local type DA structures
        for arcslides.

        """
        start_d = BraidCap(self.start).getHandlebodyByLocalDA()
        slides = Braid(self.num_strands).getArcslides(self.braid_word)
        slides = [slide.inverse() for slide in slides]
        for slide in slides:
            print "%d" % len(start_d),
            sys.stdout.flush(),
            start_d = ArcslideDA(slide).tensorD(start_d)
            start_d.reindex()
            # assert start_d.testDelta()  # remove this when more confident
            start_d.simplify()
        # First way - limited by genus. No infinity issues.
        # end_d = BraidCap(self.end).getHandlebodyByLocalDA()
        # cx = start_d.morToD(end_d)
        # Second way - has infinity issues when diagram is not admissible
        end_d = BraidCap(self.end).getHandlebodyByLocalDA().dual()
        print " -> %d, %d" % (len(start_d), len(end_d))
        cx = computeATensorD(end_d, start_d)
        cx.reindex()
        cx.checkDifferential()
        cx.simplify()
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

        start_d = BraidCap(self.start).getHandlebodyByLocalDA()
        genus = self.num_strands/2 - 1
        for twist in self.braid_word:
            print "%d" % len(start_d),
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
        # Limited by genus, no infinity issues
        end_d = BraidCap(self.end).getHandlebodyByLocalDA()
        print " -> %d, %d" % (len(start_d), len(end_d))
        cx = start_d.morToD(end_d)
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
