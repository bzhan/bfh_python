"""Handles braids and their type DD structures."""

import sys
from arcslide import Arcslide
from arcslideda import ArcslideDA
from digraph import computeATensorD, computeATensorDD, computeDATensorD
from dstructure import infTypeD, platTypeD
from pmc import linearPMC, splitPMC
from utility import memorize
from utility import PRINT_PROGRESS

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
        print "\n",
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
        assert len(matching) in (6, 8)
        self.matching = tuple(matching)
        self.num_strands = len(self.matching)
        self.genus = (len(self.matching) - 2)/2
        if len(matching) == 6:
            self.index = self.ends3.index(self.matching)
            self.fix = self.ends3_fix[self.index]
        else: # len(matching) == 8:
            self.index = self.ends4.index(self.matching)
            self.fix = self.ends4_fix[self.index]

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
        assert dstr.testDelta()
        return dstr

    ends3 = [(6,5,4,3,2,1),(6,3,2,5,4,1),(4,3,2,1,6,5),(2,1,6,5,4,3),
             (2,1,4,3,6,5)]
    ends4 = [(8,7,6,5,4,3,2,1),(8,7,4,3,6,5,2,1),(8,5,4,3,2,7,6,1),
             (8,3,2,7,6,5,4,1),(8,3,2,5,4,7,6,1),(6,5,4,3,2,1,8,7),
             (6,3,2,5,4,1,8,7),(4,3,2,1,8,7,6,5),(4,3,2,1,6,5,8,7),
             (2,1,8,7,6,5,4,3),(2,1,8,5,4,7,6,3),(2,1,6,5,4,3,8,7),
             (2,1,4,3,8,7,6,5),(2,1,4,3,6,5,8,7)]
    ends3_fix = [[3,4],[],[1,2,3,4,2,3],[1,2],[1,2,3,4]]
    ends4_fix = [[3,4,5,4,3,2],[3,4,5,6],[3,4],[5,6],[],
                 [1,2,3,4,5,6,2,3,4,5,3,4],[1,2,3,4,5,6,2,3,4,5],[1,2,3,4,2,3],
                 [1,2,3,4,5,6,2,3],[1,2,5,6],[1,2],[1,2,3,4,5,6,4,5],[1,2,3,4],
                 [1,2,3,4,5,6]]

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
        slides = [slide.inverse() for slide in slides]
        slides_dd = [slide.getDDStructure() for slide in slides]
        end_d = BraidCap(self.end).getHandlebody(True)
        start_d = composeDD(start_d, slides_dd, is_dual = False)
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
            start_d = ArcslideDA(slide).tensorD(start_d)
            start_d.reindex()
            assert start_d.testDelta()  # remove this when more confident
            start_d.simplify()
        end_d = BraidCap(self.end).getHandlebodyByLocalDA().dual()
        cx = computeATensorD(end_d, start_d)
        cx.reindex()
        cx.checkDifferential()
        cx.simplify()
        return cx

    def __str__(self):
        return str(self.name)

    def __repr__(self):
        result = "Bridge presentation for %s:\n" % str(self.name)
        result += "Start = %s\n" % str(self.start)
        result += "Braid word = %s\n" % str(self.braid_word)
        result += "End = %s\n" % str(self.end)
        return result

def readBridgePresentation(str_input):
    """Read bridge presentation from string input."""
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
