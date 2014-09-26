"""Unit test for cobordismda.py"""

from cobordismda import *
from ddstructure import identityDD
from pmc import splitPMC
import unittest

class CobordismDATest(unittest.TestCase):
    def testLeftCobordismDA(self):
        for genus, c_pair in [(2, 0), (2, 1), (2, 2), (2, 3),
                              (3, 0), (3, 1), (3, 2), (3, 3), (3, 4), (3, 5)]:
            c = Cobordism(genus, c_pair, LEFT)
            c_da = CobordismDALeft(c)
            dastr = c_da.toSimpleDAStructure()
            self.assertTrue(dastr.testDelta())
            ddstr = dastr.tensorDD(identityDD(c.end_pmc))
            ori_ddstr = c.getDDStructure()
            self.assertTrue(ddstr.compareDDStructures(ori_ddstr))

    def testRightCobordismDA(self):
        for genus, c_pair in [(2, 0), (2, 1), (2, 2), (2, 3),
                              (3, 0), (3, 1), (3, 2), (3, 3), (3, 4), (3, 5)]:
            c = Cobordism(genus, c_pair, RIGHT)
            c_da = CobordismDARight(c)
            ddstr = c_da.tensorDD(identityDD(c.end_pmc))
            ddstr.simplify()
            ori_ddstr = c.getDDStructure()
            self.assertTrue(ddstr.compareDDStructures(ori_ddstr))

class SimpleCobordismDATest(unittest.TestCase):
    def testSimpleCobordismDA(self):
        for pmc, insert_pos in [
                (splitPMC(1), 3),
                (splitPMC(1), 1),
        ]:
            c_da = SimpleCobordismDA(pmc, insert_pos)
            dastr = c_da.toSimpleDAStructure()
            self.assertTrue(dastr.testDelta())
            ddstr = dastr.tensorDD(identityDD(c_da.end_pmc))
            ddstr.simplify()
            ddstr.reindex()
            self.assertEquals(len(ddstr.getGenerators()), 2)
            self.assertEquals(
                sorted(len(gen.delta()) for gen in ddstr.getGenerators()),
                [2, 6])

    def slideSeq(self, start_pmc, slides):
        # Given a sequence of arcslides tau_1, ... tau_n specified by the
        # starting PMC of tau_1 and a list of (b1, c1), compute:
        # CFDA(tau_1) * ... CFDA(tau_n) * CFDD(Id).
        assert len(slides) > 0
        slides_da = []
        for b1, c1 in slides:
            # Find the list CFDA(tau_i)
            slides_da.append(ArcslideDA(Arcslide(start_pmc, b1, c1)))
            start_pmc = slides_da[-1].pmc2

        dd = identityDD(slides_da[-1].pmc2)
        for slide_da in reversed(slides_da):
            dd = slide_da.tensorDD(dd)
            dd.reindex()
            dd.simplify()
        return dd

    def testComposeMiddle(self):
        c_da = SimpleCobordismDA(splitPMC(1), 1)
        dd = c_da.tensorDD(self.slideSeq(
            c_da.pmc2, [(5, 4), (2, 1), (3, 2), (6, 5)]))
        dd.simplify()
        c = Cobordism(2, 1, RIGHT)
        ori_dd = c.getDDStructure()
        self.assertTrue(dd.compareDDStructures(ori_dd))

    def testComposeNextTop(self):
        c_da = SimpleCobordismDA(splitPMC(1), 3)
        dd = c_da.tensorDD(self.slideSeq(c_da.pmc2, [(7, 6)]))
        dd.simplify()
        c = Cobordism(2, 2, RIGHT)
        ori_dd = c.getDDStructure()
        self.assertTrue(dd.compareDDStructures(ori_dd))

    def testComposeTop(self):
        c_da = SimpleCobordismDA(splitPMC(1), 3)
        dd = c_da.tensorDD(self.slideSeq(c_da.pmc2, [(6, 5), (7, 6)]))
        dd.simplify()
        c = Cobordism(2, 3, RIGHT)
        ori_dd = c.getDDStructure()
        self.assertTrue(dd.compareDDStructures(ori_dd))

    def testComposeBottom(self):
        c_da = SimpleCobordismDA(splitPMC(1), 1)
        dd = c_da.tensorDD(self.slideSeq(c_da.pmc2, [(0, 1)]))
        dd.simplify()
        c = Cobordism(2, 0, RIGHT)
        ori_dd = c.getDDStructure()
        self.assertTrue(dd.compareDDStructures(ori_dd))

if __name__ == "__main__":
    unittest.main()
