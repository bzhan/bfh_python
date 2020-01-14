"""Unit test for digraph.py"""

from digraph import *
from arcslide import Arcslide
from dstructure import infTypeD, platTypeD, zeroTypeD, zeroTypeDAdm
from ddstructure import identityDD
from pmc import PMC
from pmc import splitPMC
from utility import DEFAULT_GRADING, SMALL_GRADING
import unittest

class TypeDGraphTest(unittest.TestCase):
    def testBuildTypeDGraph(self):
        graph1 = TypeDGraph(platTypeD(2))

class TypeDDGraphTest(unittest.TestCase):
    def testBuildTypeDDGraph(self):
        graph1 = TypeDDGraph(identityDD(splitPMC(1)), 2)

class TypeAAGraphTest(unittest.TestCase):
    def testBuildTypeAAGraph(self):
        graph1 = TypeAAGraph(splitPMC(1))

    def testTensorDoubleD(self):
        d1 = zeroTypeD(1)
        d2 = zeroTypeDAdm(1)
        d3 = zeroTypeD(2)
        d4 = zeroTypeDAdm(2)
        d5 = zeroTypeD(3)
        d6 = zeroTypeDAdm(3)
        d7 = zeroTypeD(4)
        d8 = zeroTypeDAdm(4)
        tests = [(d1, d2, 2), (d2, d1, 2), (d2, d2, 2),
                 (d3, d4, 4), (d4, d3, 4), (d4, d4, 4),
                 (d5, d6, 8), (d7, d8, 16)]
        for d_left, d_right, expected_len in tests:
            cx = computeATensorD(d_left.dual(), d_right)
            cx.simplify()
            self.assertEqual(len(cx), expected_len)

    def testTensorDDandD(self):
        d1 = zeroTypeD(1)
        d2 = zeroTypeDAdm(1)
        d3 = zeroTypeD(2)
        dd_id = identityDD(splitPMC(1))
        dd_id2 = identityDD(splitPMC(2))
        dd_slide1 = Arcslide(splitPMC(1), 0, 1).getDDStructure()
        dstr1 = computeDATensorD(dd_id, d1)
        dstr2 = computeDATensorD(dd_id, d2)
        dstr3 = computeDATensorD(dd_slide1, d1)
        dstr4 = computeDATensorD(dd_id2, d3)
        # Uncomment to see the structures
        # print dstr1, dstr2, dstr3, dstr4

    def testTensorDandDD(self):
        d1 = zeroTypeD(1)
        dd_id = identityDD(splitPMC(1))
        dstr1 = computeATensorDD(d1, dd_id)
        # Uncomment to see the structures
        # print dstr1

    def testTensorDoubleDD(self):
        dd_id = identityDD(splitPMC(1))
        self.assertTrue(computeDATensorDD(dd_id, dd_id).testDelta())
        dd_id2 = identityDD(splitPMC(2))
        self.assertTrue(computeDATensorDD(dd_id2, dd_id2).testDelta())

        def composeSlides(start_pmc, slides):
            cur_dd = Arcslide(start_pmc, *slides[0]).getDDStructure()
            for slide in slides[1:]:
                next_dd = Arcslide(
                    cur_dd.algebra2.pmc.opp(), *slide).getDDStructure()
                cur_dd = computeDATensorDD(cur_dd, next_dd)
                cur_dd.simplify()
                cur_dd.reindex()
            return cur_dd

        tests = [(splitPMC(1), [(0,1),(3,2)], 2),
                 (PMC([(0,3),(1,6),(2,4),(5,7)]), [(2,1),(6,5),(4,3)], 6)]
        for start_pmc, slides, result in tests:
            composed_dd = composeSlides(start_pmc, slides)
            self.assertEqual(len(composed_dd), result)
            composed_dd.checkGrading()
            if DEFAULT_GRADING == SMALL_GRADING:
                ref_gr = list(composed_dd.grading.values())[0]
                for gen, gr in list(composed_dd.grading.items()):
                    self.assertEqual(gr, ref_gr)

    def testGrading(self):
        d1 = zeroTypeD(1)
        d2 = infTypeD(1)
        dd_id = identityDD(splitPMC(1))
        cx = computeATensorD(d1.dual(), d2)
        dstr1 = computeDATensorD(dd_id, d1)
        dstr2 = computeATensorDD(d1, dd_id)

    def testAATensorDD1(self):
        pmc = splitPMC(2)
        aa_graph = getTypeAAGraph(pmc)
        ddstr1 = identityDD(pmc)
        dd_graph1 = TypeDDGraph(ddstr1, 1)
        dastr1 = aa_graph.tensorAAandDD(dd_graph1)
        self.assertTrue(dastr1.testDelta())

    def testAATensorDD2(self):
        pmc2 = PMC([(0,3),(1,6),(2,4),(5,7)])
        ddstr2 = Arcslide(pmc2, 0, 1).getDDStructure()
        dd_graph2 = TypeDDGraph(ddstr2, 1)
        aa_graph2 = getTypeAAGraph(pmc2)
        dastr2 = aa_graph2.tensorAAandDD(dd_graph2)
        self.assertTrue(dastr2.testDelta())

    def testDDTensorAA(self):
        pmc = splitPMC(1)
        aa_graph = getTypeAAGraph(pmc)
        ddstr1 = identityDD(pmc)
        dd_graph1 = TypeDDGraph(ddstr1, 2)
        dastr1 = aa_graph.tensorDDandAA(dd_graph1)
        self.assertTrue(dastr1.testDelta())

        ddstr2 = Arcslide(pmc, 0, 1).getDDStructure()
        dd_graph2 = TypeDDGraph(ddstr2, 2)
        dastr2 = aa_graph.tensorDDandAA(dd_graph2)
        self.assertTrue(dastr2.testDelta())

if __name__ == "__main__":
    unittest.main()
