"""Unit test for ddstructure.py"""

from ddstructure import *
from dstructure import infTypeD, zeroTypeD
from pmc import PMC
from pmc import linearPMC, splitPMC
from utility import DEFAULT_GRADING, SMALL_GRADING
import unittest

class DDStructureTest(unittest.TestCase):
    def testCommonIdentityDD(self):
        id2 = identityDD(splitPMC(2))
        self.assertEqual(len(id2), 6)
        self.assertTrue(id2.testDelta())

    def testMorToD(self):
        id2 = identityDD(splitPMC(2))
        d = infTypeD(2)
        d2 = id2.morToD(d)
        self.assertTrue(d2.testDelta())
        d2.simplify()
        self.assertEqual(len(d2), 1)
        d3 = zeroTypeD(2)
        d4 = id2.morToD(d3)
        self.assertTrue(d4.testDelta())
        d4.simplify()
        self.assertEqual(len(d4), 1)

    def testToDStructure(self):
        ddstr = identityDD(linearPMC(1))
        dstr = ddstr.toDStructure()
        hochchild = dstr.morToD(dstr)
        self.assertEqual(len(hochchild), 18)
        hochchild.simplify(find_homology_basis = True)
        self.assertEqual(len(hochchild), 4)
        meaning_len = [len(gen.prev_meaning)
                       for gen in hochchild.getGenerators()]
        self.assertEqual(sorted(meaning_len), [1,1,1,2])

    def testIdentityMatchDiagram(self):
        pmc_to_test = [splitPMC(2), PMC([(0,2),(1,6),(3,5),(4,7)])]
        for pmc in pmc_to_test:
            ddstr = identityDD(pmc)
        # Special check for the identity diagram: all gradings should be zero
        if DEFAULT_GRADING == SMALL_GRADING:
            for gen in ddstr.generators:
                self.assertEqual(ddstr.grading[gen], ddstr.gr_set.zero())

    def testDual(self):
        pmc = PMC([(0,2),(1,6),(3,5),(4,7)])
        ddstr = identityDD(pmc)
        ddstr_dual = ddstr.dual()
        self.assertEqual(ddstr_dual.algebra1, pmc.getAlgebra())
        self.assertEqual(ddstr_dual.algebra2, pmc.opp().getAlgebra())

if __name__ == "__main__":
    unittest.main()
