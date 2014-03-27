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

    def testHochschildCoho(self):
        #Check HH^*(A)
        #Rank 4 is right for HH^* of A in genus 1 case
        pmc = splitPMC(1)
        ddstr = identityDD(pmc)
        cx = ddstr.hochschildCochains()
        cx.simplify()
        self.assertEqual(len(cx.generators),4)
        #Rank 1 in the extremal strands grading
        pmc = splitPMC(2)
        ddstr = identityDD(pmc,0)
        cx = ddstr.hochschildCochains()
        cx.simplify()
        self.assertEqual(len(cx.generators),1)
        #Rank should be indepdendent of PMC
        pmc1 = splitPMC(2)
        pmc2 = linearPMC(2)
        pmc3 = PMC([(0,2),(1,6),(3,5),(4,7)])
        ddstr1 = identityDD(pmc1)
        ddstr2 = identityDD(pmc2)
        ddstr3 = identityDD(pmc3)
        cx1 = ddstr1.hochschildCochains()
        cx1.simplify()
        cx2 = ddstr2.hochschildCochains()
        cx2.simplify()
        cx3 = ddstr3.hochschildCochains()
        cx3.simplify()
        self.assertEqual(len(cx1.generators),len(cx2.generators))
        self.assertEqual(len(cx1.generators),len(cx3.generators))
        
        
if __name__ == "__main__":
    unittest.main()
