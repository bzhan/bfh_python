"""Unit test for minusalg.py"""

from minusalg import *
from ddstructure import SimpleDDGenerator, SimpleDDStructure
from pmc import Idempotent
import unittest

class MinusAlgTest(unittest.TestCase):
    def testGenerators(self):
        gens = MinusStrandAlgebra(F2, splitPMC(1)).getGenerators()
        self.assertEqual(len(gens), 18)

#    def testMultiply(self):
#        gens = MinusStrandAlgebra(F2, splitPMC(1)).getGenerators()
#        for gen in gens:
#            # Number of factors should be 1 + (length of strand)
#            self.assertEqual(len(gen.factor()), 1 + sum(gen.multiplicity))

    def testHochchild(self):
        pmc = splitPMC(1)
        alg = MinusStrandAlgebra(F2, pmc)
        ddstr = SimpleDDStructure(F2, alg, alg)
        # Initialize the list of generators to add to ddstr1.
        idems = {"x" : ([0], [0]),
                 "y" : ([1], [1])}
        gens = {}
        for name, (idem1, idem2) in list(idems.items()):
            gens[name] = SimpleDDGenerator(
                ddstr, Idempotent(pmc, idem1), Idempotent(pmc, idem2), name)
            ddstr.addGenerator(gens[name])
        # Now add delta
        ddstr.addDelta(gens["x"], gens["y"],
                       minusSD(pmc, [(0, 1)]), minusSD(pmc, [(2, 3)]), 1)
        ddstr.addDelta(gens["y"], gens["x"],
                       minusSD(pmc, [(1, 2)]), minusSD(pmc, [(1, 2)]), 1)
        ddstr.addDelta(gens["x"], gens["y"],
                       minusSD(pmc, [(2, 3)]), minusSD(pmc, [(0, 1)]), 1)
        ddstr.addDelta(gens["y"], gens["x"],
                       minusSD(pmc, [(3, 0)]), minusSD(pmc, [(3, 0)]), 1)
        print(ddstr)
        self.assertTrue(ddstr.testDelta())

        dstr = ddstr.toDStructure()
        print(dstr)
        self.assertTrue(dstr.testDelta())
        hochchild = dstr.morToD(dstr)
        print(hochchild)
        hochchild.simplify(find_homology_basis = True)
        print(len(hochchild))
        meaning_len = [len(gen.prev_meaning)
                       for gen in hochchild.getGenerators()]
        for gen in hochchild.getGenerators():
            print(gen.prev_meaning)

    def testLargeChainComplex(self):
        getHalfIdComplex()
        
if __name__ == "__main__":
    unittest.main()
