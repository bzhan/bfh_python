"""Unit test for algebra.py"""

from algebra import *
from pmc import splitPMC
from utility import ZZ
import unittest

class ChainComplexTest(unittest.TestCase):
    def testChainComplex(self):
        cx = SimpleChainComplex(F2)
        gens = [SimpleGenerator(cx, "gen%d"%i) for i in range(3)]
        for gen in gens:
            cx.addGenerator(gen)
        cx.addDifferential(gens[1], gens[0], 1)
        cx.addDifferential(gens[2], gens[0], 1)

        self.assertEqual(gens[0].diff(), 0)
        self.assertEqual(gens[1].diff(), 1*gens[0])
        self.assertEqual(gens[2].diff(), 1*gens[0])
        elt12 = 1*gens[1] + 1*gens[2]
        self.assertEqual(elt12.diff(), 0)
        self.assertEqual(len(cx), 3)
        cx.reindex()
        self.assertEqual(len(cx), 3)

    def testChainComplexOverZ(self):
        cx = SimpleChainComplex(ZZ)
        gens = [SimpleGenerator(cx, "gen%d"%i) for i in range(4)]
        for gen in gens:
            cx.addGenerator(gen)
        cx.addDifferential(gens[1], gens[0], 1)
        cx.addDifferential(gens[2], gens[0], 1)
        cx.addDifferential(gens[1], gens[3], 1)
        cx.addDifferential(gens[2], gens[3], -1)

        self.assertEqual(gens[0].diff(), 0)
        self.assertEqual(gens[1].diff(), 1*gens[0]+1*gens[3])
        self.assertEqual(gens[2].diff(), 1*gens[0]-1*gens[3])
        elt12 = 1*gens[1] + 1*gens[2]
        self.assertEqual(elt12.diff(), 2*gens[0])

class TensorStarTest(unittest.TestCase):
    def testTensorStarAlgebra(self):
        pmc = splitPMC(2)
        sd1 = pmc.sd([(0,1),(1,2)])
        sd2 = pmc.sd([4,(1,3)])
        sd3 = pmc.sd([4,(0,1)])
        cobarAlg = CobarAlgebra(pmc.getAlgebra())
        def formGen(*seq):
            return TensorStarGenerator(tuple(seq), cobarAlg)
        gen1 = formGen(sd1)
        gen2 = formGen(sd2)
        gen3 = formGen(sd1, sd2)
        gen4 = TensorStarGenerator(tuple(), cobarAlg, pmc.getIdempotents()[0])
        gen5 = formGen(sd3, sd3, sd3)
        # Has sd3.diff() = sd1 and sd4*sd5 = sd2
        sd4 = pmc.sd([1,(0,2)])
        sd5 = pmc.sd([4,(1,2)])
        sd6 = pmc.sd([4,(2,3)])
        self.assertEqual(gen1.diff(), 1*formGen(sd4))
        self.assertEqual(gen2.diff(), 1*formGen(sd5, sd6))
        self.assertEqual(gen3.diff(),
                         1*formGen(sd1, sd5, sd6)+1*formGen(sd4, sd2))
        self.assertEqual(gen4.diff(), 0)
        self.assertEqual(gen5.diff(), 0)

class SolveF2SystemTest(unittest.TestCase):
    def testSolveF2System(self):
        # Represents the matrix:
        # [[1,1,0,0,0],
        #  [0,0,1,1,0],
        #  [0,0,0,0,1],
        #  [0,0,1,0,0]]
        num_row, num_col = 4, 5
        entries = [(0, 0), (0, 1), (1, 2), (1, 3), (2, 4), (3, 2)]
        for vec, solution in [
                ([], []),
                ([0,1,2,3,4], [0,1,2]),
                ([0,1,4], [0,2]),
                ([0,1,2,4], [0,2,3]),
                ([0,1,3,4], [0,1,2,3]),
                ([0], None),]:
            self.assertEqual(solveOverF2(num_row, num_col, entries, vec),
                             solution)

if __name__ == "__main__":
    unittest.main()
