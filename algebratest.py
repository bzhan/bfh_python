"""Unit test for algebra.py"""

from algebra import *
from pmc import *
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

class TensorStarTest(unittest.TestCase):
    def testTensorStarAlgebra(self):
        pmc = splitPMC(2)
        sd1 = pmc.sd([(0,1),(1,2)])
        sd2 = pmc.sd([4,(1,3)])
        cobarAlg = CobarAlgebra(pmc.getAlgebra())
        def formGen(*seq):
            return TensorStarGenerator(tuple(seq), cobarAlg)
        gen1 = formGen(sd1)
        gen2 = formGen(sd2)
        gen3 = formGen(sd1, sd2)
        gen4 = formGen()
        # Has sd3.diff() = sd1 and sd4*sd5 = sd2
        sd3 = pmc.sd([1,(0,2)])
        sd4 = pmc.sd([4,(1,2)])
        sd5 = pmc.sd([4,(2,3)])
        self.assertEqual(gen1.diff(), 1*formGen(sd3))
        self.assertEqual(gen2.diff(), 1*formGen(sd4, sd5))
        self.assertEqual(gen3.diff(),
                         1*formGen(sd1, sd4, sd5)+1*formGen(sd3, sd2))
        self.assertEqual(gen4.diff(), 0)

if __name__ == "__main__":
    unittest.main()
