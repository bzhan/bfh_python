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

if __name__ == "__main__":
    unittest.main()
