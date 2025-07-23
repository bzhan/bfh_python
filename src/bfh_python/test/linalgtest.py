"""Unit test for linalg.py."""

from ..linalg import *
from ..utility import F2
import unittest

class RowSystemTest(unittest.TestCase):
    def testRowSystem(self):
        rows1 = [[2,0],[3,1]]
        sys1 = RowSystem(rows1)
        self.assertEqual(sys1.getComb([1, 1]), [-1, 1])
        self.assertEqual(sys1.getComb([0,-2]), [ 3,-2])
        self.assertEqual(sys1.getComb([1,0]), None)
        self.assertEqual(sys1.vecReduce([1,0]), ([2,-1], [0,1]))
        self.assertEqual(sys1.vecReduce([1,0], use_rational = True),
                         ([Fraction(1,2),0], [0,0]))
        self.assertEqual(sys1.reduceProfile(), [1, 2])
        self.assertEqual(sys1.shortForm([1,0]), [1])

        rows2 = [[2,0,0],[0,2,2],[1,1,1]]
        sys2 = RowSystem(rows2)
        self.assertTrue(sys2.getZeroComb() in ([[1,1,-2]],[[-1,-1,2]]))
        self.assertEqual(sys2.vecReduce([1,0,0])[1], [0,1,1]) # reduced form
        self.assertEqual(sys2.reduceProfile(), [1, 2, 0])

class F2RowSystemTest(unittest.TestCase):
    def testF2RowSystem(self):
        rows1 = [[1,1,0,0,0],
                 [0,0,1,1,0],
                 [0,0,0,0,1],
                 [0,0,1,0,0]]
        sys1 = F2RowSystem(rows1)
        self.assertEqual(sys1.getComb([1,1,1,1,1]), [1,1,1,0])
        self.assertEqual(sys1.getComb([1,1,0,0,1]), [1,0,1,0])
        self.assertEqual(sys1.getComb([1,1,1,0,1]), [1,0,1,1])
        self.assertEqual(sys1.getComb([1,1,0,1,1]), [1,1,1,1])

if __name__ == "__main__":
    unittest.main()
