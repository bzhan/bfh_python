"""Unit test for identityaa.py"""

from ..identityaa import *
from ..identityaa import _getIntervalOrdering
from ..pmc import antipodalPMC, linearPMC, splitPMC
import unittest

class HomotopyAATest(unittest.TestCase):
    def testHomotopyAA(self):
        HomotopyAA(splitPMC(1)).testHomotopy()
        HomotopyAA(splitPMC(2)).testHomotopy()

    def testGetIntervalOrdering(self):
        tests = [(splitPMC(1), [0,1,2]),
                 (splitPMC(2), [4,5,6,3,0,1,2]),
                 (antipodalPMC(2), [2,5,0,3,6,1,4]),
                 (linearPMC(2), [4,0,1,3,5,6,2])]
        for pmc, order in tests:
            self.assertEqual(_getIntervalOrdering(pmc), order)

if __name__ == "__main__":
    unittest.main()
