"""Unit test for dastructure.py"""

from dastructure import *
import unittest

class DAStructureTest(unittest.TestCase):
    def testIdentityDA(self):
        dastr = identityDA(splitPMC(2))
        self.assertTrue(dastr.testDelta())

if __name__ == "__main__":
    unittest.main()
