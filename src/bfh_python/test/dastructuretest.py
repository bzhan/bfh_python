"""Unit test for dastructure.py"""

from ..dastructure import *
from ..arcslide import Arcslide
from ..arcslideda import ArcslideDA
from ..dstructure import infTypeD, zeroTypeD
from ..pmc import splitPMC
from ..utility import DEFAULT_GRADING, SMALL_GRADING
import unittest

class DAStructureTest(unittest.TestCase):
    def testIdentityDA(self):
        dastr = identityDA(splitPMC(2))
        self.assertTrue(dastr.testDelta())

    def testIdentityDAMatchDiagram(self):
        dastr = identityDA(splitPMC(2))
        dastr.checkGrading()
        if DEFAULT_GRADING == SMALL_GRADING:
            # Special check for the identity diagram: all gradings should be
            # zero
            for gen in dastr.getGenerators():
                self.assertEqual(dastr.grading[gen], dastr.gr_set.zero())

class TensorTest(unittest.TestCase):
    def testDATensorD(self):
        dastr = identityDA(splitPMC(2))
        dstr = zeroTypeD(2)
        dstr_result = dastr.tensorD(dstr)
        cx = dstr_result.morToD(infTypeD(2))
        cx.simplify()
        # Basic check that dstr_result is still zeroTypeD(2)
        self.assertEqual(len(cx), 1)

if __name__ == "__main__":
    unittest.main()
