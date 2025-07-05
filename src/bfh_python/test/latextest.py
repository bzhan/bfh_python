"""Unit test for latex.py"""

from ..latex import *
from ..arcslide import Arcslide
from ..arcslideda import ArcslideDA
from ..localpmc import LocalPMC, LocalStrandDiagram
from ..pmc import PMC
from ..pmc import antipodalPMC, splitPMC
import unittest
import io

class LatexTest(unittest.TestCase):
    def testPrintLocalDAStructure(self):
        slide, pmc_map1, pmc_map2 = (
            # Uncomment one of the following lines
            # Arcslide(splitPMC(2), 2, 1), None, None,  # short, down
            # Arcslide(splitPMC(2), 2, 3), None, None,  # short, up
            # Arcslide(antipodalPMC(2), 2, 1), [0,1,2,3,4,6,7], [0,1,3,4,5,6,7],
            Arcslide(antipodalPMC(2), 4, 5), [0,1,3,4,5,6,7], [0,1,2,3,4,6,7],
        )
        local_dastr = ArcslideDA(slide).getLocalDAStructure()
        f = io.StringIO()
        f.write(beginDoc())
        f.write(showDAStructure(local_dastr, pmc_map1, pmc_map2))
        f.write(endDoc())
        f.close()

if __name__ == "__main__":
    unittest.main()
