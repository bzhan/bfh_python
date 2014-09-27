"""Unit test for hdiagram.py"""

from hdiagram import *
from hdiagram import _Point, _Segment, _OrientedSegment, _Path, _Cell, \
    _Domain, _OneChain
from arcslide import Arcslide
from cobordism import Cobordism
from cobordism import LEFT
from pmc import antipodalPMC, linearPMC, splitPMC
import unittest

class OrientedSegmentTest(unittest.TestCase):
    def testOrientedSegment(self):
        """Testing functions in OrientedSegment."""
        p1 = _Point(1)
        p2 = _Point(2)
        seg = _Segment(12, p1, p2)
        segp = seg.oseg()
        segn = segp.opp()
        self.assertEqual(segp.start, p1)
        self.assertEqual(segp.end, p2)
        self.assertEqual(segn.start, p2)
        self.assertEqual(segn.end, p1)
        self.assertEqual(segp.toOneChain(), {seg : 1})
        self.assertEqual(segn.toOneChain(), {seg : -1})

class PathTest(unittest.TestCase):
    def setUp(self):
        self.n = 10
        self.pts = [_Point(i) for i in range(self.n)]
        self.segs = [_Segment(i, self.pts[i], self.pts[(i+1)%self.n])
                     for i in range(self.n)]
        self.osegs = [seg.oseg() for seg in self.segs]
        self.osegsr = [seg.oseg().opp() for seg in self.segs]
        self.loopseg = _Segment("loop", self.pts[0], self.pts[0])
        self.oloopseg = self.loopseg.oseg()

        # These should construct valid paths
        self.path1 = _Path([], "empty")
        self.path2 = _Path([self.oloopseg], "small_loop", True)
        self.path3 = _Path(self.osegs[0:-1], "straight")
        self.path4 = _Path(self.osegs, "loop", True)
        self.path5 = _Path(self.osegs*2, "twoloops")
        self.path6 = _Path(reversed(self.osegsr), "opploop", True)

    def testPathValidity(self):
        """Testing the validity checks in Path."""
        # Some FAIL tests
        self.assertRaises(ValueError, _Path,
                          self.osegs[0:1] + self.osegs[2:])
        self.assertRaises(ValueError, _Path,
                          self.osegs[0:-1], "notloop", True)
        self.assertRaises(ValueError, _Path, self.osegsr)

    def testOpp(self):
        """Testing the opp function in path."""
        self.assertEqual(self.path1.opp("empty"), self.path1)
        self.assertNotEqual(self.path2.opp("small_loop"), self.path2)
        self.assertEqual(self.path4.opp("opploop"), self.path6)

    def testToOneChain(self):
        """Testing toOneChain function in Path."""
        self.assertEqual(self.path1.toOneChain().values(), [])
        self.assertEqual(self.path5.toOneChain().values(), [2]*self.n)
        self.assertEqual(self.path6.toOneChain().values(), [-1]*self.n)

class CellTest(unittest.TestCase):
    def testCell(self):
        """Testing functions in Cell."""
        pts = [_Point(i) for i in range(3)]
        segs = [_Segment("0", pts[0], pts[1]),
                _Segment("1", pts[1], pts[0]),
                _Segment("2", pts[2], pts[2])]
        cell1 = _Cell("c1", _Path([segs[0].oseg(), segs[1].oseg()],
                                  iscycle = True))
        cell2 = _Cell("c2", [_Path([segs[0].oseg(), segs[1].oseg()],
                                   iscycle = True),
                             _Path([segs[2].oseg()], iscycle = True)])
        self.assertEqual(cell1.toDomain(), {cell1 : 1})
        self.assertEqual(cell1.bdOneChain(), {segs[0]:1, segs[1]:1})
        self.assertEqual(cell2.bdOneChain(), {segs[0]:1, segs[1]:1, segs[2]:1})

class DomainTest(unittest.TestCase):
    def testDiff(self):
        """Testing diff function in Domain."""
        pts = [_Point(i) for i in range(4)]
        segs = [_Segment("0", pts[0], pts[1]),
                _Segment("1", pts[1], pts[2]),
                _Segment("2", pts[2], pts[0]),
                _Segment("3", pts[1], pts[3]),
                _Segment("4", pts[3], pts[0])]
        cell1 = _Cell("c1", _Path([seg.oseg() for seg in segs[0:3]],
                                      iscycle = True))
        cell2 = _Cell("c2", _Path([seg.oseg() for seg in segs[0:1]+segs[3:]],
                                     iscycle = True).opp())
        domain1 = cell1.toDomain() + cell2.toDomain()
        self.assertEqual(domain1.diff(),
                          {segs[1]:1, segs[2]:1, segs[3]:-1, segs[4]:-1})

class DiagramBuildTest(unittest.TestCase):
    def testDiagramFromCycleInfo(self):
        """Testing the function diagramFromCycleInfo."""
        # A standard diagram for solid torus
        diagram1 = diagramFromCycleInfo("Solid torus",
            num_interior_point = 1, length_border = [4],
            alpha_arcs = [[(0,0),0,(0,2)], [(0,1),(0,3)]],
            beta_cycles = [[0]], crossing_orientation = [-1])

        # A standard diagram for the identity cobordism
        diagram2 = diagramFromCycleInfo("Identity cobordism",
            num_interior_point = 4, length_border = [4,4],
            alpha_arcs = [[(0,0),3,(0,2)], [(0,1),1,(0,3)],
                          [(1,0),0,(1,2)], [(1,1),2,(1,3)]],
            beta_cycles = [[0,1],[2,3]], crossing_orientation = [1,-1,1,-1])

        # A standard diagram for the anti-braid resolution
        diagram3 = diagramFromCycleInfo("Antibraid resolution",
            num_interior_point = 3, length_border = [4,4],
            alpha_arcs = [[(0,0),2,(0,2)], [(0,1),(0,3)],
                          [(1,0),(1,2)], [(1,1),0,1,(1,3)]],
            beta_cycles = [[0], [1,2]], crossing_orientation = [1,1,-1])

        # Uncomment to see full printout of diagrams
        # print repr(diagram1), repr(diagram2), repr(diagram3)

class CommonDiagramsTest(unittest.TestCase):
    def testIdentityDiagram(self):
        pmc_to_test = [splitPMC(1), splitPMC(2), linearPMC(2), antipodalPMC(2)]
        pmc_to_test.append(PMC([(0,2),(1,6),(3,5),(4,7)]))
        pmc_to_test += [splitPMC(3), antipodalPMC(4)]
        genus_to_size = [2, 6, 20, 70]
        for pmc in pmc_to_test:
            diagram = getIdentityDiagram(pmc)
            self.assertEqual(diagram.getPMCs(), [pmc.opp(), pmc])
            expected_size = genus_to_size[pmc.genus-1]
            self.assertEqual(len(diagram.getHFGenerators()), expected_size)
            self.assertEqual(len(diagram.getPeriodicDomains()), pmc.genus*2)

    def testArcslideDiagram(self):
        slide_to_test = [Arcslide(splitPMC(2), 4, 3),
                         Arcslide(splitPMC(2), 2, 3)]
        for slide in slide_to_test:
            diagram = getArcslideDiagram(slide)
            self.assertEqual(diagram.getPMCs(),
                              [slide.start_pmc.opp(), slide.end_pmc])
            self.assertEqual(len(diagram.getHFGenerators()), 8)
            periodic_domains = diagram.getPeriodicDomains()
            self.assertEqual(len(periodic_domains), 4)
            for domain in periodic_domains:
                alpha_bd = diagram.restrictOneChain(domain.diff(), ALPHA)
                self.assertEqual(diagram.restrictZeroChain(alpha_bd.diff()), 0)

    def testHandlebodyDiagram(self):
        diagram = getInfFrameDiagram(2)
        self.assertEqual(diagram.getPMCs(), [splitPMC(2)])
        self.assertEqual(len(diagram.getHFGenerators()), 1)
        periodic_domains = diagram.getPeriodicDomains()
        self.assertEqual(len(periodic_domains), 2)
        for domain in periodic_domains:
            alpha_bd = diagram.restrictOneChain(domain.diff(), ALPHA)
            self.assertEqual(diagram.restrictZeroChain(alpha_bd.diff()), 0)
        diagram2 = getPlatDiagram(2)
        # Uncomment to see full printout of diagrams
        # print repr(diagram), repr(diagram2)

    def testAdmHandlebodyDiagram(self):
        diagram = getZeroFrameAdmDiagram(2)
        # Uncomment to see full printout of diagrams
        # print repr(diagram)

    def testSimpleCobordismDiagram(self):
        diagram = getSimpleCobordismDiagram(splitPMC(1), 1)
        # Uncomment to see full printout of diagrams
        # print repr(diagram)

    def testGetCobordismDiagramLeft(self):
        for c_pair in [0, 2, 3]:
            diagram = getCobordismDiagramLeft(Cobordism(2, c_pair, LEFT))
            # Uncomment to see full printout of diagrams
            # print repr(diagram)

    def testConnectingDomain(self):
        diagram = getIdentityDiagram(splitPMC(2))
        gens = diagram.getHFGenerators()
        for x in gens:
            for y in gens:
                domain = diagram.getConnectingDomain(x, y)
                self.assertTrue(domain is not None)
                # Alpha curves go from x to y
                alpha_bd = diagram.restrictOneChain(domain.diff(), ALPHA)
                self.assertEqual(diagram.restrictZeroChain(alpha_bd.diff()),
                                 y.toZeroChain() - x.toZeroChain())
                self.assertEqual(diagram.getMaslov(domain, x, y), 0)
                gr = diagram.getBigGrading(domain, x, y)
                self.assertEqual(gr[0] * gr[1].Ropp(), 0)

class ComputeGradingTest(unittest.TestCase):
    def testDDBigGrading(self):
        # Just check it will compute grading without raising errors.
        # Correctness of grading is crossed checked with DD structures.
        diagram = getIdentityDiagram(splitPMC(2))
        base_gen = diagram.getHFGenerators()[0]
        gr_set, gr_vals = diagram.computeDDGrading(base_gen)

if __name__ == "__main__":
    unittest.main()
