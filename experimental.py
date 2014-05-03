"""Try different things here by adding test cases. Tests added here are not
included in testmod.

"""

from braid import *
from dehntwist import *
from digraph import computeDATensorDD
from dstructure import SimpleDStructure, SimpleDGenerator
from dstructure import zeroTypeD
from ddstructure import SimpleDDGenerator, SimpleDDStructure
from ddstructure import DDStrFromDStr
from dastructure import DAStrFromChords
from utility import DEFAULT_GRADING, F2, SMALL_GRADING
import itertools
import unittest

class ExperimentalTest(unittest.TestCase):
    def testDehnTwist(self):
        slides = Braid(8).getArcslides(-5)
        assert len(slides) == 2
        print "Getting DD Structures"
        slides_dd = [slide.getDDStructure() for slide in slides]
        print "Tensoring"
        dehn_twist = computeDATensorDD(*slides_dd)
        print "Cleaning up and checks"
        dehn_twist.reindex()
        dehn_twist.checkGrading()
        self.assertTrue(dehn_twist.testDelta())

        twist = DehnTwist(3, 4, NEG)
        print "Getting DD from dehntwist"
        twist_dd = twist.getDDStructure()
        print "Comparing"
        print twist_dd.compareDDStructures(dehn_twist)

    def testAbsoluteGrading(self):
        assert DEFAULT_GRADING == SMALL_GRADING
        dd_abs_info = 1
        gr_info = [4]
        d1 = zeroTypeD(1, is_dual = False, abs_gr_info = gr_info)
        d1d = zeroTypeD(1, is_dual = True, abs_gr_info = gr_info)
        d2 = infTypeD(1, is_dual = False, abs_gr_info = gr_info)
        d2d = infTypeD(1, is_dual = True, abs_gr_info = gr_info)
        cases = [(d1d, [], d2),
                 (d2d, [], d1),
                 (d1d, [(1,0)], d2), # Dehn twist for d1d
                 (d2d, [(2,1)], d1), # Dehn twist for d2d
                 (d1d, [(3,2)], d1), # d1d -> d2d
                 (d2d, [(0,1)], d2), # d2d -> d1d
                 (d2d, [(1,0)], d1), # Dehn twist for d1
                 (d1d, [(2,1)], d2), # Dehn twist for d2
                 (d2d, [(3,2)], d2), # d2 -> d1
                 (d1d, [(0,1)], d1), # d1 -> d2
                 (d1d, [(3,2),(3,2)], d2), # Case 0 and 1
                 (d1d, [(2,1),(2,1)], d1), # Hopf link
                 (d1d, [(2,3)]*3, d1), # Trefoil
                 (d1d, [(2,3),(1,0),(2,3)], d1), # Trefoil 2
                 (d1d, [(2,3),(1,2)]*3, d2), # Boundary dehn twist
                 (d2d, [(2,3),(1,2)]*3, d1), # Boundary dehn twist, #2
                 (d2d, [(3,2)]*3, d1), # ?
                 (d1d, [(2,1),(1,0),(1,2),(2,3),(2,3)], d1), # Unknot
                 (d1d, [(3,2),(0,1)], d2)
                 ]
        for start, slides, end in cases[0:10]:
            slides_dd = [Arcslide(splitPMC(1), b1, c1).\
                             getDDStructure(dd_abs_info) for b1, c1 in slides]
            d_mid = start
            # print "start grading ", d_mid.grading
            for dd in slides_dd:
                # print dd.gr_set
                # for gen, gr in dd.grading.items():
                    # print gen, gr
                d_mid = computeATensorDD(d_mid, dd)
                d_mid.simplify()
                # print "mid grading ", d_mid.grading
            # print d_mid.gr_set.simplifiedSet()
            # for gen in d_mid.getGenerators():
                # print gen, d_mid.grading[gen].simplifiedElt()
            cur_cx = computeATensorD(d_mid, end)
            cur_cx.simplify()
            # print cur_cx.gr_set, cur_cx.grading
            # Alternate way of computing
            # d_mid = end
            # for dd in reversed(slides_dd):
            #     d_mid = computeDATensorD(dd, d_mid)
            # cur_cx = computeATensorD(start, d_mid)
            cur_abs_gr = cur_cx.getAbsGradingInfo()
            print [str(n) for n in cur_abs_gr]

    def testGenus2AbsoluteGrading(self):
        dd_abs_info = 0
        gr_info = [0,0]
        print gr_info
        d1 = zeroTypeD(2, is_dual = False, abs_gr_info = gr_info)
        d1d = zeroTypeD(2, is_dual = True, abs_gr_info = gr_info)
        d2 = infTypeD(2, is_dual = False, abs_gr_info = gr_info)
        d2d = infTypeD(2, is_dual = True, abs_gr_info = gr_info)
        cases = [(d1d, [], d2),
                 (d2d, [], d1),
                 (d1d, [(1,0)], d2), # Dehn twist for d1d
                 (d2d, [(2,1)], d1), # Dehn twist for d2d
                 (d1d, [(3,2),(7,6)], d1), # d1d -> d2d
                 (d2d, [(0,1),(4,5)], d2), # d2d -> d1d
                 (d1d, [(7,6),(0,1)], d1), # mixed
                 (d2d, [(1,0)], d1), # Dehn twist for d1
                 (d1d, [(2,1)], d2), # Dehn twist for d2
                 (d2d, [(3,2),(7,6)], d2), # d2 -> d1
                 (d1d, [(0,1),(4,5)], d1), # d1 -> d2
                 (d1d, [(3,2),(7,6),(3,4),(6,7),(5,6),(1,2),(3,4),(3,4),(5,6),
                        (5,6),(6,7),(4,3),(5,4),(6,5)], d2),
                 (d2d, [(3,4),(6,7),(5,6),(6,5),(4,3),(4,3),(2,1),(2,1),(1,0),
                        (4,3),(5,4),(6,5)], d2),
                 (d2d, [(3,4),(6,7),(5,6),(6,7),(5,6),(5,6),(3,4),(3,4),(1,2),
                        (4,3),(5,4),(6,5)], d2),
                 (d1d, [(3,4),(6,7),(5,6)] + \
                      [(6,5),(4,3),(4,3),(2,1),(2,1),(1,0)]*5 + \
                      [(4,3),(5,4),(6,5)], d2),
                 ]
        for start, slides, end in cases[0:10]:
            cur_pmc = splitPMC(2)
            slides_dd = []
            for b1, c1 in slides:
                arcslide = Arcslide(cur_pmc, b1, c1)
                cur_pmc = arcslide.end_pmc
                slides_dd.append(arcslide.getDDStructure(dd_abs_info))

            d_mid = start
            for dd in slides_dd:
                d_mid = computeATensorDD(d_mid, dd)
                d_mid.simplify()
            d_mid.reindex()
            # print d_mid
            # print d_mid.gr_set.simplifiedSet()
            # for gen in d_mid.getGenerators():
                # print gen, d_mid.grading[gen].simplifiedElt()
            cur_cx = computeATensorD(d_mid, end)

            # dd_mid = slides_dd[0]
            # for dd in slides_dd[1:]:
            #     dd_mid = computeDATensorDD(dd_mid, dd)
            #     dd_mid.simplify()
            # dd_mid.reindex()
            # # print dd_mid
            # print dd_mid.gr_set.simplifiedSet()
            # for gen in dd_mid.getGenerators():
            #     print gen, dd_mid.grading[gen].simplifiedElt()

            # Alternate way of computing
            # d_mid = end
            # for dd in reversed(slides_dd):
            #     d_mid = computeDATensorD(dd, d_mid)
            # cur_cx = computeATensorD(start, d_mid)
            cur_abs_gr = cur_cx.getAbsGradingInfo()
            print [str(n) for n in cur_abs_gr]

    def testBraidAbsoluteGrading(self):
        # getHF() does not yet support absolute grading. Test using the other
        # functions first
        std_cap = [6,3,2,5,4,1]
        to_test = [[std_cap, [], [2,1,4,3,6,5]],
                   [std_cap, [1,-1], [2,1,4,3,6,5]],
                   [std_cap, [2], [2,1,4,3,6,5]],
                   [std_cap, [-4], [2,1,4,3,6,5]]]
        for start_cap, braid_word, end_cap in to_test:
            br = BridgePresentation("br_test", start_cap, braid_word, end_cap)
            cx = br.getHF()
            # print cx.gr_set
            # for gen, gr in cx.grading.items():
            # print gen, gr
            abs_gr = cx.getAbsGradingInfo()
            print [str(n) for n in abs_gr]

    def testAlgSize(self):
        print len(splitPMC(3).getAlgebra().getGenerators())

    def testTypeDInvariant(self):
        d_start = infTypeD(2, is_dual = True, abs_gr_info = [2,2])
        # (b_1, c_1) for arcslides
        cases = [
            # Original
            [],
            # Twisting a handle
            [(2,1)],
            # Twisting a knob (half twist)
            [(2,3),(1,2)]*3,
            # Interchanging two knobs
            [(3,4),(6,7),(5,6),(4,5),(2,3),(5,6),(4,5),(3,4),
             (1,2),(4,5),(3,4),(2,3),(0,1),(3,4),(2,3),(1,2)],
            # Slide1
            [(4,3),(1,0),(1,2),(5,4),(6,5)],
            # Slide2
            [(0,1),(3,4),(6,7),(6,5),(2,3),(1,2),(3,2)]
            ]
        for slides in cases:
            # Convert (b_1, c_1) into arcslides and then DD structures
            cur_pmc = splitPMC(2)
            slides_dd = []
            for b1, c1 in slides:
                arcslide = Arcslide(cur_pmc, b1, c1)
                cur_pmc = arcslide.end_pmc
                slides_dd.append(arcslide.getDDStructure())

            # Tensor each of the DD structures onto d_start
            d_mid = d_start
            for dd in slides_dd:
                d_mid = computeATensorDD(d_mid, dd)
                d_mid.simplify()
            d_mid.reindex()

            print "Case: %s" % slides
            print d_mid
            # Rough check that this equals the original
            self.assertEquals(len(d_mid), 1)
            self.assertEquals(len(d_mid.getGenerators()[0].delta()), 2)

    def testTrefoilSurgery(self):
        """Computes HF for +1 and -1 surgery on left-handed trefoil. """
        # Everything is over the PMC of genus 1
        pmc = splitPMC(1)
        algebra = pmc.getAlgebra()
        # Two idempotents
        i0 = pmc.idem([0])
        i1 = pmc.idem([1])
        # Some algebra elements
        rho1 = pmc.sd([(0,1)])
        rho2 = pmc.sd([(1,2)])
        rho3 = pmc.sd([(2,3)])
        rho23 = pmc.sd([(1,3)])
        rho123 = pmc.sd([(0,3)])
        # Now CFD(H_+1)
        d_p1 = SimpleDStructure(F2, algebra)
        a = SimpleDGenerator(d_p1, i1, "a")
        b = SimpleDGenerator(d_p1, i0, "b")
        d_p1.addGenerator(a)
        d_p1.addGenerator(b)
        d_p1.addDelta(a, b, rho2, 1)
        d_p1.addDelta(b, a, rho123, 1)
        print "CFD(H_+1): ", d_p1
        # and CFD(H_-1)
        d_p2 = SimpleDStructure(F2, algebra)
        a = SimpleDGenerator(d_p2, i1, "a")
        b = SimpleDGenerator(d_p2, i0, "b")
        d_p2.addGenerator(a)
        d_p2.addGenerator(b)
        d_p2.addDelta(b, a, rho1, 1)
        d_p2.addDelta(b, a, rho3, 1)
        print "CFD(H_-1): ", d_p2
        # CFD(trefoil)
        d_trefoil = SimpleDStructure(F2, algebra)
        x = SimpleDGenerator(d_trefoil, i0, "x")
        y = SimpleDGenerator(d_trefoil, i0, "y")
        z = SimpleDGenerator(d_trefoil, i0, "z")
        k = SimpleDGenerator(d_trefoil, i1, "k")
        l = SimpleDGenerator(d_trefoil, i1, "l")
        mu1 = SimpleDGenerator(d_trefoil, i1, "mu1")
        mu2 = SimpleDGenerator(d_trefoil, i1, "mu2")
        for gen in [x, y, z, k, l, mu1, mu2]:
            d_trefoil.addGenerator(gen)
        d_trefoil.addDelta(x, k, rho1, 1)
        d_trefoil.addDelta(y, k, rho123, 1)
        d_trefoil.addDelta(mu2, x, rho2, 1)
        d_trefoil.addDelta(mu1, mu2, rho23, 1)
        d_trefoil.addDelta(z, mu1, rho123, 1)
        d_trefoil.addDelta(l, y, rho2, 1)
        d_trefoil.addDelta(z, l, rho3, 1)
        print "CFD(trefoil): ", d_trefoil
        # Compute the Mor complexes
        cx1 = d_p1.morToD(d_trefoil)
        # cx1 = computeATensorD(d_p1, d_trefoil)
        cx1.simplify()
        print "First result: ", cx1
        cx2 = d_p2.morToD(d_trefoil)
        # cx2 = computeATensorD(d_p2, d_trefoil)
        cx2.simplify()
        print "Second result: ", cx2

    def testLinkComplement(self):
        """Computes type DD structure associated to the complement of a certain
        link. Sequence of arcslides provided by Adam Levine.

        """
        twist1 = [(7,6),(7,6),(7,6),(7,6),(7,6),(7,6)]
        twist2 = [(4,3),(1,0),(2,1),(3,2),
                  (5,4),(2,1),(3,2),(4,3),
                  (6,5),(3,2),(4,3),(1,2),(0,1),(6,5)]
        start_pmc = splitPMC(2)
        twist_slides = {}
        twist_slides[1] = []
        twist_slides[2] = []

        cur_pmc = start_pmc
        for (b1, c1) in twist1:
            arcslide = Arcslide(cur_pmc, b1, c1)
            cur_pmc = arcslide.end_pmc
            twist_slides[1].append(arcslide)

        cur_pmc = start_pmc
        for (b1, c1) in twist2:
            arcslide = Arcslide(cur_pmc, b1, c1)
            cur_pmc = arcslide.end_pmc
            twist_slides[2].append(arcslide)

        twist_slides[-1] = [slide.inverse()
                            for slide in reversed(twist_slides[1])]
        twist_slides[-2] = [slide.inverse()
                            for slide in reversed(twist_slides[2])]

        # seq = [-2, -2, 1, 2]
        # seq = [2, 1, -2, -2]
        # seq = [-2]
        seq = [1]
        slides_total = []
        for twist in seq:
            slides_total.extend(twist_slides[twist])

        d_mid = infTypeD(2)
        # This shows our choice of starting type D structure is correct.
        # Doing this Dehn twist only should not change the starting type D
        # structure.
        # slides_total = [Arcslide(start_pmc, 6, 7)]
        # slides_total = [Arcslide(start_pmc, 7, 6)]
        for slide in slides_total:
            print slide
            print d_mid
            slide_dd = slide.getDDStructure(0)
            d_mid = slide_dd.morToD(d_mid)
            d_mid.simplify()
            d_mid.reindex()

        print d_mid
        dd_final = DDStrFromDStr(d_mid, 1)
        dd_final.testDelta()
        print dd_final
        dd_final.simplify()
        dd_final.reindex()
        print dd_final

    def testTwoStrandGenus1(self):
        # Just code to print out differential and multiplication for an algebra.
        gens = splitPMC(1).getAlgebra(
            idem_size = 2, mult_one = False).getGenerators()
        for g in gens:
            print "d(%s) = %s" % (g, g.diff())
        for g1, g2 in itertools.product(gens, gens):
            if g1.isIdempotent() or g2.isIdempotent():
                continue
            if g1 * g2 != 0:
                print "%s * %s = %s" % (g1, g2, g1*g2)

    def testT4nTorus(self):
        # Computation for torus links T(4,n).
        for n in range(1, 15):
            knot = BridgePresentation("T4_%d" % n, (8,7,6,5,4,3,2,1),
                                      [1,2,3]*n, (8,7,6,5,4,3,2,1))
            print knot.name, len(knot.getHFByLocalDA())

    def testDDStructureDelta(self):
        # Construct type DD structures, and test whether d^2 = 0 holds.

        # PMC on both sides are genus 1 split PMC.
        pmc = splitPMC(1)
        # Strand algebra corresponding to pmc.
        alg = pmc.getAlgebra()
        # Initialize type DD structure over field F_2, with (left-left) action
        # by the genus 1 strand algebra. Intend to make this type DD bimodule
        # for identity.
        ddstr1 = SimpleDDStructure(F2, alg, alg)
        # Initialize the list of generators to add to ddstr1.
        # The generators have "complementary" idempotents. However, since the
        # PMCs are in opposite direction on both sides, the vector specifying
        # idempotents are the same.
        idems = {"x" : ([0], [0]),
                 "y" : ([1], [1])}
        gens = {}
        for name, (idem1, idem2) in idems.items():
            gens[name] = SimpleDDGenerator(
                ddstr1, Idempotent(pmc, idem1), Idempotent(pmc, idem2), name)
            ddstr1.addGenerator(gens[name])
        # Now add delta
        ddstr1.addDelta(gens["x"], gens["y"],
                        pmc.sd([(0, 1)]), pmc.sd([(2, 3)]), 1)
        ddstr1.addDelta(gens["y"], gens["x"],
                        pmc.sd([(1, 2)]), pmc.sd([(1, 2)]), 1)
        ddstr1.addDelta(gens["x"], gens["y"],
                        pmc.sd([(2, 3)]), pmc.sd([(0, 1)]), 1)
        # This already satisfies d^2 = 0
        self.assertTrue(ddstr1.testDelta())
        # However, one more arrow to finish the bimodule
        ddstr1.addDelta(gens["x"], gens["y"],
                        pmc.sd([(0, 3)]), pmc.sd([(0, 3)]), 1)
        # This is now the identity bimodule, of course satisfying d^2 = 0.
        self.assertTrue(ddstr1.testDelta())

        # Second example, showing failure of testDelta()
        ddstr2 = SimpleDDStructure(F2, alg, alg)
        # Add the same generators as before
        gens = {}
        for name, (idem1, idem2) in idems.items():
            gens[name] = SimpleDDGenerator(
                ddstr2, Idempotent(pmc, idem1), Idempotent(pmc, idem2), name)
            ddstr2.addGenerator(gens[name])
        # Now add delta
        ddstr2.addDelta(gens["x"], gens["y"],
                        pmc.sd([(0, 1)]), pmc.sd([(0, 1)]), 1)
        ddstr2.addDelta(gens["y"], gens["x"],
                        pmc.sd([(1, 2)]), pmc.sd([(1, 2)]), 1)
        # Prints the type DD structure. Note the code already checks that
        # idempotent matches in all added arrows (throws an error if they don't
        # match).
        print ddstr2
        # However, testDelta() fails. Prints a term in d^2(x).
        self.assertFalse(ddstr2.testDelta())

if __name__ == "__main__":
    unittest.main()
