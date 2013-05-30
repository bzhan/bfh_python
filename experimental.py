"""Try different things here by adding test cases. Tests added here are not
included in testmod.

"""

from braid import *
from dehntwist import *
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
        for slides in cases[0:1]+cases[3:4]:
            # Convert (b_1, c_1) into arcslides and then DD structures
            cur_pmc = splitPMC(2)
            slides_dd = []
            for b1, c1 in slides:
                arcslide = Arcslide(cur_pmc, b1, c1)
                cur_pmc = arcslide.end_pmc
                slides_dd.append(arcslide.getDDStructure(0))

            # Tensor each of the DD structures onto d_start
            d_mid = d_start
            for dd in slides_dd:
                d_mid = computeATensorDD(d_mid, dd)
                d_mid.simplify()
            d_mid.reindex()

            print "Case: %s" % slides
            print d_mid
            print d_mid.gr_set.simplifiedSet()
            for gen in d_mid.getGenerators():
                print gen, d_mid.grading[gen].simplifiedElt()

    def testTrefoilSurgery(self):
        """Computes HF for +1 and -1 surgery on left-handed trefoil. Currently
        DOES NOT APPEAR TO BE CORRECT.

        """
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
        # CFD(H_-1)
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

if __name__ == "__main__":
    unittest.main()
