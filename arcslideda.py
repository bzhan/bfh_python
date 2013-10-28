"""Producing type DA structures for arcslides, using local actions."""

from algebra import E0
from dastructure import SimpleDAGenerator, SimpleDAStructure
from dastructure import AddChordToDA
from linalg import F2RowSystem
from localpmc import LocalIdempotent, LocalStrandAlgebra, LocalStrandDiagram
from localpmc import restrictPMC, restrictStrandDiagram
from pmc import Strands, StrandDiagram
from utility import subset
from utility import F2
import itertools
from Queue import Queue

class ArcslideDA:
    """Responsible for producing a type DA structure for an arcslide, using
    local actions.

    """
    def __init__(self, slide):
        """Specifies the arcslide to use. slide should be of type Arcslide.
        In addition to recording slide, construct the following:
        local_pmc1, mapping1 - restriction of starting pmc to location of slide.
        outer_pmc1, outer_mapping1 - complement of slide in starting pmc.
        local_pmc2, mapping2, outer_pmc2, outer_mapping2
          - same, but for ending pmc.
        arrow_patterns - key is tuple of LocalStrandDiagram specifying what the
          A-side inputs look like in the local PMC. Value is a list of possible
          local D-side outputs.

        """
        self.slide = slide

        pmc1, pmc2 = slide.start_pmc, slide.end_pmc
        n = pmc1.n
        b1, c1, c2 = slide.b1, slide.c1, slide.c2

        if b1 == 1 and c1 == 0 and c2 == 2:  # c1-b1-c2-*
            local_cut1, outer_cut1, local_cut2, outer_cut2 = (
                [(0, 2)], [(3, n-1)],
                [(0, 2)], [(3, n-1)])
            patterns_raw = ArcslideDA._short_underslide_down_bottom()
        elif b1 == n-2 and c1 == n-3 and c2 == n-1:  # *-c1-b1-c2
            local_cut1, outer_cut1, local_cut2, outer_cut2 = (
                [(n-3, n-1)], [(0, n-4)],
                [(n-3, n-1)], [(0, n-4)])
            patterns_raw = ArcslideDA._short_underslide_down_top()
        elif b1 == c1 + 1 and c2 == c1 + 2:  # *-c1-b1-c2-*
            assert c1 != 0 and c2 != n-1
            local_cut1, outer_cut1, local_cut2, outer_cut2 = (
                [(c1, c2)], [(0, c1-1), (c2+1, n-1)],
                [(c1, c2)], [(0, c1-1), (c2+1, n-1)])
            patterns_raw = ArcslideDA._short_underslide_down_middle()
        elif b1 == 1 and c1 == 0 and c2 == 3:  # c1-b1-x-c2-*
            local_cut1, outer_cut1, local_cut2, outer_cut2 = (
                [(0, 3)], [(4, n-1)],
                [(0, 3)], [(4, n-1)])
            patterns_raw = ArcslideDA._length3_underslide_down_bottom()
        # elif b1 == c1 + 1 and c2 == c1 + 3:  # *-c1-b1-x-c2-*
        #     assert c1 != 0 and c2 != n-1
        #     local_cut1, outer_cut1, local_cut2, outer_cut2 = (
        #         [(c1, c2)], [(0, c1-1), (c2+1, n-1)],
        #         [(c1, c2)], [(0, c1-1), (c2+1, n-1)])
        #     patterns_raw = ArcslideDA._length3_underslide_down_middle()
        elif b1 == c1 + 1 and c2 > b1:  # *-c1-b1-*, *-c2-*
            assert c1 > 0 and c2 != n-1 and c2 - c1 >= 3
            b1p = c2 - 1
            local_cut1, outer_cut1, local_cut2, outer_cut2 = (
                [(c1, b1), (c2, c2)], [(0, c1-1), (b1+1, c2-1), (c2+1, n-1)],
                [(c1, c1), (b1p, c2)], [(0, c1-1), (c1+1, b1p-1), (c2+1, n-1)])
            patterns_raw = ArcslideDA._general_underslide_down_middle()
        else:
            raise NotImplementedError(
                "This slide pattern is not yet implemented.")

        self.local_pmc1, self.mapping1 = restrictPMC(pmc1, local_cut1)
        self.outer_pmc1, self.outer_mapping1 = restrictPMC(pmc1, outer_cut1)
        self.local_pmc2, self.mapping2 = restrictPMC(pmc2, local_cut2)
        self.outer_pmc2, self.outer_mapping2 = restrictPMC(pmc2, outer_cut2)

        # Required so the left to right transition on the outside can proceed.
        assert self.outer_pmc1 == self.outer_pmc2

        self.arrow_patterns = {}
        for pattern in patterns_raw:
            key = []
            for i in range(len(pattern)-1):
                key.append(self.local_pmc2.sd(pattern[i]))
            key = tuple(key)
            if key not in self.arrow_patterns:
                self.arrow_patterns[key] = []
            self.arrow_patterns[key].append(self.local_pmc1.sd(pattern[-1]))

    @staticmethod
    def idemMatchDA(x, y, coeff_d, coeffs_a):
        """Tests whether idempotent matches in the potential arrow
        x * coeffs_a -> coeff_d * y.

        """
        if x.idem1 != coeff_d.left_idem or y.idem1 != coeff_d.right_idem:
            return False
        if len(coeffs_a) == 0:
            return x.idem2 == y.idem2
        else:
            return x.idem2 == coeffs_a[0].left_idem and \
                y.idem2 == coeffs_a[-1].right_idem

    def getDAStructure(self):
        """Returns the type DA structure corresponding to slide."""
        dd_idems = self.slide.getIdems()
        da_idems = [(l_idem, r_idem.opp().comp())
                    for l_idem, r_idem in dd_idems]
        pmc1, pmc2 = self.slide.start_pmc, self.slide.end_pmc
        n = pmc1.n
        alg1, alg2 = pmc1.getAlgebra(), pmc2.getAlgebra()
        dastr = SimpleDAStructure(F2, alg1, alg2)
        for i in range(len(da_idems)):
            l_idem, r_idem = da_idems[i]
            dastr.addGenerator(
                SimpleDAGenerator(dastr, l_idem, r_idem, "%d" % i))

        alg1_gens, alg2_gens = alg1.getGenerators(), alg2.getGenerators()
        mod_gens = dastr.getGenerators()

        # Add action with zero algebra input
        short_chord = [tuple(sorted([self.slide.b1, self.slide.c1]))]
        AddChordToDA(dastr, Strands(pmc1, short_chord), [])

        def softMatch(coeffs_a, coeffs_a_ref, prefix = False):
            """Returns True if and only if coeffs_a and coeffs_a_ref differ only
            in single idempotents, with coeffs_a having more single idempotents,
            and at least one coeffs_a agrees with the corresponding
            coeffs_a_ref.

            If prefix is True, then no need for at least one coeffs_a to agree
            with the corresponding coeffs_a_ref (used for prefix matching).

            """
            if len(coeffs_a) != len(coeffs_a_ref):
                return False
            for i in range(len(coeffs_a)):
                if (coeffs_a[i].strands != coeffs_a_ref[i].strands or \
                    coeffs_a[i].double_hor != coeffs_a_ref[i].double_hor or \
                    not all([idem in coeffs_a[i].single_hor for idem in
                             coeffs_a_ref[i].single_hor])):
                    return False

            return prefix or any([coeffs_a[i] == coeffs_a_ref[i]
                                  for i in range(len(coeffs_a))])

        def adjustSingleHors(coeffs_a):
            """Given a tuple of A-side inputs, add or delete inputs in a best
            effort to make the idempotents match.

            """
            coeffs_a = list(coeffs_a)
            # Two cases: delete single-hor in coeffs_a[i] (preferred) or add
            # single-hor to coeffs_a[i+1].
            for i in range(len(coeffs_a)-1):
                if coeffs_a[i].right_idem == coeffs_a[i+1].left_idem:
                    continue
                for idem in coeffs_a[i].right_idem:
                    if idem not in coeffs_a[i+1].left_idem and \
                       idem in coeffs_a[i].single_hor:
                        coeffs_a[i] = LocalStrandDiagram(
                            coeffs_a[i].parent,
                            [j for j in coeffs_a[i].left_idem if j != idem],
                            coeffs_a[i].strands)
                for idem in coeffs_a[i+1].left_idem:
                    if idem not in coeffs_a[i].right_idem and \
                       idem in coeffs_a[i+1].single_hor:
                        coeffs_a[i+1] = LocalStrandDiagram(
                            coeffs_a[i+1].parent,
                            [j for j in coeffs_a[i+1].left_idem if j != idem],
                            coeffs_a[i+1].strands)
            return tuple(coeffs_a)

        def search(cur_list, cur_list_local, cur_prod_d):
            """Find arrows matching one of the local actions by recursively
            searching through possible lists of algebra inputs. The parameters
            are:
            - cur_list: current list of algebra inputs.
            - cur_list_local: restrictions of current list of algebra inputs to
              the local PMC. Must match the prefix of one of the local patterns.
            - cur_prod_d: product of the restrictions of algebra generators to
              the outside local PMC. So named since it equals the restriction of
              the D-side output to the outside local PMC. Must not be None,
              except at the beginning (when cur_list and cur_list_local are
              empty).

            """
            for cur_a in alg2_gens:
                if cur_a.isIdempotent():
                    continue
                if len(cur_list) > 0 and \
                   cur_list[-1].right_idem != cur_a.left_idem:
                    continue
                new_list = cur_list + (cur_a,)
                new_list_local = cur_list_local + (restrictStrandDiagram(
                    pmc2, cur_a, self.local_pmc2, self.mapping2),)
                outer_a = restrictStrandDiagram(
                    pmc2, cur_a, self.outer_pmc2, self.outer_mapping2)
                # Compute product on the outside
                if cur_prod_d is None:
                    new_prod_d = 1 * outer_a
                else:
                    new_prod_d = cur_prod_d.parent.multiplyGeneral(
                        cur_prod_d, outer_a, False)  # strict_idems = False
                if new_prod_d == 0:
                    continue
                new_prod_d = new_prod_d.getElt()
                # Local patterns match exactly one of the arrows
                for pattern in self.arrow_patterns:
                    if not softMatch(adjustSingleHors(new_list_local), pattern):
                        continue
                    for local_d in self.arrow_patterns[pattern]:
                        alg_d = local_d.join(new_prod_d.removeSingleHor(), pmc1,
                                             self.mapping1, self.outer_mapping1)
                        if alg_d is None:
                            continue
                        for x, y in itertools.product(mod_gens, mod_gens):
                            if ArcslideDA.idemMatchDA(x, y, alg_d, new_list):
                                dastr.addDelta(x, y, alg_d, new_list, 1)
                # Local patterns match the prefix of one of the arrows
                if any([softMatch(new_list_local, pattern[0:len(new_list)],
                                  prefix = True)
                        for pattern in self.arrow_patterns
                        if len(pattern) > len(new_list)]):
                    search(new_list, new_list_local, new_prod_d)

        search((), (), None)
        return dastr

    def getLocalDAStructure(self):
        """Returns the local type DA structure associated to slide. Mainly for
        testing (that it satisfies the type DA structure equations.

        """
        alg1 = LocalStrandAlgebra(F2, self.local_pmc1)
        alg2 = LocalStrandAlgebra(F2, self.local_pmc2)
        local_dastr = SimpleDAStructure(F2, alg1, alg2)

        # Mappings between local starting and ending PMC.
        slide = self.slide
        local_to_r = dict()
        for i in range(slide.start_pmc.n):
            if i in self.mapping1:
                # to_r[i] must be in mapping2
                local_to_r[self.mapping1[i]] = self.mapping2[slide.to_r[i]]
        local_pair_to_r = dict()
        for i in range(self.local_pmc1.n):
            if i not in self.local_pmc1.endpoints:
                local_pair_to_r[self.local_pmc1.pairid[i]] \
                    = self.local_pmc2.pairid[local_to_r[i]]
            
        b1, c1 = self.slide.b1, self.slide.c1
        local_b1, local_c1 = self.mapping1[b1], self.mapping1[c1]
        b_pair1 = self.local_pmc1.pairid[local_b1]
        c_pair1 = self.local_pmc1.pairid[local_c1]

        da_idems = []
        num_pair = self.local_pmc1.num_pair
        for idem in subset(range(num_pair)):
            da_idems.append((list(idem), [local_pair_to_r[p] for p in idem]))
        for idem in subset([p for p in range(num_pair)
                            if p != b_pair1 and p != c_pair1]):
            da_idems.append((list(idem) + [c_pair1],
                             [local_pair_to_r[p]
                              for p in (list(idem) + [b_pair1])]))
        for i in range(len(da_idems)):
            l_idem, r_idem = da_idems[i]
            local_dastr.addGenerator(SimpleDAGenerator(
                local_dastr, LocalIdempotent(self.local_pmc1, l_idem),
                LocalIdempotent(self.local_pmc2, r_idem), "%d" % i))
        mod_gens = local_dastr.getGenerators()
        for coeffs_a in self.arrow_patterns.keys():
            if len(coeffs_a) == 0 or coeffs_a[0].isIdempotent():
                continue
            for coeff_d in self.arrow_patterns[coeffs_a]:
                arrow_used = False
                for x in mod_gens:
                    for y in mod_gens:
                        if x.idem1 == coeff_d.left_idem and \
                           x.idem2 == coeffs_a[0].left_idem and \
                           y.idem1 == coeff_d.right_idem and \
                           y.idem2 == coeffs_a[-1].right_idem:
                            local_dastr.addDelta(x, y, coeff_d, coeffs_a, 1)
                            arrow_used = True
                if not arrow_used:
                    print "Warning: unused arrow: ", coeffs_a, coeff_d
                    pass
        # Now add the arrows with no A-side inputs.
        local_short_chord = [tuple(sorted([local_b1, local_c1]))]
        empty_idems = [i for i in range(self.local_pmc1.num_pair)
                       if i != b_pair1 and i != c_pair1]
        for idems_to_add in subset(empty_idems):
            alg_d = LocalStrandDiagram(alg1, [c_pair1] + list(idems_to_add),
                                       local_short_chord)
            for x in mod_gens:
                for y in mod_gens:
                    if x.idem2 == y.idem2 and \
                       x.idem1 == alg_d.left_idem and \
                       y.idem1 == alg_d.right_idem:
                        local_dastr.addDelta(x, y, alg_d, [], 1)
        return local_dastr

    @staticmethod
    def _restrict_local_PMC():
        pass

    # The next series of functions specify the local arrows. The format is as
    # follows:
    # All but the last element of the tuple is a list to be passed to the sd()
    # function of the A-side local PMC, specifying the A-side inputs. The last
    # element of the tuple is a list to be passed to the sd() function of the
    # D-side local PMC, specifying the D-side output.
    @staticmethod
    def _short_underslide_down_bottom():
        """Short underslide going up, on the bottom of PMC."""
        # Local PMC is 0-1-2-3*, with 0 and 2 paired.
        patterns_raw = [
            #### Single patterns
            # () -> ()
            ([], []),  ([0], [0]), ([1], [1]), ([1], [0]),
            # (1, 2) -> ()
            ([(1, 2)], [0]),
            # (0, 1) -> (0, 2)
            ([(0, 1)], [(0, 2)]),
            # (0, 2) -> (0, 2)
            ([(0, 2)], [(0, 2)]),
            ([1, (0, 2)], [1, (0, 2)]),
            # (0, 1)-(1, 2) -> (0, 1)-(1, 2)
            ([(0, 1),(1, 2)], [(0, 1),(1, 2)]),
            # (2, 3) -> (2, 3) extension of null
            ([(2, 3)], [(2, 3)]),
            ([1, (2, 3)], [1, (2, 3)]),
            # (1, 3) -> (2, 3) extension of (1, 2) -> ()
            ([(1, 3)], [(2, 3)]),
            # (0, 3) -> (0, 3) extension of (0, 2) -> (0, 2)
            ([(0, 3)], [(0, 3)]),
            ([1, (0, 3)], [1, (0, 3)]),
            # (0, 1)-(1, 3) -> (0, 1)-(1, 3) extension of (0, 1)-(1, 2) -> ...
            ([(0, 1),(1, 3)], [(0, 1),(1, 3)]),

            # Double patterns
            # (1, 2), (0, 1) -> (1, 2)
            ([(1, 2)], [(0, 1)], [(1, 2)]),
            # (1, 2), (0, 2) -> (1, 2)
            ([(1, 2)], [(0, 2)], [(1, 2)]),
            # *** Extensions of (1, 2), (0, 1) -> (1, 2) ***
            # (1, 2)-(2, 3), (0, 1) -> (1, 2)-(2, 3)
            ([(1, 2),(2, 3)], [(0, 1)], [(1, 2),(2, 3)]),
            # (1, 3), (0, 1) -> (1, 3)
            ([0, (1, 3)], [(0, 1)], [0, (1, 3)]),
            # *** Extensions of (1, 2), (0, 2) -> (1, 2) ***
            # (1, 2), (0, 3) -> (1, 3)
            ([(1, 2)], [(0, 3)], [(1, 3)]),
            # (1, 2)-(2, 3), (0, 2) -> (1, 2)-(2, 3)
            ([(1, 2),(2, 3)], [(0, 2)], [(1, 2),(2, 3)]),
            # (1, 3), (0, 2) -> (1, 3)
            ([0, (1, 3)], [(0, 2)], [0, (1, 3)]),
        ]
        return patterns_raw

    @staticmethod
    def _short_underslide_down_top():
        """Short underslide going up, on the top of PMC."""
        # Local PMC is 0*-1-2-3, with 1 and 2 paired.
        patterns_raw = [
            #### Single patterns
            # () -> ()
            ([], []), ([1], [1]), ([2], [2]), ([2], [1]),
            # (2, 3) -> ()
            ([(2, 3)], [1]),
            # (1, 2) -> (1, 3)
            ([(1, 2)], [(1, 3)]),
            # (1, 3) -> (1, 3)
            ([(1, 3)], [(1, 3)]),
            ([2, (1, 3)], [2, (1, 3)]),
            # (1, 2)-(2, 3) -> (1, 2)-(2, 3)
            ([(1, 2),(2, 3)], [(1, 2),(2, 3)]),
            # (0, 1) -> (0, 1) extension of null
            ([(0, 1)], [(0, 1)]),
            ([2, (0, 1)], [2, (0, 1)]),
            # (0, 1) -> (0, 2) extension of null -> (1, 2)
            ([2, (0, 1)], [1, (0, 2)]),
            # (0, 2) -> (0, 3) extension of (1, 2) -> (1, 3)
            ([(0, 2)], [(0, 3)]),
            # (0, 3) -> (0, 3) extension of (1, 3) -> (1, 3)
            ([(0, 3)], [(0, 3)]),
            ([2, (0, 3)], [2, (0, 3)]),
            # (0, 2)-(2, 3) -> (0, 2)-(2, 3) extension of (1, 2)-(2, 3) -> ...
            ([(0, 2),(2, 3)], [(0, 2),(2, 3)]),

            # Double patterns
            # (2, 3), (1, 2) -> (2, 3)
            ([(2, 3)], [(1, 2)], [(2, 3)]),
            # (2, 3), (1, 3) -> (2, 3)
            ([(2, 3)], [(1, 3)], [(2, 3)]),
            # Combining (2, 3), (1, 2) -> (2, 3) with (0, 1) -> (0, 2)
            # (2, 3), (0, 1)-(1, 2) -> (0, 2)-(2, 3)
            ([(2, 3)], [(0, 1),(1, 2)], [(0, 2),(2, 3)]),
        ]
        return patterns_raw

    @staticmethod
    def _short_underslide_down_middle():
        """Short underslide going up, in the middle of PMC."""
        # Local PMC is 0*-1-2-3-4*, with 1 and 3 paired.
        patterns_raw = [
            #### Single patterns
            # () -> ()
            ([], []), ([1], [1]), ([2], [2]), ([2], [1]),
            # (2, 3) -> ()
            ([(2, 3)], [1]),
            # (1, 2) -> (1, 3)
            ([(1, 2)], [(1, 3)]),
            # (1, 3) -> (1, 3)
            ([(1, 3)], [(1, 3)]),
            ([2, (1, 3)], [2, (1, 3)]),
            # (1, 2)-(2, 3) -> (1, 2)-(2, 3)
            ([(1, 2),(2, 3)], [(1, 2),(2, 3)]),
            # *** Lower extension ***
            # (0, 1) -> (0, 1) extension of null
            ([(0, 1)], [(0, 1)]),
            ([2, (0, 1)], [2, (0, 1)]),
            # (0, 1) -> (0, 2) extension of null -> (1, 2)
            ([2, (0, 1)], [1, (0, 2)]),
            # (0, 2) -> (0, 3) extension of (1, 2) -> (1, 3)
            ([(0, 2)], [(0, 3)]),
            # (0, 3) -> (0, 3) extension of (1, 3) -> (1, 3)
            ([(0, 3)], [(0, 3)]),
            ([2, (0, 3)], [2, (0, 3)]),
            # (0, 2)-(2, 3) -> (0, 2)-(2, 3) extension of (1, 2)-(2, 3) -> ...
            ([(0, 2),(2, 3)], [(0, 2),(2, 3)]),
            # *** Upper extension ***
            # (3, 4) -> (3, 4) extension of null
            ([(3, 4)], [(3, 4)]),
            ([2, (3, 4)], [2, (3, 4)]),
            # (2, 4) -> (3, 4) extension of (2, 3) -> ()
            ([(2, 4)], [(3, 4)]),
            # (1, 4) -> (1, 4) extension of (1, 3) -> (1, 3)
            ([(1, 4)], [(1, 4)]),
            ([2, (1, 4)], [2, (1, 4)]),
            # (1, 2)-(2, 4) -> (1, 2)-(2, 4) extension of (1, 2)-(2, 3) -> ...
            ([(1, 2),(2, 4)], [(1, 2),(2, 4)]),
            # *** Both extensions ***
            # (0, 1)-(3, 4) -> (0, 1)-(3, 4) extension of null
            ([(0, 1),(3, 4)], [(0, 1),(3, 4)]),
            ([2, (0, 1),(3, 4)], [2, (0, 1),(3, 4)]),
            # (0, 2)-(3, 4) -> (0, 3)-(3, 4) extension of (0, 2) -> (0, 3)
            ([(0, 2),(3, 4)], [(0, 3),(3, 4)]),
            # (0, 4) -> (0, 4) extension of (1, 3) -> (1, 3)
            ([(0, 4)], [(0, 4)]),
            ([1, (0, 4)], [1, (0, 4)]),
            ([2, (0, 4)], [2, (0, 4)]),
            ([2, (0, 4)], [1, (0, 4)]),
            ([1, 2, (0, 4)], [1, 2, (0, 4)]),
            # (0, 1)-(1, 4) -> (0, 1)-(1, 4) extension of (1, 3) -> (1, 3)
            ([(0, 1),(1, 4)], [(0, 1),(1, 4)]),
            ([2, (0, 1),(1, 4)], [2, (0, 1),(1, 4)]),
            ([(0, 1),(1, 2),(2, 4)], [(0, 1),(1, 2),(2, 4)]),
            # (0, 3)-(3, 4) -> (0, 3)-(3, 4) extension of (1, 3) -> (1, 3)
            ([(0, 3),(3, 4)], [(0, 3),(3, 4)]),
            ([2, (0, 3),(3, 4)], [2, (0, 3),(3, 4)]),
            ([(0, 2),(2, 3),(3, 4)], [(0, 2),(2, 3),(3, 4)]),
            # (0, 1)-(2, 4) -> (0, 1)-(3, 4) extension of (2, 4) -> (3, 4)
            ([(0, 1),(2, 4)], [(0, 1),(3, 4)]),
            # (0, 2)-(2, 4) -> (0, 3)-(3, 4) rather strange
            ([(0, 2),(2, 4)], [(0, 3),(3, 4)]),
            ([1, (0, 2),(2, 4)], [1, (0, 2),(2, 4)]),

            #### Double or triple patterns
            # (2, 3), (1, 2) -> (2, 3)
            ([(2, 3)], [(1, 2)], [(2, 3)]),
            # (2, 3), (1, 3) -> (2, 3)
            ([(2, 3)], [(1, 3)], [(2, 3)]),
            # *** Lower extension ***
            # Combining (2, 3), (1, 2) -> (2, 3) with (0, 1) -> (0, 2)
            # (2, 3), (0, 1)-(1, 2) -> (0, 2)-(2, 3)
            ([(2, 3)], [(0, 1),(1, 2)], [(0, 2),(2, 3)]),
            # *** Upper extension ***
            # Extension of (2, 3), (1, 2) -> (2, 3)
            # (2, 3)-(3, 4), (1, 2) -> (2, 3)-(3, 4)
            ([(2, 3),(3, 4)], [(1, 2)], [(2, 3),(3, 4)]),
            # (2, 4), (1, 2) -> (2, 4)
            ([1, (2, 4)], [(1, 2)], [1, (2, 4)]),
            # Extension of (2, 3), (0, 2) -> (2, 3)
            # (2, 3), (1, 4) -> (2, 4)
            ([(2, 3)], [(1, 4)], [(2, 4)]),
            # (2, 3)-(3, 4), (1, 3) -> (2, 3)-(3, 4)
            ([(2, 3),(3, 4)], [(1, 3)], [(2, 3),(3, 4)]),
            # (2, 4), (1, 3) -> (2, 4)
            ([1, (2, 4)], [(1, 3)], [1, (2, 4)]),
            # *** Both extensions ***
            # Output covers (0, 1)-(3, 4)
            ([2, (0, 1)], [2, (3, 4)], [(0, 1),(3, 4)]),
            ([2, (0, 1)], [(2, 3),(3, 4)], [(0, 1),(3, 4)]),
            # Output is (0, 1)-(2, 4):
            # -- Input covers (0, 4)
            ([2, (0, 4)], [(0, 1),(2, 4)]),
            ([2, (0, 1)], [(1, 2),(2, 4)], [(0, 1),(2, 4)]),
            ([(0, 2),(2, 3)], [2, (3, 4)], [(0, 1),(2, 4)]),
            ([(2, 3)], [(0, 1),(1, 2)], [2, (3, 4)], [(0, 1),(2, 4)]),
            # -- Input covers (0, 4)-(2, 3)
            ([(2, 3)], [(0, 1),(1, 4)], [(0, 1),(2, 4)]),
            ([(2, 3)], [(0, 3),(3, 4)], [(0, 1),(2, 4)]),
            ([(2, 4)], [(0, 3)], [(0, 1),(2, 4)]),
            ([(0, 1),(2, 4)], [(1, 3)], [(0, 1),(2, 4)]),
            ([(0, 2),(2, 4)], [(2, 3)], [(0, 1),(2, 4)]),
            ([(2, 3)], [(0, 1),(1, 2)], [(2, 3),(3, 4)], [(0, 1),(2, 4)]),
            # Output covers the whole interval
            # -- Input covers (0, 2)-(3, 4)
            ([1, (0, 2)], [2, (3, 4)], [1, (0, 4)]),
            ([(0, 1),(1, 2)], [2, (3, 4)], [(0, 1),(1, 4)]),
            # -- Input covers (0, 4)
            ([1, (2, 4)], [(0, 1),(1, 2)], [1, (0, 2),(2, 4)]),
            ([(2, 3),(3, 4)], [(0, 1),(1, 2)], [(0, 2),(2, 3),(3, 4)]),
            ([(0, 1),(1, 2)], [(2, 3),(3, 4)], [(0, 1),(1, 4)]),
            ([2, (0, 1)], [(1, 2),(2, 4)], [1, (0, 4)]),
            ([1, (0, 2)], [(2, 3),(3, 4)], [1, (0, 4)]),
            ([(0, 2),(2, 3)], [2, (3, 4)], [1, (0, 4)]),
            # -- Input covers (0, 4)-(2, 3)
            ([(2, 3)], [(0, 3),(3, 4)], [1, (0, 4)]),
            ([(2, 4)], [(0, 3)], [1, (0, 4)]),
            ([(0, 1),(2, 4)], [(1, 3)], [1, (0, 4)]),
            ([(0, 2),(2, 4)], [(2, 3)], [1, (0, 4)]),
        ]
        return patterns_raw

    @staticmethod
    def _length3_underslide_down_bottom():
        """Length 3 underslide going up, on the bottom of PMC."""
        # Local PMC is 0-1-2-3-4*, with 0 and 3 paired.
        patterns_raw = [
            #### Zero patterns
            ([(0, 1)],),
            ([2, (0, 1)],),

            #### Single patterns
            # Seeds
            ([], []), ([0], [0]), ([1], [2]), ([2], [1]), ([1, 2], [2, 1]),
            ([0, 2], [0, 1]), ([0, 1], [0, 2]), ([0, 1, 2], [0, 2, 1]),
            ([2], [0]), ([1, 2], [2, 0]),

            ([(2, 3)], [0]),
            ([1, (2, 3)], [0, 2]),
            ([(1, 2)], [(2, 3)]),
            ([2, (0, 1)], [0, (1, 2)]),

            # First round
            ([(0, 1), (2, 3)], [0, (1, 2)]),
            ([(1, 3)], [(2, 3)]),
            ([2, (1, 3)], [1, (2, 3)]),
            ([2, (0, 1)], [1, (0, 2)]),
            ([(0, 1)], [(0, 2)]),
            ([(2, 3)], [(0, 1)], [(1, 2)]),
            ([(2, 3)], [(0, 3)], [(1, 3)]),
            ([1, (2, 3)], [1, (0, 3)], [2, (1, 3)]),
            ([(0, 3)], [(0, 3)]),
            ([2, (0, 3)], [1, (0, 3)]),
            ([1, (0, 3)], [2, (0, 3)]),
            ([1, 2, (0, 3)], [1, 2, (0, 3)]),
            ([1, (2, 3)], [(0, 1), (1, 3)], [(1, 2), (2, 3)]),
            ([(2, 3)], [(0, 2)], [(1, 3)]),
            ([1, (2, 3)], [1, (0, 2)], [2, (1, 3)]),
            ([(0, 2)], [(0, 3)]),
            ([1, (0, 2)], [2, (0, 3)]),
            ([1, (2, 3)], [(0, 1), (1, 2)], [(1, 2), (2, 3)]),
            ([(0, 2), (2, 3)], [(0, 1), (1, 3)]),
            ([1, (0, 2), (2, 3)], [2, (0, 1), (1, 3)]),
            ([(1, 2), (2, 3)], [(0, 1), (2, 3)]),

            # Need to add
            ([(0, 1), (1, 3)], [(0, 2), (2, 3)]),
            ([2, (0, 1), (1, 3)], [1, (0, 2), (2, 3)]),
            ([(0, 1), (1, 2)], [(0, 2), (2, 3)]),
            ([(0, 1), (1, 2), (2, 3)], [(0, 1), (1, 2), (2, 3)]),

            # Optional
            ([0, (1, 2)], [(0, 1), (2, 3)]),
            ([(1, 2), (2, 3)], [2, (0, 1)], [(1, 2), (2, 3)]),
            ([(1, 2), (2, 3)], [(0, 1), (2, 3)], [(1, 2), (2, 3)]),

            # *** Upper extension ***
            # First run
            ([(3, 4)], [(3, 4)]),  # seed
            ([2, (3, 4)], [1, (3, 4)]),
            ([1, (3, 4)], [2, (3, 4)]),
            ([1, 2, (3, 4)], [1, 2, (3, 4)]),

            # Added by linear algebra
            ([1, (2, 3), (3, 4)], [1, (0, 2)], [2, (1, 3), (3, 4)]),
            ([(2, 3), (3, 4)], [(0, 2)], [(1, 3), (3, 4)]),
            ([(1, 2), (2, 4)], [(2, 3), (3, 4)]),
            ([1, (2, 3)], [(1, 3), (3, 4)], [0, (2, 4)]),
            ([0, (1, 4)], [0, (2, 4)]),
            ([0, 2, (1, 4)], [0, 1, (2, 4)]),
            ([(0, 1), (2, 4)], [(1, 2), (3, 4)]),
            ([2, (0, 1)], [(1, 2), (2, 4)], [0, (1, 4)]),
            ([1, 2, (0, 4)], [1, 2, (0, 4)]),
            ([1, (0, 4)], [2, (0, 4)]),
            ([(0, 4)], [(0, 4)]),
            ([2, (0, 4)], [1, (0, 4)]),
            ([(2, 3)], [(0, 4)], [(1, 4)]),
            ([1, (2, 3)], [1, (0, 4)], [2, (1, 4)]),
            ([(0, 1), (2, 3)], [(1, 3), (3, 4)], [0, (1, 4)]),
            ([(2, 3), (3, 4)], [(0, 1)], [(1, 2), (3, 4)]),
            ([0, 1, (2, 4)], [1, (0, 3)], [0, 2, (1, 4)]),
            ([0, (2, 4)], [(0, 3)], [0, (1, 4)]),
            ([0, 1, (2, 4)], [(0, 1), (1, 3)], [0, (1, 2), (2, 4)]),
            ([2, (0, 1), (1, 4)], [1, (0, 2), (2, 4)]),
            ([(0, 1), (1, 4)], [(0, 2), (2, 4)]),
            ([0, (2, 4)], [(0, 2)], [0, (1, 4)]),
            ([0, 1, (2, 4)], [1, (0, 2)], [0, 2, (1, 4)]),
            ([1, (2, 3)], [(1, 3), (3, 4)], [(0, 1)], [(1, 2), (2, 4)]),
            ([(1, 2), (2, 4)], [(2, 3)], [0, (2, 4)]),
            ([(2, 3), (3, 4)], [(0, 3)], [(1, 3), (3, 4)]),
            ([1, (2, 3), (3, 4)], [1, (0, 3)], [2, (1, 3), (3, 4)]),
            ([0, (1, 2), (2, 4)], [2, (0, 1)], [0, (1, 2), (2, 4)]),
            ([2, (1, 4)], [1, (2, 4)]),
            ([(1, 4)], [(2, 4)]),
            ([(0, 2), (2, 4)], [(0, 1), (1, 4)]),
            ([1, (0, 2), (2, 4)], [2, (0, 1), (1, 4)]),
            ([2, (1, 4)], [0, (2, 4)]),
            ([(1, 2), (2, 3), (3, 4)], [2, (0, 1)], [1, (2, 3)], [0, (1, 2), (2, 4)]),
            ([(1, 2), (2, 3), (3, 4)], [2, (0, 1)], [(1, 2), (2, 3), (3, 4)]),
            ([1, (2, 3), (3, 4)], [(0, 1), (1, 2)], [(1, 2), (2, 3), (3, 4)]),
            ([1, 2, (3, 4)], [(1, 2), (2, 3)], [(0, 1), (2, 3)], [0, (1, 2), (2, 4)]),
            ([1, (2, 4)], [(1, 3)], [0, (2, 4)]),
            ([2, (1, 3), (3, 4)], [(0, 1), (2, 3)], [0, (1, 2), (2, 4)]),
            ([(1, 2), (2, 4)], [(2, 3)], [(0, 1)], [(1, 2), (2, 4)]),
            ([(1, 2), (2, 3)], [2, (3, 4)], [0, (2, 4)]),
            ([(1, 2), (3, 4)], [(2, 3), (3, 4)]),
            ([(2, 4)], [(3, 4)]),
            ([1, (2, 4)], [2, (3, 4)]),
            ([(0, 1), (2, 4)], [(1, 2)], [0, (1, 4)]),
            ([(1, 2), (2, 3)], [(2, 3), (3, 4)], [(0, 1)], [(1, 2), (2, 4)]),
            ([(1, 3), (3, 4)], [(2, 3), (3, 4)]),
            ([2, (1, 3), (3, 4)], [1, (2, 3), (3, 4)]),
            ([1, (2, 3), (3, 4)], [(0, 1), (1, 3)], [(1, 2), (2, 3), (3, 4)]),
            ([0, 1, (2, 4)], [(0, 1), (1, 2)], [0, (1, 2), (2, 4)]),
            ([(0, 1), (2, 4)], [(1, 3)], [0, (1, 4)]),
            ([(0, 1), (2, 3)], [(1, 2), (3, 4)], [0, (1, 4)]),
            ([1, (2, 3), (3, 4)], [0, (1, 2)], [(0, 1), (2, 3)], [0, (1, 2), (2, 4)]),
            ([0, (1, 2)], [2, (3, 4)], [0, (2, 4)]),
            ([0, (1, 2)], [(2, 3), (3, 4)], [0, (2, 4)]),
            ([(0, 1), (1, 2), (2, 4)], [(0, 1), (1, 2), (2, 4)]),
            ([(1, 2), (2, 3)], [(0, 1), (2, 4)], [(1, 2), (2, 4)]),
            ([1, (2, 4)], [(1, 3)], [(0, 1)], [(1, 2), (2, 4)]),
            ([1, (2, 3)], [(0, 1), (1, 4)], [(1, 2), (2, 4)]),
        ]
        return patterns_raw

    @staticmethod
    def _length3_underslide_down_middle():
        """Length 3 underslide going up, in the middle of PMC."""
        # Local PMC is 0*-1-2-3-4-5*, with 1 and 4 paired.
        patterns_raw = [
            #### Zero patterns
            ([(1, 2)],),
            ([3, (1, 2)],),

            #### Single patterns
            # Seeds
            ([], []), ([1], [1]), ([2], [3]), ([3], [2]), ([2, 3], [3, 2]),
            ([1, 3], [1, 2]), ([1, 2], [1, 3]), ([1, 2, 3], [1, 3, 2]),
            ([3], [1]), ([2, 3], [3, 1]),

            ([(3, 4)], [1]),
            ([2, (3, 4)], [1, 3]),
            ([(2, 3)], [(3, 4)]),
            ([3, (1, 2)], [1, (2, 3)]),

            # First round
            ([(1, 2), (3, 4)], [1, (2, 3)]),
            ([(2, 4)], [(3, 4)]),
            ([3, (2, 4)], [2, (3, 4)]),
            ([3, (1, 2)], [2, (1, 3)]),
            ([(1, 2)], [(1, 3)]),
            ([(3, 4)], [(1, 2)], [(2, 3)]),
            ([(3, 4)], [(1, 4)], [(2, 4)]),
            ([2, (3, 4)], [2, (1, 4)], [3, (2, 4)]),
            ([(1, 4)], [(1, 4)]),
            ([3, (1, 4)], [2, (1, 4)]),
            ([2, (1, 4)], [3, (1, 4)]),
            ([2, 3, (1, 4)], [2, 3, (1, 4)]),
            ([2, (3, 4)], [(1, 2), (2, 4)], [(2, 3), (3, 4)]),
            ([(3, 4)], [(1, 3)], [(2, 4)]),
            ([2, (3, 4)], [2, (1, 3)], [3, (2, 4)]),
            ([(1, 3)], [(1, 4)]),
            ([2, (1, 3)], [3, (1, 4)]),
            ([2, (3, 4)], [(1, 2), (2, 3)], [(2, 3), (3, 4)]),
            ([(1, 3), (3, 4)], [(1, 2), (2, 4)]),
            ([2, (1, 3), (3, 4)], [3, (1, 2), (2, 4)]),
            ([(2, 3), (3, 4)], [(1, 2), (3, 4)]),

            # Need to add
            ([(1, 2), (2, 4)], [(1, 3), (3, 4)]),
            ([3, (1, 2), (2, 4)], [2, (1, 3), (3, 4)]),
            ([(1, 2), (2, 3)], [(1, 3), (3, 4)]),
            ([(1, 2), (2, 3), (3, 4)], [(1, 2), (2, 3), (3, 4)]),

            # Optional
            ([1, (2, 3)], [(1, 2), (3, 4)]),
            ([(2, 3), (3, 4)], [3, (1, 2)], [(2, 3), (3, 4)]),
            ([(2, 3), (3, 4)], [(1, 2), (3, 4)], [(2, 3), (3, 4)]),

            # *** Upper extension ***
            # First run
            ([(4, 5)], [(4, 5)]),  # seed
            ([3, (4, 5)], [2, (4, 5)]),
            ([2, (4, 5)], [3, (4, 5)]),
            ([2, 3, (4, 5)], [2, 3, (4, 5)]),

            # Added by linear algebra
            ([2, (3, 4), (4, 5)], [2, (1, 3)], [3, (2, 4), (4, 5)]),
            ([(3, 4), (4, 5)], [(1, 3)], [(2, 4), (4, 5)]),
            ([(2, 3), (3, 5)], [(3, 4), (4, 5)]),
            ([2, (3, 4)], [(2, 4), (4, 5)], [1, (3, 5)]),
            ([1, (2, 5)], [1, (3, 5)]),
            ([1, 3, (2, 5)], [1, 2, (3, 5)]),
            ([(1, 2), (3, 5)], [(2, 3), (4, 5)]),
            ([3, (1, 2)], [(2, 3), (3, 5)], [1, (2, 5)]),
            ([2, 3, (1, 5)], [2, 3, (1, 5)]),
            ([2, (1, 5)], [3, (1, 5)]),
            ([(1, 5)], [(1, 5)]),
            ([3, (1, 5)], [2, (1, 5)]),
            ([(3, 4)], [(1, 5)], [(2, 5)]),
            ([2, (3, 4)], [2, (1, 5)], [3, (2, 5)]),
            ([(1, 2), (3, 4)], [(2, 4), (4, 5)], [1, (2, 5)]),
            ([(3, 4), (4, 5)], [(1, 2)], [(2, 3), (4, 5)]),
            ([1, 2, (3, 5)], [2, (1, 4)], [1, 3, (2, 5)]),
            ([1, (3, 5)], [(1, 4)], [1, (2, 5)]),
            ([1, 2, (3, 5)], [(1, 2), (2, 4)], [1, (2, 3), (3, 5)]),
            ([3, (1, 2), (2, 5)], [2, (1, 3), (3, 5)]),
            ([(1, 2), (2, 5)], [(1, 3), (3, 5)]),
            ([1, (3, 5)], [(1, 3)], [1, (2, 5)]),
            ([1, 2, (3, 5)], [2, (1, 3)], [1, 3, (2, 5)]),
            ([2, (3, 4)], [(2, 4), (4, 5)], [(1, 2)], [(2, 3), (3, 5)]),
            ([(2, 3), (3, 5)], [(3, 4)], [1, (3, 5)]),
            ([(3, 4), (4, 5)], [(1, 4)], [(2, 4), (4, 5)]),
            ([2, (3, 4), (4, 5)], [2, (1, 4)], [3, (2, 4), (4, 5)]),
            ([1, (2, 3), (3, 5)], [3, (1, 2)], [1, (2, 3), (3, 5)]),
            ([3, (2, 5)], [2, (3, 5)]),
            ([(2, 5)], [(3, 5)]),
            ([(1, 3), (3, 5)], [(1, 2), (2, 5)]),
            ([2, (1, 3), (3, 5)], [3, (1, 2), (2, 5)]),
            ([3, (2, 5)], [1, (3, 5)]),
            ([(2, 3), (3, 4), (4, 5)], [3, (1, 2)], [2, (3, 4)], [1, (2, 3), (3, 5)]),
            ([(2, 3), (3, 4), (4, 5)], [3, (1, 2)], [(2, 3), (3, 4), (4, 5)]),
            ([2, (3, 4), (4, 5)], [(1, 2), (2, 3)], [(2, 3), (3, 4), (4, 5)]),
            ([2, 3, (4, 5)], [(2, 3), (3, 4)], [(1, 2), (3, 4)], [1, (2, 3), (3, 5)]),
            ([2, (3, 5)], [(2, 4)], [1, (3, 5)]),
            ([3, (2, 4), (4, 5)], [(1, 2), (3, 4)], [1, (2, 3), (3, 5)]),
            ([(2, 3), (3, 5)], [(3, 4)], [(1, 2)], [(2, 3), (3, 5)]),
            ([(2, 3), (3, 4)], [3, (4, 5)], [1, (3, 5)]),
            ([(2, 3), (4, 5)], [(3, 4), (4, 5)]),
            ([(3, 5)], [(4, 5)]),
            ([2, (3, 5)], [3, (4, 5)]),
            ([(1, 2), (3, 5)], [(2, 3)], [1, (2, 5)]),
            ([(2, 3), (3, 4)], [(3, 4), (4, 5)], [(1, 2)], [(2, 3), (3, 5)]),
            ([(2, 4), (4, 5)], [(3, 4), (4, 5)]),
            ([3, (2, 4), (4, 5)], [2, (3, 4), (4, 5)]),
            ([2, (3, 4), (4, 5)], [(1, 2), (2, 4)], [(2, 3), (3, 4), (4, 5)]),
            ([1, 2, (3, 5)], [(1, 2), (2, 3)], [1, (2, 3), (3, 5)]),
            ([(1, 2), (3, 5)], [(2, 4)], [1, (2, 5)]),
            ([(1, 2), (3, 4)], [(2, 3), (4, 5)], [1, (2, 5)]),
            ([2, (3, 4), (4, 5)], [1, (2, 3)], [(1, 2), (3, 4)], [1, (2, 3), (3, 5)]),
            ([1, (2, 3)], [3, (4, 5)], [1, (3, 5)]),
            ([1, (2, 3)], [(3, 4), (4, 5)], [1, (3, 5)]),
            ([(1, 2), (2, 3), (3, 5)], [(1, 2), (2, 3), (3, 5)]),
            ([(2, 3), (3, 4)], [(1, 2), (3, 5)], [(2, 3), (3, 5)]),
            ([2, (3, 5)], [(2, 4)], [(1, 2)], [(2, 3), (3, 5)]),
            ([2, (3, 4)], [(1, 2), (2, 5)], [(2, 3), (3, 5)]),

            ([(0, 1)], [(0, 1)]),  # seed
            ([3, (0, 1)], [2, (0, 1)]),
            ([2, (0, 1)], [3, (0, 1)]),
            ([2, 3, (0, 1)], [2, 3, (0, 1)]),

            # Added by linear algebra
            ([2, (3, 4), (4, 5)], [2, (0, 1), (1, 3)], [3, (0, 1), (1, 2), (2, 5)]),
            ([(3, 4), (4, 5)], [(0, 1), (1, 3)], [(0, 1), (1, 2), (2, 5)]),
            ([2, (3, 4)], [(0, 1), (1, 2), (2, 3)], [2, (3, 4), (4, 5)], [(0, 1), (2, 3), (3, 5)]),
            ([2, 3, (4, 5)], [2, 3, (0, 4)], [1, 2, 3, (0, 5)]),
            ([2, (4, 5)], [2, (0, 4)], [1, 3, (0, 5)]),
            ([(4, 5)], [(0, 4)], [1, (0, 5)]),
            ([3, (4, 5)], [3, (0, 4)], [1, 2, (0, 5)]),
            ([(2, 4)], [(0, 1), (1, 2)], [(0, 3), (3, 4)]),
            ([3, (2, 4)], [3, (0, 1), (1, 2)], [2, (0, 3), (3, 4)]),
            ([2, (3, 5)], [(2, 4)], [(0, 1), (1, 2)], [1, (0, 3), (3, 5)]),
            ([3, (2, 4), (4, 5)], [1, 3, (0, 2)], [1, 2, (0, 3), (3, 5)]),
            ([(2, 4), (4, 5)], [1, (0, 2)], [1, (0, 3), (3, 5)]),
            ([(0, 3), (3, 4), (4, 5)], [(0, 2), (2, 4), (4, 5)]),
            ([2, (0, 3), (3, 4), (4, 5)], [3, (0, 2), (2, 4), (4, 5)]),
            ([2, 3, (0, 4)], [2, 3, (0, 4)]),
            ([2, (0, 4)], [3, (0, 4)]),
            ([(0, 4)], [(0, 4)]),
            ([3, (0, 4)], [2, (0, 4)]),
            ([2, (0, 4), (4, 5)], [3, (0, 1), (1, 5)]),
            ([2, 3, (0, 4), (4, 5)], [2, 3, (0, 1), (1, 5)]),
            ([3, (0, 4), (4, 5)], [2, (0, 1), (1, 5)]),
            ([(0, 4), (4, 5)], [(0, 1), (1, 5)]),
            ([2, (3, 4)], [(0, 1), (1, 2), (2, 3)], [(0, 2), (2, 3), (3, 4)]),
            ([(3, 4)], [(0, 1), (1, 3)], [(0, 2), (2, 4)]),
            ([2, (3, 4)], [2, (0, 1), (1, 3)], [3, (0, 2), (2, 4)]),
            ([2, (3, 4)], [(0, 2), (2, 4), (4, 5)], [1, (0, 3), (3, 5)]),
            ([(3, 5)], [(0, 4)], [(0, 1), (2, 5)]),
            ([2, (3, 5)], [2, (0, 4)], [3, (0, 1), (2, 5)]),
            ([2, (3, 4)], [2, (0, 4), (4, 5)], [3, (0, 1), (2, 5)]),
            ([(3, 4)], [(0, 4), (4, 5)], [(0, 1), (2, 5)]),
            ([(1, 2), (2, 3)], [2, 3, (0, 1)], [2, (3, 4), (4, 5)], [1, (0, 3), (3, 5)]),
            ([2, 3, (0, 1), (4, 5)], [2, 3, (1, 4)], [1, 2, 3, (0, 5)]),
            ([2, (0, 1), (4, 5)], [2, (1, 4)], [1, 3, (0, 5)]),
            ([(0, 1), (4, 5)], [(1, 4)], [1, (0, 5)]),
            ([3, (0, 1), (4, 5)], [3, (1, 4)], [1, 2, (0, 5)]),
            ([(0, 2), (3, 5)], [(0, 3), (4, 5)]),
            ([3, (0, 1), (2, 5)], [2, (0, 1), (3, 5)]),
            ([(0, 1), (2, 5)], [(0, 1), (3, 5)]),
            ([(0, 1), (1, 2), (2, 5)], [(0, 1), (1, 3), (3, 5)]),
            ([3, (0, 1), (1, 2), (2, 5)], [2, (0, 1), (1, 3), (3, 5)]),
            ([2, (3, 4)], [(0, 2), (2, 3), (4, 5)], [1, (0, 3), (3, 5)]),
            ([2, (3, 4)], [2, (0, 1), (1, 5)], [3, (0, 1), (2, 5)]),
            ([(3, 4)], [(0, 1), (1, 5)], [(0, 1), (2, 5)]),
            ([2, (3, 4), (4, 5)], [(0, 1), (1, 2), (2, 3)], [(0, 1), (1, 2), (2, 3), (3, 5)]),
            ([3, (0, 4), (4, 5)], [1, (0, 2), (2, 5)]),
            ([2, 3, (0, 4), (4, 5)], [1, 3, (0, 2), (2, 5)]),
            ([(2, 3), (3, 4)], [3, (0, 1), (1, 2)], [2, 3, (4, 5)], [(0, 1), (2, 3), (3, 5)]),
            ([3, (0, 5)], [(0, 1), (2, 5)]),
            ([2, 3, (0, 5)], [3, (0, 1), (2, 5)]),
            ([2, (0, 3), (3, 5)], [2, (3, 4)], [1, 3, (0, 5)]),
            ([(0, 3), (3, 5)], [(3, 4)], [1, (0, 5)]),
            ([2, (0, 1), (1, 3)], [2, 3, (4, 5)], [3, (0, 1), (1, 5)]),
            ([(0, 1), (1, 3)], [3, (4, 5)], [(0, 1), (1, 5)]),
            ([2, 3, (4, 5)], [2, (3, 4)], [2, (1, 3)], [2, 3, (0, 1)], [1, 2, 3, (0, 5)]),
            ([3, (4, 5)], [(3, 4)], [(1, 3)], [3, (0, 1)], [1, 2, (0, 5)]),
            ([2, (0, 3)], [3, (0, 4)]),
            ([(0, 3)], [(0, 4)]),
            ([(0, 1), (2, 3)], [3, (4, 5)], [(0, 1), (3, 5)]),
            ([(0, 3), (3, 5)], [(3, 4)], [(0, 1), (2, 5)]),
            ([2, (0, 3), (3, 5)], [2, (3, 4)], [3, (0, 1), (2, 5)]),
            ([3, (0, 1), (4, 5)], [2, (0, 1), (4, 5)]),
            ([(0, 1), (4, 5)], [(0, 1), (4, 5)]),
            ([2, (0, 1), (4, 5)], [3, (0, 1), (4, 5)]),
            ([2, 3, (0, 1), (4, 5)], [2, 3, (0, 1), (4, 5)]),
            ([(0, 1), (2, 3), (3, 5)], [(0, 1), (1, 2), (3, 5)]),
            ([3, (0, 1)], [(1, 3), (3, 5)], [(0, 1), (2, 5)]),
            ([2, 3, (0, 1)], [2, (1, 3), (3, 5)], [3, (0, 1), (2, 5)]),
            ([(0, 1), (1, 5)], [(0, 1), (1, 5)]),
            ([3, (0, 1), (1, 5)], [2, (0, 1), (1, 5)]),
            ([2, 3, (0, 1), (1, 5)], [2, 3, (0, 1), (1, 5)]),
            ([2, (0, 1), (1, 5)], [3, (0, 1), (1, 5)]),
            ([(3, 4), (4, 5)], [(1, 3)], [3, (0, 1)], [1, (0, 2), (2, 5)]),
            ([2, (3, 4), (4, 5)], [2, (1, 3)], [2, 3, (0, 1)], [1, 3, (0, 2), (2, 5)]),
            ([3, (4, 5)], [(3, 4)], [(0, 1), (1, 3)], [2, (0, 4), (4, 5)]),
            ([2, 3, (4, 5)], [2, (3, 4)], [2, (0, 1), (1, 3)], [2, 3, (0, 4), (4, 5)]),
            ([(3, 4)], [(0, 1), (1, 3)], [(3, 4), (4, 5)], [(0, 1), (2, 5)]),
            ([2, (3, 4)], [2, (0, 1), (1, 3)], [2, (3, 4), (4, 5)], [3, (0, 1), (2, 5)]),
            ([(3, 4)], [(0, 1), (1, 2)], [1, (0, 3)]),
            ([(3, 4)], [(0, 1), (1, 3)], [3, (4, 5)], [(0, 1), (2, 5)]),
            ([2, (3, 4)], [2, (0, 1), (1, 3)], [2, 3, (4, 5)], [3, (0, 1), (2, 5)]),
            ([2, (3, 4)], [(0, 1), (1, 2), (2, 3)], [2, (3, 4), (4, 5)], [1, (0, 3), (3, 5)]),
            ([(0, 2), (2, 4), (4, 5)], [(0, 3), (3, 4), (4, 5)]),
            ([3, (0, 2), (2, 4), (4, 5)], [2, (0, 3), (3, 4), (4, 5)]),
            ([2, (4, 5)], [(2, 4)], [(0, 1), (1, 2)], [1, (0, 3), (3, 5)]),
            ([2, 3, (4, 5)], [3, (2, 4)], [3, (0, 1), (1, 2)], [1, 2, (0, 3), (3, 5)]),
            ([2, 3, (0, 1)], [1, 3, (0, 2)]),
            ([3, (0, 1)], [1, (0, 2)]),
            ([(3, 4)], [(0, 1), (1, 2)], [(2, 3), (4, 5)], [1, (0, 5)]),
            ([(0, 3), (4, 5)], [(3, 4)], [1, (0, 5)]),
            ([2, (0, 3), (4, 5)], [2, (3, 4)], [1, 3, (0, 5)]),
            ([2, 3, (0, 5)], [2, 3, (0, 5)]),
            ([2, (0, 5)], [3, (0, 5)]),
            ([(0, 5)], [(0, 5)]),
            ([3, (0, 5)], [2, (0, 5)]),
            ([3, (0, 1), (1, 2)], [2, 3, (4, 5)], [2, (0, 3), (4, 5)]),
            ([(0, 1), (1, 2)], [2, (4, 5)], [(0, 3), (4, 5)]),
            ([(0, 1), (1, 2), (2, 3)], [2, 3, (4, 5)], [(0, 1), (1, 3), (3, 5)]),
            ([3, (0, 2), (2, 5)], [1, (0, 3), (3, 5)]),
            ([(3, 4), (4, 5)], [1, (0, 2)], [1, (2, 3)], [1, 2, (0, 5)]),
            ([1, 2, (3, 5)], [(0, 1), (1, 2), (2, 3)], [1, (0, 2), (2, 3), (3, 5)]),
            ([(1, 3), (3, 4)], [3, (0, 1), (4, 5)], [1, (0, 2), (2, 5)]),
            ([2, (1, 3), (3, 4)], [2, 3, (0, 1), (4, 5)], [1, 3, (0, 2), (2, 5)]),
            ([3, (0, 1), (1, 2)], [2, 3, (4, 5)], [(2, 3), (3, 4)], [1, (0, 2), (2, 5)]),
            ([(0, 1), (3, 5)], [(1, 4)], [(0, 1), (2, 5)]),
            ([2, (0, 1), (3, 5)], [2, (1, 4)], [3, (0, 1), (2, 5)]),
            ([(2, 3), (3, 4)], [(3, 4), (4, 5)], [(0, 1), (1, 2)], [(0, 1), (2, 3), (3, 5)]),
            ([(0, 2), (2, 3), (3, 4)], [(0, 2), (2, 3), (3, 4)]),
            ([2, 3, (0, 1)], [(1, 2), (2, 3), (3, 5)], [1, (0, 3), (3, 5)]),
            ([2, (0, 1), (1, 3)], [2, (3, 4), (4, 5)], [1, (0, 3), (3, 5)]),
            ([(2, 4), (4, 5)], [(1, 2)], [2, (0, 1)], [1, (0, 3), (3, 5)]),
            ([3, (2, 4), (4, 5)], [3, (1, 2)], [2, 3, (0, 1)], [1, 2, (0, 3), (3, 5)]),
            ([(3, 4), (4, 5)], [(1, 2)], [(0, 1), (2, 3)], [1, (0, 2), (2, 5)]),
            ([2, (0, 3), (4, 5)], [3, (0, 4), (4, 5)]),
            ([(0, 3), (4, 5)], [(0, 4), (4, 5)]),
            ([2, 3, (4, 5)], [2, (3, 4)], [(0, 1), (1, 2), (2, 3)], [1, (0, 2), (2, 3), (3, 5)]),
            ([(0, 2), (3, 4), (4, 5)], [(0, 1), (2, 3), (4, 5)]),
            ([1, 3, (0, 2)], [(2, 3), (3, 4), (4, 5)], [1, (0, 2), (2, 5)]),
            ([2, (3, 4), (4, 5)], [1, (2, 3)], [3, (0, 1), (1, 2)], [1, (0, 2), (2, 3), (3, 5)]),
            ([3, (0, 1)], [(3, 4), (4, 5)], [(0, 1), (4, 5)]),
            ([2, 3, (0, 1)], [2, (3, 4), (4, 5)], [3, (0, 1), (4, 5)]),
            ([(2, 3), (4, 5)], [(3, 4)], [(0, 1), (1, 2)], [1, (0, 3), (3, 5)]),
            ([2, (3, 4), (4, 5)], [(1, 2), (2, 3)], [2, 3, (0, 1)], [1, (0, 2), (2, 3), (3, 5)]),
            ([(0, 2), (2, 3), (3, 4), (4, 5)], [(0, 2), (2, 3), (3, 4), (4, 5)]),
            ([2, 3, (0, 1)], [2, 3, (4, 5)], [3, (0, 1), (4, 5)]),
            ([3, (0, 1)], [3, (4, 5)], [(0, 1), (4, 5)]),
            ([2, (4, 5)], [(0, 1), (2, 3)], [1, (0, 2), (3, 5)]),
            ([(0, 1), (1, 2), (2, 3)], [1, 2, (3, 5)], [1, (0, 3), (3, 5)]),
            ([2, 3, (4, 5)], [2, (3, 4)], [2, (0, 1), (1, 3)], [1, 2, (0, 3), (3, 5)]),
            ([(0, 1), (1, 2)], [1, (2, 3)], [(3, 4), (4, 5)], [1, (0, 5)]),
            ([(0, 2), (2, 3), (3, 5)], [(0, 3), (3, 4), (4, 5)]),
            ([2, (3, 4)], [(0, 1), (1, 2), (2, 3)], [2, 3, (4, 5)], [(0, 1), (2, 3), (3, 5)]),
            ([3, (0, 1), (1, 5)], [1, (0, 2), (2, 5)]),
            ([2, 3, (0, 1), (1, 5)], [1, 3, (0, 2), (2, 5)]),
            ([(2, 3), (3, 5)], [(3, 4)], [(0, 1), (1, 2)], [1, (0, 3), (3, 5)]),
            ([(3, 5)], [(0, 4)], [1, (0, 5)]),
            ([2, (3, 5)], [2, (0, 4)], [1, 3, (0, 5)]),
            ([2, (3, 5)], [(0, 1), (2, 3)], [1, (0, 2), (3, 5)]),
            ([(1, 3), (3, 5)], [3, (0, 1)], [1, (0, 2), (2, 5)]),
            ([2, (1, 3), (3, 5)], [2, 3, (0, 1)], [1, 3, (0, 2), (2, 5)]),
            ([3, (0, 5)], [1, (0, 5)]),
            ([2, 3, (0, 5)], [1, 3, (0, 5)]),
            ([(0, 1), (1, 2), (2, 3)], [2, 3, (4, 5)], [2, (3, 4)], [1, (0, 3), (3, 5)]),
            ([3, (2, 4), (4, 5)], [3, (0, 1), (1, 2)], [1, (0, 2), (2, 3), (3, 5)]),
            ([3, (1, 2)], [2, 3, (0, 1)], [1, 2, (0, 3)]),
            ([(1, 2)], [2, (0, 1)], [1, (0, 3)]),
            ([2, (0, 1), (4, 5)], [1, (2, 3)], [1, (0, 2), (3, 5)]),
            ([(0, 2), (3, 4)], [(0, 1), (2, 3)]),
            ([3, (4, 5)], [(3, 4)], [1, (0, 3)], [1, 2, (0, 5)]),
            ([2, 3, (4, 5)], [2, (3, 4)], [1, 2, (0, 3)], [1, 2, 3, (0, 5)]),
            ([2, 3, (0, 1)], [(2, 3), (3, 4), (4, 5)], [(1, 2), (3, 4)], [(0, 1), (2, 3), (3, 5)]),
            ([1, (2, 3)], [(0, 2), (3, 4), (4, 5)], [1, (0, 3), (3, 5)]),
            ([1, (0, 2)], [(2, 3), (4, 5)], [1, (0, 5)]),
            ([1, 2, (0, 3)], [2, 3, (4, 5)], [1, 3, (0, 5)]),
            ([1, (0, 3)], [3, (4, 5)], [1, (0, 5)]),
            ([3, (1, 2)], [(0, 1), (2, 3), (3, 5)], [1, (0, 2), (2, 5)]),
            ([(0, 3), (3, 4)], [3, (4, 5)], [(0, 1), (2, 5)]),
            ([2, (0, 3), (3, 4)], [2, 3, (4, 5)], [3, (0, 1), (2, 5)]),
            ([3, (0, 2), (2, 5)], [(0, 1), (2, 3), (3, 5)]),
            ([(0, 3), (3, 4)], [(0, 2), (2, 4)]),
            ([2, (0, 3), (3, 4)], [3, (0, 2), (2, 4)]),
            ([(0, 1), (2, 3), (4, 5)], [(1, 2), (3, 4)], [1, (0, 3), (3, 5)]),
            ([(2, 3), (3, 4)], [3, (0, 2), (4, 5)], [1, (0, 3), (3, 5)]),
            ([(0, 2), (2, 3), (3, 5)], [2, (3, 4)], [(0, 1), (2, 3), (3, 5)]),
            ([(2, 3), (4, 5)], [3, (0, 1)], [1, (0, 2), (3, 5)]),
            ([3, (0, 1), (1, 2)], [3, (2, 4), (4, 5)], [1, (0, 2), (2, 5)]),
            ([3, (4, 5)], [(3, 4)], [(1, 2)], [(0, 1), (2, 3)], [1, 2, (0, 5)]),
            ([(0, 1), (2, 3), (4, 5)], [(0, 1), (1, 2), (3, 5)]),
            ([2, 3, (0, 1)], [(1, 2), (2, 3), (3, 5)], [(0, 1), (2, 3), (3, 5)]),
            ([(3, 4), (4, 5)], [(0, 1), (1, 2)], [1, (2, 3)], [2, (0, 4), (4, 5)]),
            ([3, (0, 1)], [(1, 3), (3, 5)], [1, (0, 5)]),
            ([2, 3, (0, 1)], [2, (1, 3), (3, 5)], [1, 3, (0, 5)]),
            ([(0, 2), (2, 3), (4, 5)], [(0, 1), (1, 3), (3, 5)]),
            ([1, 2, (3, 5)], [2, (0, 1), (1, 3)], [1, 3, (0, 2), (2, 5)]),
            ([1, (3, 5)], [(0, 1), (1, 3)], [1, (0, 2), (2, 5)]),
            ([2, (3, 4)], [(0, 1), (1, 2), (2, 5)], [(0, 1), (2, 3), (3, 5)]),
            ([2, 3, (0, 1)], [(2, 3), (3, 4), (4, 5)], [(0, 2), (3, 4), (4, 5)]),
            ([2, 3, (0, 4), (4, 5)], [1, 2, (0, 3), (3, 5)]),
            ([2, (0, 4), (4, 5)], [1, (0, 3), (3, 5)]),
            ([2, 3, (4, 5)], [2, (3, 4)], [2, (0, 1), (1, 3)], [2, 3, (0, 1), (1, 5)]),
            ([3, (4, 5)], [(3, 4)], [(0, 1), (1, 3)], [2, (0, 1), (1, 5)]),
            ([1, 3, (0, 2), (2, 5)], [1, 2, (0, 3), (3, 5)]),
            ([1, (0, 2), (2, 5)], [1, (0, 3), (3, 5)]),
            ([(0, 2), (4, 5)], [(0, 3), (4, 5)]),
            ([3, (0, 2), (4, 5)], [2, (0, 3), (4, 5)]),
            ([(0, 3), (3, 5)], [(0, 4), (4, 5)]),
            ([2, (0, 3), (3, 5)], [3, (0, 4), (4, 5)]),
            ([(2, 3), (3, 4)], [(0, 2), (3, 4), (4, 5)], [(0, 1), (2, 3), (3, 5)]),
            ([2, (4, 5)], [(0, 2), (2, 3)], [1, (0, 3), (3, 5)]),
            ([(1, 2), (3, 4)], [(0, 1), (2, 3), (4, 5)], [1, (0, 2), (2, 5)]),
            ([(2, 3), (3, 5)], [3, (0, 1)], [1, (0, 2), (3, 5)]),
            ([2, (3, 4), (4, 5)], [1, 2, (0, 3)], [1, 3, (0, 2), (2, 5)]),
            ([(3, 4), (4, 5)], [1, (0, 3)], [1, (0, 2), (2, 5)]),
            ([3, (0, 2)], [2, (0, 3)]),
            ([(0, 2)], [(0, 3)]),
            ([(0, 1), (1, 2)], [2, (4, 5)], [(2, 3)], [1, (0, 5)]),
            ([3, (0, 1), (1, 2)], [1, (2, 3), (3, 5)], [1, (0, 2), (2, 5)]),
            ([2, (0, 1), (3, 5)], [2, (1, 4)], [1, 3, (0, 5)]),
            ([(0, 1), (3, 5)], [(1, 4)], [1, (0, 5)]),
            ([(3, 4), (4, 5)], [(0, 1), (1, 2)], [1, (2, 3)], [2, (0, 1), (1, 5)]),
            ([(0, 2), (3, 4)], [1, (0, 3)]),
            ([1, (3, 5)], [(0, 1), (1, 2)], [1, (2, 3)], [1, 2, (0, 5)]),
            ([1, (2, 3)], [3, (0, 1), (1, 2)], [2, (3, 4), (4, 5)], [1, (0, 3), (3, 5)]),
            ([(2, 3), (3, 5)], [(0, 2), (3, 4)], [(0, 1), (2, 3), (3, 5)]),
            ([(0, 2), (3, 4)], [(2, 4), (4, 5)], [(0, 1), (2, 5)]),
            ([2, (0, 1), (1, 3), (3, 5)], [3, (0, 2), (2, 4), (4, 5)]),
            ([(0, 1), (1, 3), (3, 5)], [(0, 2), (2, 4), (4, 5)]),
            ([3, (0, 1), (1, 2)], [2, 3, (4, 5)], [(0, 1), (2, 3), (4, 5)]),
            ([(0, 1), (1, 2), (3, 5)], [(0, 1), (2, 3), (4, 5)]),
            ([2, (3, 4)], [(2, 4), (4, 5)], [(0, 1), (1, 2)], [(0, 1), (2, 3), (3, 5)]),
            ([3, (0, 1), (1, 2)], [3, (2, 4), (4, 5)], [2, (0, 4), (4, 5)]),
            ([(0, 1), (1, 2)], [(2, 4), (4, 5)], [(0, 4), (4, 5)]),
            ([3, (0, 1), (1, 2)], [(2, 3), (3, 4), (4, 5)], [(0, 1), (1, 2), (2, 5)]),
            ([(0, 2), (2, 3), (3, 4)], [2, 3, (4, 5)], [(0, 1), (2, 3), (3, 5)]),
            ([3, (0, 2)], [(0, 1), (2, 3)]),
            ([2, (3, 5)], [(0, 2), (2, 4)], [1, (0, 3), (3, 5)]),
            ([1, (0, 2)], [(2, 4), (4, 5)], [1, (0, 5)]),
            ([1, 3, (0, 2)], [3, (2, 4), (4, 5)], [1, 2, (0, 5)]),
            ([2, (3, 5)], [(0, 2), (2, 4)], [(0, 1), (2, 3), (3, 5)]),
            ([(2, 3), (3, 4), (4, 5)], [3, (0, 1), (1, 2)], [(0, 2), (2, 3), (3, 4), (4, 5)]),
            ([3, (0, 2)], [1, (0, 3)]),
            ([(0, 1), (1, 2)], [(2, 3), (4, 5)], [(0, 1), (1, 5)]),
            ([(2, 3), (3, 4)], [3, (0, 1), (1, 2)], [2, (3, 4), (4, 5)], [(0, 1), (2, 3), (3, 5)]),
            ([(0, 1), (1, 2), (2, 3), (3, 5)], [(0, 2), (2, 3), (3, 4), (4, 5)]),
            ([(0, 2), (2, 3), (3, 5)], [2, (3, 4)], [1, (0, 3), (3, 5)]),
            ([(0, 2), (2, 3), (3, 4)], [2, 3, (4, 5)], [1, (0, 3), (3, 5)]),
            ([(0, 1), (1, 3)], [(3, 4), (4, 5)], [(0, 1), (1, 5)]),
            ([2, (0, 1), (1, 3)], [2, (3, 4), (4, 5)], [3, (0, 1), (1, 5)]),
            ([(2, 3), (3, 5)], [(0, 2), (3, 4)], [1, (0, 3), (3, 5)]),
            ([(3, 4), (4, 5)], [(0, 1), (1, 2)], [(0, 1), (2, 3), (4, 5)]),
            ([2, (0, 3), (3, 4)], [2, 3, (4, 5)], [1, 3, (0, 5)]),
            ([(0, 3), (3, 4)], [3, (4, 5)], [1, (0, 5)]),
            ([(2, 3), (3, 4)], [3, (0, 1), (1, 2)], [(0, 2), (2, 3), (3, 4)]),
            ([(0, 2), (2, 3), (4, 5)], [2, (3, 4)], [1, (0, 3), (3, 5)]),
            ([2, (3, 4)], [(0, 2), (2, 4), (4, 5)], [(0, 1), (2, 3), (3, 5)]),
            ([(0, 1), (2, 3), (3, 5)], [(1, 2), (3, 4)], [(0, 1), (2, 3), (3, 5)]),
            ([(3, 4), (4, 5)], [(1, 2)], [2, (0, 1)], [1, (2, 3)], [1, 2, (0, 5)]),
            ([(0, 1), (1, 2)], [(2, 3), (4, 5)], [(3, 4)], [1, (0, 5)]),
            ([(0, 2), (4, 5)], [(2, 4)], [1, (0, 5)]),
            ([3, (0, 2), (4, 5)], [3, (2, 4)], [1, 2, (0, 5)]),
            ([(0, 1), (2, 3), (3, 5)], [3, (1, 2)], [1, (0, 3), (3, 5)]),
            ([2, (3, 4), (4, 5)], [1, (0, 2), (2, 3)], [1, (0, 2), (2, 3), (3, 5)]),
            ([3, (2, 4), (4, 5)], [3, (0, 1), (1, 2)], [2, (0, 1), (1, 3), (3, 5)]),
            ([(2, 4), (4, 5)], [(0, 1), (1, 2)], [(0, 1), (1, 3), (3, 5)]),
            ([2, (0, 1), (3, 5)], [3, (0, 1), (4, 5)]),
            ([(0, 1), (3, 5)], [(0, 1), (4, 5)]),
            ([(0, 1), (2, 3)], [(3, 4), (4, 5)], [(0, 1), (3, 5)]),
            ([2, (0, 1), (3, 5)], [(1, 2), (2, 4)], [(0, 1), (2, 3), (3, 5)]),
            ([3, (0, 1), (1, 2)], [1, 3, (2, 5)], [1, 2, (0, 5)]),
            ([(0, 1), (1, 2)], [1, (2, 5)], [1, (0, 5)]),
            ([(2, 3), (4, 5)], [3, (0, 2)], [1, (0, 3), (3, 5)]),
            ([3, (0, 1), (1, 2), (2, 5)], [1, (0, 2), (2, 3), (3, 5)]),
            ([(0, 2), (3, 5)], [(2, 4)], [1, (0, 5)]),
            ([2, (3, 5)], [(2, 4)], [(0, 1), (1, 2)], [(0, 1), (2, 3), (3, 5)]),
            ([(0, 2), (2, 5)], [(0, 3), (3, 5)]),
            ([3, (0, 2), (2, 5)], [2, (0, 3), (3, 5)]),
            ([2, 3, (0, 1)], [(2, 3), (3, 4), (4, 5)], [(1, 2), (3, 4)], [1, (0, 3), (3, 5)]),
            ([(0, 1), (1, 2), (2, 3)], [2, (3, 4), (4, 5)], [(0, 3), (3, 4), (4, 5)]),
            ([2, (0, 1), (3, 5)], [1, (2, 3)], [1, (0, 2), (3, 5)]),
            ([(1, 2)], [2, (0, 1)], [(2, 3), (4, 5)], [1, (0, 5)]),
            ([(2, 3), (3, 5)], [(3, 4)], [(0, 1), (1, 2)], [(0, 1), (2, 3), (3, 5)]),
            ([(0, 2), (2, 3)], [(0, 3), (3, 4)]),
            ([(3, 4)], [(0, 1), (1, 2)], [(0, 1), (2, 3)]),
            ([2, (0, 3), (4, 5)], [1, (0, 3), (3, 5)]),
            ([3, (1, 2)], [2, 3, (0, 1)], [(2, 3), (3, 4), (4, 5)], [1, (0, 2), (2, 5)]),
            ([(0, 2), (3, 5)], [(2, 4)], [(0, 1), (2, 5)]),
            ([(1, 2), (2, 3), (3, 4)], [2, 3, (0, 1), (4, 5)], [1, (0, 2), (2, 3), (3, 5)]),
            ([(2, 3), (3, 4)], [(0, 1), (1, 2), (3, 5)], [(0, 1), (2, 3), (3, 5)]),
            ([2, (0, 1), (4, 5)], [(1, 2), (2, 3)], [1, (0, 3), (3, 5)]),
            ([3, (0, 2), (4, 5)], [(0, 1), (2, 3), (4, 5)]),
            ([3, (0, 1), (1, 2)], [2, (3, 4), (4, 5)], [1, (2, 3)], [1, (0, 2), (2, 5)]),
            ([2, 3, (0, 1)], [1, (2, 3), (3, 5)], [1, (0, 2), (3, 5)]),
            ([3, (0, 2), (2, 4)], [2, (0, 3), (3, 4)]),
            ([(0, 2), (2, 4)], [(0, 3), (3, 4)]),
            ([(1, 2), (2, 3), (3, 5)], [2, 3, (0, 1)], [1, (0, 2), (2, 3), (3, 5)]),
            ([(0, 1), (2, 3)], [(0, 2), (3, 4)]),
            ([(0, 2), (3, 4)], [(2, 4), (4, 5)], [1, (0, 5)]),
            ([3, (0, 1), (1, 2)], [2, (3, 4), (4, 5)], [(0, 1), (2, 3), (4, 5)]),
            ([2, (0, 1), (3, 5)], [(1, 2), (2, 4)], [1, (0, 3), (3, 5)]),
        ]
        return patterns_raw

    @staticmethod
    def _general_underslide_down_middle():
        """Underslide of length > 3 going up, in the middle of PMC."""
        # Local PMC at left (D-side) is 0*-1-2-3*, 4*-5-6*, with 1 and 5 paired.
        # Local PMC at right (A-side) is 0*-1-2*, 3*-4-5-6*, with 1 and 5
        # paired.
        patterns_raw = [
            #### Initial patterns
            ([(1, 2)],),
            ([], []), ([1], [1]), ([4], [2]), ([4], [1]), ([1, 4], [1, 2]),
            ([(4, 5)], [1]),

            #### Seed for top
            ([(5, 6)], [(5, 6)]),
            ([4, (5, 6)], [2, (5, 6)]),
            # Linear algebra produces
            ([(4, 6)], [(5, 6)]),

            #### Seed for bottom
            ([(0, 1)], [(0, 1)]),
            ([4, (0, 1)], [2, (0, 1)]),
            # Linear algebra produces
            ([4, (0, 1)], [4, (5, 6)], [(0, 1), (5, 6)]),
            ([4, (0, 1)], [(4, 5), (5, 6)], [(0, 1), (5, 6)]),
            ([(0, 1), (5, 6)], [(0, 1), (5, 6)]),
            ([4, (0, 1), (5, 6)], [2, (0, 1), (5, 6)]),
            ([(0, 1), (4, 6)], [(0, 1), (5, 6)]),
            ([4, (0, 1)], [1, (0, 2)]),

            #### Seed for upper middle
            ([(3, 4)], [(4, 5)]),
            # Linear algebra produces
            ([1, (3, 4)], [4, (0, 1), (5, 6)], [1, (0, 2), (4, 6)]),
            ([(3, 4), (4, 6)], [(4, 5), (5, 6)]),
            ([1, (3, 4)], [(4, 5), (5, 6)], [1, (4, 6)]),
            ([(4, 5)], [(3, 5), (5, 6)], [1, (4, 6)]),
            ([(0, 1), (3, 4), (5, 6)], [(0, 2), (4, 5), (5, 6)]),
            ([(0, 1), (3, 4)], [(0, 2), (4, 5)]),
            ([1, (3, 4)], [4, (5, 6)], [1, (4, 6)]),
            ([(4, 6)], [(3, 5)], [1, (4, 6)]),
            ([(0, 1), (3, 4)], [(4, 5), (5, 6)], [(0, 1), (4, 6)]),
            ([(3, 4), (5, 6)], [(4, 5), (5, 6)]),
            ([4, (0, 1), (3, 6)], [2, (0, 1), (4, 6)]),
            ([(0, 1), (3, 6)], [(0, 1), (4, 6)]),
            ([4, (3, 6)], [1, (4, 6)]),
            ([(3, 4), (4, 5)], [4, (5, 6)], [1, (4, 6)]),
            ([4, (0, 1)], [(3, 4), (4, 5), (5, 6)], [(0, 2), (4, 5), (5, 6)]),
            ([1, (3, 4)], [(1, 2), (4, 5)]),
            ([(0, 1), (3, 4), (4, 6)], [(0, 2), (4, 5), (5, 6)]),
            ([(0, 1), (3, 4)], [4, (5, 6)], [(0, 1), (4, 6)]),
            ([1, 4, (3, 6)], [1, 2, (4, 6)]),
            ([1, (3, 6)], [1, (4, 6)]),
            ([(3, 5), (5, 6)], [(4, 5), (5, 6)]),
            ([4, (3, 5), (5, 6)], [2, (4, 5), (5, 6)]),
            ([(3, 6)], [(4, 6)]),
            ([4, (3, 6)], [2, (4, 6)]),
            ([(3, 4), (4, 5)], [4, (0, 1), (5, 6)], [1, (0, 2), (4, 6)]),
            ([(3, 4), (4, 5)], [(1, 2), (4, 5)]),
            ([(3, 4), (4, 6)], [(4, 5)], [1, (4, 6)]),
            ([4, (0, 1), (3, 6)], [1, (0, 2), (4, 6)]),
            ([4, (3, 5)], [2, (4, 5)]),
            ([(3, 5)], [(4, 5)]),

            #### Seed for lower middle
            ([4, (1, 2)], [1, (2, 3)]),
            # Linear algebra produces
            ([4, (1, 2), (3, 6)], [2, (1, 3), (4, 6)]),
            ([(1, 2), (3, 6)], [(1, 3), (4, 6)]),
            ([(0, 2), (4, 5)], [(3, 4), (5, 6)], [1, (0, 3), (4, 6)]),
            ([4, (0, 1), (1, 2)], [4, (5, 6)], [2, (0, 3), (5, 6)]),
            ([(0, 1), (1, 2)], [(5, 6)], [(0, 3), (5, 6)]),
            ([1, (3, 4)], [(4, 5), (5, 6)], [(0, 1), (1, 2)], [1, (0, 3), (4, 6)]),
            ([(0, 1), (1, 2), (4, 6)], [(0, 1), (2, 3), (5, 6)]),
            ([(4, 5)], [(1, 2), (3, 6)], [(2, 3), (4, 6)]),
            ([(4, 6)], [(0, 2), (3, 4)], [1, (0, 3), (4, 6)]),
            ([4, (0, 1), (1, 2)], [4, (5, 6)], [(3, 4), (4, 5)], [1, (0, 2), (2, 3), (4, 6)]),
            ([(4, 5)], [(1, 2), (3, 4)], [(2, 3), (4, 5)]),
            ([1, (0, 2)], [(3, 5), (5, 6)], [1, (0, 3), (4, 6)]),
            ([1, 4, (0, 2)], [4, (3, 5), (5, 6)], [1, 2, (0, 3), (4, 6)]),
            ([(3, 4), (4, 6)], [(0, 2), (4, 5)], [(0, 1), (2, 3), (4, 6)]),
            ([1, (3, 6)], [(0, 1), (1, 2)], [1, (0, 3), (4, 6)]),
            ([1, 4, (3, 6)], [4, (0, 1), (1, 2)], [1, 2, (0, 3), (4, 6)]),
            ([4, (0, 2), (3, 5)], [2, (0, 3), (4, 5)]),
            ([(0, 2), (3, 5)], [(0, 3), (4, 5)]),
            ([(1, 2), (3, 4), (4, 6)], [4, (0, 1)], [1, (0, 2), (2, 3), (4, 6)]),
            ([(5, 6)], [(0, 2), (3, 5)], [1, (0, 3), (4, 6)]),
            ([4, (5, 6)], [4, (0, 2), (3, 5)], [1, 2, (0, 3), (4, 6)]),
            ([(0, 2), (3, 5), (5, 6)], [(0, 1), (1, 3), (4, 6)]),
            ([4, (0, 2), (3, 5), (5, 6)], [2, (0, 1), (1, 3), (4, 6)]),
            ([(1, 2), (4, 6)], [(3, 5)], [1, (2, 3), (4, 6)]),
            ([(4, 5)], [(0, 2), (3, 4), (5, 6)], [1, (0, 3), (4, 6)]),
            ([(4, 5)], [(1, 2)], [(2, 3)]),
            ([(1, 2), (3, 4), (4, 6)], [(4, 5)], [1, (2, 3), (4, 6)]),
            ([(0, 2), (3, 4), (4, 5)], [(0, 2), (2, 3), (4, 5)]),
            ([(4, 5), (5, 6)], [(0, 1), (1, 2)], [1, (3, 4)], [1, (0, 2), (2, 3), (4, 6)]),
            ([(3, 4), (4, 5), (5, 6)], [4, (1, 2)], [(2, 3), (4, 5), (5, 6)]),
            ([4, (1, 2)], [4, (0, 1)], [(3, 4), (4, 5), (5, 6)], [1, (0, 2), (2, 3), (4, 6)]),
            ([(4, 6)], [(0, 2), (3, 5)], [1, (0, 3), (4, 6)]),
            ([4, (0, 1), (1, 2)], [(3, 4), (4, 5), (5, 6)], [(0, 1), (1, 2), (2, 3), (4, 6)]),
            ([1, (0, 2), (3, 4)], [4, (5, 6)], [1, (0, 3), (4, 6)]),
            ([(4, 5)], [(0, 1), (1, 2), (3, 4)], [(0, 2), (2, 3), (4, 5)]),
            ([(4, 5), (5, 6)], [(1, 2), (3, 5)], [(2, 3), (4, 5), (5, 6)]),
            ([(4, 5), (5, 6)], [(0, 1), (1, 2), (3, 4)], [(0, 1), (1, 2), (2, 3), (4, 6)]),
            ([(0, 2), (3, 4), (5, 6)], [(0, 3), (4, 5), (5, 6)]),
            ([(0, 1), (1, 2)], [1, (3, 4)], [(4, 5), (5, 6)], [1, (0, 3), (4, 6)]),
            ([(0, 1), (3, 4), (4, 6)], [4, (1, 2)], [1, (0, 3), (4, 6)]),
            ([(4, 5), (5, 6)], [(1, 2)], [(2, 3), (5, 6)]),
            ([4, (0, 2), (5, 6)], [2, (0, 3), (5, 6)]),
            ([(0, 2), (5, 6)], [(0, 3), (5, 6)]),
            ([(0, 1), (1, 2), (3, 4)], [4, (5, 6)], [(0, 1), (1, 3), (4, 6)]),
            ([(1, 2), (3, 4), (4, 5)], [(1, 2), (2, 3), (4, 5)]),
            ([4, (0, 2)], [(3, 4), (4, 6)], [1, (0, 3), (4, 6)]),
            ([4, (1, 2)], [2, (1, 3)]),
            ([(1, 2)], [(1, 3)]),
            ([(4, 6)], [(3, 5)], [(0, 1), (1, 2)], [(0, 1), (2, 3), (4, 6)]),
            ([(3, 4), (4, 5)], [4, (0, 1), (1, 2)], [(4, 5), (5, 6)], [(0, 1), (2, 3), (4, 6)]),
            ([4, (0, 2)], [2, (0, 3)]),
            ([(0, 2)], [(0, 3)]),
            ([(1, 2), (4, 5)], [1, (2, 3)]),
            ([4, (0, 1), (1, 2)], [1, 4, (3, 6)], [1, 2, (0, 3), (4, 6)]),
            ([(0, 1), (1, 2)], [1, (3, 6)], [1, (0, 3), (4, 6)]),
            ([(1, 2), (3, 4), (4, 5)], [4, (5, 6)], [1, (2, 3), (4, 6)]),
            ([(1, 2), (3, 4), (4, 6)], [(2, 3), (4, 5), (5, 6)]),
            ([(1, 2), (4, 6)], [(2, 3), (5, 6)]),
            ([(3, 4), (4, 5)], [(0, 2), (4, 5), (5, 6)], [(0, 1), (2, 3), (4, 6)]),
            ([(4, 6)], [(0, 2), (3, 5)], [(0, 1), (2, 3), (4, 6)]),
            ([(1, 2), (4, 5)], [(0, 1), (3, 4), (5, 6)], [1, (0, 2), (2, 3), (4, 6)]),
            ([(0, 2), (4, 5)], [1, (0, 3)]),
            ([(3, 4), (4, 5)], [(1, 2), (4, 5)], [(2, 3), (4, 5)]),
            ([4, (0, 1), (1, 2)], [4, (3, 5), (5, 6)], [2, (0, 3), (4, 5), (5, 6)]),
            ([(0, 1), (1, 2)], [(3, 5), (5, 6)], [(0, 3), (4, 5), (5, 6)]),
            ([(4, 5)], [(0, 1), (1, 2)], [1, (0, 3)]),
            ([(0, 2), (4, 5)], [(3, 5), (5, 6)], [(0, 1), (2, 3), (4, 6)]),
            ([1, (3, 4), (4, 6)], [(1, 2), (4, 5)], [1, (2, 3), (4, 6)]),
            ([(1, 2), (3, 4)], [(1, 3), (4, 5)]),
            ([4, (0, 2), (5, 6)], [4, (3, 5)], [1, 2, (0, 3), (4, 6)]),
            ([(0, 2), (5, 6)], [(3, 5)], [1, (0, 3), (4, 6)]),
            ([(4, 5), (5, 6)], [(1, 2), (3, 4)], [4, (0, 1)], [1, (0, 2), (2, 3), (4, 6)]),
            ([(4, 6)], [(3, 5)], [(0, 1), (1, 2)], [1, (0, 3), (4, 6)]),
            ([(0, 1), (4, 6)], [(1, 2), (3, 4)], [1, (0, 3), (4, 6)]),
            ([4, (5, 6)], [(0, 2), (3, 4), (4, 5)], [1, (0, 2), (2, 3), (4, 6)]),
            ([(4, 5)], [(3, 5), (5, 6)], [(0, 1), (1, 2)], [1, (0, 3), (4, 6)]),
            ([(3, 4), (4, 5)], [(0, 1), (1, 2), (4, 6)], [(0, 1), (2, 3), (4, 6)]),
            ([(0, 2), (4, 5)], [(0, 1), (2, 3)]),
            ([4, (0, 1), (1, 2)], [4, (5, 6)], [(0, 1), (2, 3), (5, 6)]),
            ([(0, 2), (3, 4), (4, 5), (5, 6)], [(0, 1), (1, 2), (2, 3), (4, 6)]),
            ([(3, 4), (4, 5)], [(1, 2), (4, 6)], [(2, 3), (4, 6)]),
            ([(1, 2), (4, 6)], [(0, 1), (3, 4)], [1, (0, 2), (2, 3), (4, 6)]),
            ([4, (0, 1), (1, 2)], [4, (3, 5), (5, 6)], [1, (0, 2), (2, 3), (4, 6)]),
            ([(4, 5)], [(0, 1), (1, 2), (3, 4)], [(4, 5), (5, 6)], [(0, 1), (2, 3), (4, 6)]),
            ([4, (0, 2), (3, 6)], [(0, 1), (2, 3), (4, 6)]),
            ([(3, 5), (5, 6)], [(0, 1), (1, 2)], [(0, 3), (4, 5), (5, 6)]),
            ([4, (3, 5), (5, 6)], [4, (0, 1), (1, 2)], [2, (0, 3), (4, 5), (5, 6)]),
            ([4, (0, 2), (3, 6)], [2, (0, 3), (4, 6)]),
            ([(0, 2), (3, 6)], [(0, 3), (4, 6)]),
            ([(0, 2), (3, 4)], [(0, 3), (4, 5)]),
            ([(0, 2), (4, 5)], [(3, 5), (5, 6)], [1, (0, 3), (4, 6)]),
            ([4, (1, 2)], [4, (0, 1)], [1, 2, (0, 3)]),
            ([(1, 2)], [(0, 1)], [1, (0, 3)]),
            ([(4, 5)], [(1, 2), (3, 5)], [(2, 3), (4, 5)]),
            ([(4, 5)], [(0, 2), (3, 5), (5, 6)], [(0, 1), (2, 3), (4, 6)]),
            ([(4, 5), (5, 6)], [(0, 1), (1, 2)], [(0, 1), (2, 3), (5, 6)]),
            ([4, (0, 2)], [1, (0, 3)]),
            ([(0, 1), (4, 6)], [(1, 2), (3, 5)], [(0, 1), (2, 3), (4, 6)]),
            ([(3, 4), (4, 5)], [4, (1, 2)], [(2, 3), (4, 5)]),
            ([(4, 5)], [(3, 5), (5, 6)], [(0, 1), (1, 2)], [(0, 1), (2, 3), (4, 6)]),
            ([(3, 4), (4, 6)], [(4, 5)], [(1, 2)], [(2, 3), (4, 6)]),
            ([(3, 4), (4, 6)], [(4, 5)], [(0, 1), (1, 2)], [1, (0, 3), (4, 6)]),
            ([(0, 2), (3, 4), (5, 6)], [(4, 5)], [1, (0, 3), (4, 6)]),
            ([(0, 2), (4, 6)], [(0, 3), (5, 6)]),
            ([4, (5, 6)], [(4, 5)], [(0, 1), (1, 2), (3, 4)], [1, (0, 2), (2, 3), (4, 6)]),
            ([4, (3, 5)], [4, (0, 1), (1, 2)], [2, (0, 3), (4, 5)]),
            ([(3, 5)], [(0, 1), (1, 2)], [(0, 3), (4, 5)]),
            ([(3, 4), (4, 5), (5, 6)], [(1, 2), (4, 5)], [(2, 3), (4, 5), (5, 6)]),
            ([(3, 4), (4, 5), (5, 6)], [4, (1, 2)], [4, (0, 1)], [1, (0, 2), (2, 3), (4, 6)]),
            ([(0, 2), (4, 5), (5, 6)], [1, (3, 4)], [1, (0, 2), (2, 3), (4, 6)]),
            ([(0, 1), (4, 6)], [(1, 2), (3, 5)], [1, (0, 3), (4, 6)]),
            ([4, (0, 1), (1, 2)], [(4, 5), (5, 6)], [1, (3, 4)], [1, (0, 2), (2, 3), (4, 6)]),
            ([4, (0, 1)], [(1, 2), (3, 4), (4, 6)], [1, (0, 3), (4, 6)]),
            ([4, (0, 2), (5, 6)], [(3, 4), (4, 5)], [1, (0, 2), (2, 3), (4, 6)]),
            ([(4, 5), (5, 6)], [(1, 2), (3, 4)], [(2, 3), (4, 5), (5, 6)]),
            ([(3, 4), (4, 5)], [4, (0, 1), (1, 2)], [4, (5, 6)], [(0, 1), (2, 3), (4, 6)]),
            ([4, (5, 6)], [(3, 4), (4, 5)], [4, (0, 1), (1, 2)], [1, (0, 2), (2, 3), (4, 6)]),
            ([(0, 1), (1, 2), (3, 4)], [(4, 5), (5, 6)], [(0, 1), (1, 3), (4, 6)]),
            ([(0, 1), (1, 2)], [1, (3, 4)], [4, (5, 6)], [1, (0, 3), (4, 6)]),
            ([(3, 4), (4, 5)], [(4, 5), (5, 6)], [(1, 2)], [(2, 3), (4, 6)]),
            ([4, (1, 2)], [(0, 1), (3, 4), (4, 6)], [1, (0, 2), (2, 3), (4, 6)]),
            ([(0, 1), (3, 4), (4, 6)], [(1, 2), (4, 5)], [(0, 1), (2, 3), (4, 6)]),
            ([(0, 2), (3, 4), (4, 6)], [(4, 5)], [(0, 1), (2, 3), (4, 6)]),
            ([(3, 4), (4, 5)], [(4, 5), (5, 6)], [(0, 1), (1, 2)], [(0, 1), (2, 3), (4, 6)]),
            ([4, (0, 1), (5, 6)], [(1, 2), (3, 4), (4, 5)], [1, (0, 2), (2, 3), (4, 6)]),
            ([(0, 1), (3, 4), (4, 6)], [(1, 2), (4, 5)], [1, (0, 3), (4, 6)]),
            ([(3, 4), (4, 5), (5, 6)], [4, (0, 1), (1, 2)], [(0, 1), (1, 2), (2, 3), (4, 6)]),
            ([4, (1, 2), (3, 5)], [2, (1, 3), (4, 5)]),
            ([(1, 2), (3, 5)], [(1, 3), (4, 5)]),
            ([(3, 4), (4, 6)], [(4, 5)], [(0, 1), (1, 2)], [(0, 1), (2, 3), (4, 6)]),
            ([(1, 2), (4, 5)], [(3, 5), (5, 6)], [1, (2, 3), (4, 6)]),
            ([(0, 2), (3, 4), (4, 5)], [4, (5, 6)], [(0, 1), (2, 3), (4, 6)]),
            ([(0, 2), (3, 4), (4, 6)], [(0, 1), (1, 3), (4, 6)]),
            ([(0, 1), (3, 4), (5, 6)], [(1, 2), (4, 5)], [1, (0, 3), (4, 6)]),
            ([4, (1, 2), (3, 6)], [1, (2, 3), (4, 6)]),
            ([4, (0, 1), (1, 2)], [(4, 5), (5, 6)], [(0, 1), (2, 3), (5, 6)]),
            ([1, (3, 4)], [4, (0, 2), (5, 6)], [1, (0, 3), (4, 6)]),
            ([(0, 1), (1, 2), (3, 6)], [(0, 1), (1, 3), (4, 6)]),
            ([4, (0, 1), (1, 2), (3, 6)], [2, (0, 1), (1, 3), (4, 6)]),
            ([(0, 2), (4, 6)], [(3, 5)], [1, (0, 3), (4, 6)]),
            ([(3, 4), (4, 5)], [4, (0, 1), (1, 2)], [(0, 2), (2, 3), (4, 5)]),
            ([(0, 1), (1, 2), (3, 4), (4, 6)], [(0, 1), (1, 2), (2, 3), (4, 6)]),
            ([4, (0, 2), (5, 6)], [(0, 1), (2, 3), (5, 6)]),
            ([(0, 1), (1, 2), (4, 6)], [1, (3, 4)], [1, (0, 2), (2, 3), (4, 6)]),
            ([(3, 4), (5, 6)], [(0, 2), (4, 5)], [1, (0, 3), (4, 6)]),
            ([1, (4, 6)], [(1, 2), (3, 5)], [1, (2, 3), (4, 6)]),
            ([(0, 1), (1, 2)], [(3, 4), (5, 6)], [(0, 3), (4, 5), (5, 6)]),
            ([(4, 6)], [(3, 5)], [(1, 2)], [(2, 3), (4, 6)]),
            ([(3, 5), (5, 6)], [1, (0, 2)], [1, (0, 3), (4, 6)]),
            ([4, (3, 5), (5, 6)], [1, 4, (0, 2)], [1, 2, (0, 3), (4, 6)]),
            ([4, (0, 2)], [(0, 1), (2, 3)]),
            ([4, (3, 5), (5, 6)], [4, (0, 1), (1, 2)], [1, (0, 2), (2, 3), (4, 6)]),
            ([4, (0, 1), (5, 6)], [4, (1, 2), (3, 5)], [1, 2, (0, 3), (4, 6)]),
            ([(0, 1), (5, 6)], [(1, 2), (3, 5)], [1, (0, 3), (4, 6)]),
            ([(3, 4), (4, 6)], [4, (0, 2)], [1, (0, 3), (4, 6)]),
            ([(0, 2), (4, 5), (5, 6)], [(0, 1), (2, 3), (5, 6)]),
            ([(3, 4), (4, 6)], [(0, 2), (4, 5)], [1, (0, 3), (4, 6)]),
            ([(4, 5)], [(3, 5), (5, 6)], [(1, 2)], [(2, 3), (4, 6)]),
            ([(4, 5)], [(0, 1), (1, 2), (3, 4)], [4, (5, 6)], [(0, 1), (2, 3), (4, 6)]),
            ([(4, 5), (5, 6)], [(1, 2)], [(0, 1), (3, 4)], [1, (0, 2), (2, 3), (4, 6)]),
            ([(0, 2), (4, 6)], [(3, 4)], [1, (0, 3), (4, 6)]),
            ([(0, 2), (4, 6)], [(3, 5)], [(0, 1), (2, 3), (4, 6)]),
            ([(4, 5)], [(0, 1), (1, 2), (3, 6)], [(0, 1), (2, 3), (4, 6)]),
            ([(4, 5)], [(0, 1), (1, 2)], [(0, 1), (2, 3)]),
        ]
        return patterns_raw

    @staticmethod
    def getDerivedTwoStepArrows(arrows_base, arrow_new, mod_gens):
        """Find all ways of deriving two-step arrows from arrow_new, including
        by combining it with elements of arrows_base.

        """
        result = []
        coeff_d, coeffs_a = arrow_new
        if len(coeffs_a) > 0 and coeffs_a[0].isIdempotent():
            return result
        # Take anti-differential of one of coeffs_a
        for i in range(len(coeffs_a)):
            for anti_da, coeff in coeffs_a[i].antiDiff().items():
                result.append(
                    (coeff_d, coeffs_a[:i] + (anti_da,) + coeffs_a[i+1:]))
            for (a, b), coeff in coeffs_a[i].factor().items():
                if a.isIdempotent() or b.isIdempotent():
                    continue
                result.append((coeff_d, coeffs_a[:i] + (a, b) + coeffs_a[i+1:]))
        # Take differential of coeff_d
        for dd, coeff in coeff_d.diff().items():
            result.append((dd, coeffs_a))
        # Multiply two together
        for coeff_d2, coeffs_a2 in arrows_base:
            if len(coeffs_a2) > 0 and coeffs_a2[0].isIdempotent():
                continue
            # One direction
            if coeff_d * coeff_d2 != E0:
                if any([ArcslideDA.idemMatchDA(x, y, coeff_d, coeffs_a) and
                        ArcslideDA.idemMatchDA(y, z, coeff_d2, coeffs_a2)
                        for x, y, z in itertools.product(
                                mod_gens, mod_gens, mod_gens)]):
                    result.append(((coeff_d * coeff_d2).getElt(),
                                   coeffs_a + coeffs_a2))
            # Other direction
            if coeff_d2 * coeff_d != E0:
                if any([ArcslideDA.idemMatchDA(x, y, coeff_d2, coeffs_a2) and
                        ArcslideDA.idemMatchDA(y, z, coeff_d, coeffs_a)
                        for x, y, z in itertools.product(
                                mod_gens, mod_gens, mod_gens)]):
                    result.append(((coeff_d2 * coeff_d).getElt(),
                                   coeffs_a2 + coeffs_a))
        return result

    @staticmethod
    def getAltFactorizations(arrows_base, arrow_new, mod_gens):
        result = []
        coeff_d, coeffs_a = arrow_new
        # Take differential of one of coeffs_a
        for i in range(len(coeffs_a)):
            for da, coeff in coeffs_a[i].diff().items():
                result.append((coeff_d, coeffs_a[:i] + (da,) + coeffs_a[i+1:]))
            if i > 0 and coeffs_a[i-1] * coeffs_a[i] != E0:
                result.append((coeff_d, coeffs_a[:i-1] +
                               ((coeffs_a[i-1] * coeffs_a[i]).getElt(),) +
                               coeffs_a[i+1:]))
        # Take anti-differential of coeff_d
        for anti_dd, coeff in coeff_d.antiDiff().items():
            result.append((anti_dd, coeffs_a))
        # Split into two sequences
        for (a, b), coeff in coeff_d.factor().items():
            # Number of A-inputs in the first sequence
            for i in range(len(coeffs_a)):
                if (a, coeffs_a[:i]) in arrows_base:
                    arrow_used = False
                    for x, y in itertools.product(mod_gens, mod_gens):
                        if ArcslideDA.idemMatchDA(x, y, a, coeffs_a[:i]):
                            arrow_used = True
                    if arrow_used:
                        result.append((b, coeffs_a[i:]))
                if (b, coeffs_a[i:]) in arrows_base:
                    result.append((a, coeffs_a[:i]))
        return result

    @staticmethod
    def getAltIdempotents(arrows, single_idems):
        """Arrows is a list of tuples (coeff_d, coeffs_a).
        single_idems is a list of tuples (idem_d, idem_a), specifying the ID of
        single idempotents on the D-side and A-side.
        Returns a list of arrows that are different from those in the input by
        only the single idempotents.

        """
        def uses_idempotent(alg_gen, idem):
            return idem in alg_gen.left_idem or idem in alg_gen.right_idem

        def has_singlehor(alg_gen, idem):
            return idem in alg_gen.single_hor

        def add_singlehor(alg_gen, idem):
            return LocalStrandDiagram(
                alg_gen.parent, [idem] + list(alg_gen.left_idem),
                alg_gen.strands)

        def remove_singlehor(alg_gen, idem):
            return LocalStrandDiagram(
                alg_gen.parent, [i for i in alg_gen.left_idem if i != idem],
                alg_gen.strands)

        results = []
        for idem_d, idem_a in single_idems:
            for coeff_d, coeffs_a in results + arrows:
                if has_singlehor(coeff_d, idem_d) and \
                   all([has_singlehor(coeff, idem_a) for coeff in coeffs_a]):
                    new_arrow = (remove_singlehor(coeff_d, idem_d),
                                 tuple([remove_singlehor(coeff, idem_a)
                                        for coeff in coeffs_a]))
                    if not new_arrow in results + arrows:
                        results.append(new_arrow)
                if (not uses_idempotent(coeff_d, idem_d)) and \
                   all([not uses_idempotent(coeff, idem_a)
                        for coeff in coeffs_a]):
                    new_arrow = (add_singlehor(coeff_d, idem_d),
                                 tuple([add_singlehor(coeff, idem_a)
                                        for coeff in coeffs_a]))
                    if not new_arrow in results + arrows:
                        results.append(new_arrow)
        return results

    @staticmethod
    def autoCompleteByLinAlg(arrows_base, arrows_new, single_idems, mod_gens):
        """Auto-complete arrows by solving a system of linear equations mod 2.
        Use only when it is clear that no two added arrows can be composed
        together.

        Returns the list of suggested arrows.

        """
        # First step: find all possible arrows, and all possible two-step
        # arrows.
        one_step_arrows = set()
        two_step_arrows = set()
        one_step_arrows_queue = Queue()
        two_step_arrows_queue = Queue()
        for arrow in arrows_new:
            one_step_arrows_queue.put(arrow)
            one_step_arrows.add(arrow)

        while not one_step_arrows_queue.empty() or \
              not two_step_arrows_queue.empty():
            if not one_step_arrows_queue.empty():
                cur_arrow = one_step_arrows_queue.get()
                for arrow in ArcslideDA.getDerivedTwoStepArrows(
                        arrows_base, cur_arrow, mod_gens):
                    if arrow not in two_step_arrows:
                        two_step_arrows.add(arrow)
                        two_step_arrows_queue.put(arrow)
                for arrow in ArcslideDA.getAltIdempotents(
                        [cur_arrow], single_idems):
                    if arrow not in one_step_arrows:
                        one_step_arrows.add(arrow)
                        one_step_arrows_queue.put(arrow)
            else:
                cur_arrow = two_step_arrows_queue.get()
                for arrow in ArcslideDA.getAltFactorizations(
                        arrows_base, cur_arrow, mod_gens):
                    if arrow not in one_step_arrows:
                        one_step_arrows.add(arrow)
                        one_step_arrows_queue.put(arrow)
        print len(one_step_arrows), len(two_step_arrows)

        # Combine some of the one-step arrows
        combined_one_step_arrows = []  # list of lists of arrows
        while len(one_step_arrows) != 0:
            arrow = one_step_arrows.pop()

            # Temporary measure: forbid cases that don't agree on the center
            # idempotent (4 on the A-side and 2 (idem 1) on the D-side).
            middle_a, middle_d = 4, 1
            coeff_d, coeffs_a = arrow
            start1 = [id for id in range(len(coeffs_a))
                      if any([middle_a == p for p, q in coeffs_a[id].strands])]
            end1 = [id for id in range(len(coeffs_a))
                    if any([middle_a == q for p, q in coeffs_a[id].strands])]
            if len(start1) == 1 and len(end1) == 1 and \
               start1[0] < end1[0] and middle_d in coeff_d.single_hor:
                continue

            alt_idems = ArcslideDA.getAltIdempotents([arrow], single_idems)
            for alt_idem in alt_idems:
                one_step_arrows.remove(alt_idem)
            combined_one_step_arrows.append([arrow] + alt_idems)

        print len(combined_one_step_arrows)

        # Generate the matrix mapping from one-step arrows to two-step arrows
        two_step_arrows = list(two_step_arrows)
        matrix_map = [[0] * len(two_step_arrows)
                      for i in range(len(combined_one_step_arrows))]
        target_vec = [0] * len(two_step_arrows)
        for i in range(len(combined_one_step_arrows)):
            for one_step_arrow in combined_one_step_arrows[i]:
                derived_two_steps = ArcslideDA.getDerivedTwoStepArrows(
                    arrows_base, one_step_arrow, mod_gens)
                for j in range(len(two_step_arrows)):
                    if two_step_arrows[j] in derived_two_steps:
                        if one_step_arrow in arrows_new:
                            target_vec[j] += 1
                            target_vec[j] %= 2
                        else:
                            matrix_map[i][j] += 1
                            matrix_map[i][j] %= 2
        print "Constructing system..."
        lin_sys = F2RowSystem(matrix_map)
        comb, reduced_vec = lin_sys.vecReduce(target_vec)
        for i in range(len(reduced_vec)):
            if reduced_vec[i] == 1:
                print "Missing: ", two_step_arrows[i]

        result = []
        for i in range(len(combined_one_step_arrows)):
            if comb[i] != 0 and \
               combined_one_step_arrows[i][0] not in arrows_new:
                result.extend(combined_one_step_arrows[i])
        return result

    @staticmethod
    def autoCompleteArrows(arrows, single_idems, mod_gens):
        """Arrows is a list of tuples (coeff_d, coeffs_a)."""
        # First see if there are idempotent inconsistencies
        alt_idempotents = ArcslideDA.getAltIdempotents(arrows, single_idems)
        if len(alt_idempotents) > 0:
            return alt_idempotents[0]

        # Mapping (coeff_d, coeffs_a) to the number of times it appears as a two
        # step arrow.
        two_step_arrows = {}
        def add_to_step(coeff_d, coeffs_a):
            if (coeff_d, coeffs_a) not in two_step_arrows:
                two_step_arrows[(coeff_d, coeffs_a)] = 0
            two_step_arrows[(coeff_d, coeffs_a)] += 1

        for i in range(len(arrows)):
            for coeff_d, coeffs_a in ArcslideDA.getDerivedTwoStepArrows(
                    arrows[:i], arrows[i], mod_gens):
                add_to_step(coeff_d, coeffs_a)

        def findAlternateFactorizations(coeff_d, coeffs_a):
            return [(alt_d, alt_a) for alt_d, alt_a in
                    ArcslideDA.getAltFactorizations(
                        arrows, (coeff_d, coeffs_a), mod_gens)
                    if (alt_d, alt_a) not in arrows
                    # and alt_d.multiplicity[3] != 0  # temporary measure
                ]

        problems = []
        # Mapping from suggested (coeff_d, coeffs_a) to number of times it is
        # suggested.
        suggestions_count = {}
        for (coeff_d, coeffs_a), count in two_step_arrows.items():
            if count % 2 == 0:
                continue
            problems.append((coeff_d, coeffs_a))
            alts = findAlternateFactorizations(coeff_d, coeffs_a)
            assert len(alts) > 0, "No way out of %s %s" % \
                (str(coeffs_a), str(coeff_d))
            for alt in alts:  # coeff_d, coeffs_a
                if alt not in suggestions_count:
                    suggestions_count[alt] = 0
                suggestions_count[alt] += 1
            if len(alts) == 1:
                print "Resolving: ", coeff_d, coeffs_a
                return alts[0]

        if len(problems) == 0:
            return (True,)

        # Add the most-suggested arrow
        # suggestions = []
        # for (coeff_d, coeffs_a), times in suggestions_count.items():
        #     suggestions.append((times, coeff_d, coeffs_a))
        # times, coeff_d, coeffs_a = max(suggestions)
        # print "Taking chances with: ", times, coeff_d, coeffs_a
        # return (coeff_d, coeffs_a)
        return (len(problems) == 0, problems, suggestions_count)
