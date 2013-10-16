"""Producing type DA structures for arcslides, using local actions."""

from dastructure import SimpleDAGenerator, SimpleDAStructure
from dastructure import AddChordToDA
from localpmc import LocalIdempotent, LocalStrandAlgebra, LocalStrandDiagram
from localpmc import restrictPMC, restrictStrandDiagram
from pmc import Strands
from utility import subset
from utility import F2

class ArcslideDA:
    """Responsible for producing a type DA structure for an arcslide, using
    local actions.

    """
    def __init__(self, slide):
        """Specifies the arcslide to use. slide should be of type Arcslide.
        In addition to recording slide, construct the following:
        short_chord - chord on the D-side for the action with no A-side input.
        short_chord_local - corresponding local strand diagram on the D-side.
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

        if b1 == 1 and c1 == 2 and c2 == 0:
            local_cut1, outer_cut1, local_cut2, outer_cut2 = (
                [(0, 2)], [(3, n-1)],
                [(0, 2)], [(3, n-1)])
            self.short_chord = [(0, 1)]
            short_chord_local_raw = [(0, 1)]
            patterns_raw = ArcslideDA._short_underslide_up_bottom()
        elif b1 == n-2 and c1 == n-1 and c2 == n-3:
            local_cut1, outer_cut1, local_cut2, outer_cut2 = (
                [(n-3, n-1)], [(0, n-4)],
                [(n-3, n-1)], [(0, n-4)])
            self.short_chord = [(n-3, n-2)]
            short_chord_local_raw = [(1, 2)]
            patterns_raw = ArcslideDA._short_underslide_up_top()
        elif b1 == c1 - 1 and c2 == c1 - 2:  # *-c2-b1-c1-*
            assert c2 != 0 and c1 != n-1
            local_cut1, outer_cut1, local_cut2, outer_cut2 = (
                [(c2, c1)], [(0, c2-1), (c1+1, n-1)],
                [(c2, c1)], [(0, c2-1), (c1+1, n-1)])
            self.short_chord = [(c2, b1)]
            short_chord_local_raw = [(1, 2)]
            patterns_raw = ArcslideDA._short_underslide_up_middle()
        else:
            raise NotImplementedError(
                "This slide pattern is not yet implemented.")

        self.local_pmc1, self.mapping1 = restrictPMC(pmc1, local_cut1)
        self.outer_pmc1, self.outer_mapping1 = restrictPMC(pmc1, outer_cut1)
        self.local_pmc2, self.mapping2 = restrictPMC(pmc2, local_cut2)
        self.outer_pmc2, self.outer_mapping2 = restrictPMC(pmc2, outer_cut2)
        self.short_chord_local = self.local_pmc2.sd(short_chord_local_raw)

        # Required so the left to right transition on the outside can proceed.
        assert self.outer_pmc1 == self.outer_pmc2

        self.arrow_patterns = {}
        for pattern in patterns_raw:
            key = []
            for i in range(len(pattern)-1):
                key.append(self.local_pmc1.sd(pattern[i]))
            key = tuple(key)
            if key not in self.arrow_patterns:
                self.arrow_patterns[key] = []
            self.arrow_patterns[key].append(self.local_pmc2.sd(pattern[-1]))

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

        alg1_gens = alg1.getGenerators()
        mod_gens = dastr.getGenerators()

        # Add action with zero algebra input
        AddChordToDA(dastr, Strands(pmc2, self.short_chord), [])

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
            for cur_a in alg1_gens:
                if cur_a.isIdempotent():
                    continue
                if len(cur_list) > 0 and \
                   cur_list[-1].right_idem != cur_a.left_idem:
                    continue
                new_list = cur_list + (cur_a,)
                new_list_local = cur_list_local + (
                    restrictStrandDiagram(
                        pmc1, cur_a, self.local_pmc1, self.mapping1),)
                outer_a = restrictStrandDiagram(
                    pmc1, cur_a, self.outer_pmc1, self.outer_mapping1)
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
                if new_list_local in self.arrow_patterns:
                    for local_d in self.arrow_patterns[new_list_local]:
                        alg_d = local_d.join(new_prod_d.removeSingleHor(), pmc2,
                                             self.mapping2, self.outer_mapping2)
                        if alg_d is None:
                            continue
                        for x in mod_gens:
                            for y in mod_gens:
                                if x.idem1 == alg_d.left_idem and \
                                   x.idem2 == new_list[0].left_idem and \
                                   y.idem1 == alg_d.right_idem and \
                                   y.idem2 == new_list[-1].right_idem:
                                    dastr.addDelta(x, y, alg_d, new_list, 1)
                # Local patterns match the prefix of one of the arrows
                if any([new_list_local == pattern[0:len(new_list)]
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
            if coeffs_a[0].isIdempotent():
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
                    # print "Warning: unused arrow: ", coeffs_a, coeff_d
                    pass
        for x in mod_gens:
            for y in mod_gens:
                if x.idem2 == y.idem2 and \
                   x.idem1 == self.short_chord_local.left_idem and \
                   y.idem1 == self.short_chord_local.right_idem:
                    local_dastr.addDelta(x, y, self.short_chord_local, [], 1)
        return local_dastr

    # The next series of functions specify the local arrows. The format is as
    # follows:
    # All but the last element of the tuple is a list to be passed to the sd()
    # function of the A-side local PMC, specifying the A-side inputs. The last
    # element of the tuple is a list to be passed to the sd() function of the
    # D-side local PMC, specifying the D-side output.
    @staticmethod
    def _short_underslide_up_bottom():
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
            ([(1, 2)], [1, (0, 2)], [(1, 2)]),
            # *** Extensions of (1, 2), (0, 1) -> (1, 2) ***
            # (1, 2)-(2, 3), (0, 1) -> (1, 2)-(2, 3)
            ([(1, 2),(2, 3)], [(0, 1)], [(1, 2),(2, 3)]),
            # (1, 3), (0, 1) -> (1, 3)
            ([0, (1, 3)], [(0, 1)], [0, (1, 3)]),
            # *** Extensions of (1, 2), (0, 2) -> (1, 2) ***
            # (1, 2), (0, 3) -> (1, 3)
            ([(1, 2)], [(0, 3)], [(1, 3)]),
            ([(1, 2)], [1, (0, 3)], [(1, 3)]),
            # (1, 2)-(2, 3), (0, 2) -> (1, 2)-(2, 3)
            ([(1, 2),(2, 3)], [(0, 2)], [(1, 2),(2, 3)]),
            ([(1, 2),(2, 3)], [1, (0, 2)], [(1, 2),(2, 3)]),
            # (1, 3), (0, 2) -> (1, 3)
            ([0, (1, 3)], [(0, 2)], [0, (1, 3)]),
            ([0, (1, 3)], [1, (0, 2)], [0, (1, 3)]),
        ]
        return patterns_raw

    @staticmethod
    def _short_underslide_up_top():
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
            ([(2, 3)], [2, (1, 3)], [(2, 3)]),
            # Combining (2, 3), (1, 2) -> (2, 3) with (0, 1) -> (0, 2)
            # (2, 3), (0, 1)-(1, 2) -> (0, 2)-(2, 3)
            ([(2, 3)], [(0, 1),(1, 2)], [(0, 2),(2, 3)]),
        ]
        return patterns_raw

    @staticmethod
    def _short_underslide_up_middle():
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
            ([(2, 3)], [2, (1, 3)], [(2, 3)]),
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
            ([(2, 3)], [2, (1, 4)], [(2, 4)]),
            # (2, 3)-(3, 4), (1, 3) -> (2, 3)-(3, 4)
            ([(2, 3),(3, 4)], [(1, 3)], [(2, 3),(3, 4)]),
            ([(2, 3),(3, 4)], [2, (1, 3)], [(2, 3),(3, 4)]),
            # (2, 4), (1, 3) -> (2, 4)
            ([1, (2, 4)], [(1, 3)], [1, (2, 4)]),
            ([1, (2, 4)], [2, (1, 3)], [1, (2, 4)]),
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
            ([(2, 4)], [2, (0, 3)], [(0, 1),(2, 4)]),
            ([(0, 1),(2, 4)], [(1, 3)], [(0, 1),(2, 4)]),
            ([(0, 1),(2, 4)], [2, (1, 3)], [(0, 1),(2, 4)]),
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
            ([(2, 4)], [2, (0, 3)], [1, (0, 4)]),
            ([(0, 1),(2, 4)], [(1, 3)], [1, (0, 4)]),
            ([(0, 1),(2, 4)], [2, (1, 3)], [1, (0, 4)]),
            ([(0, 2),(2, 4)], [(2, 3)], [1, (0, 4)]),
        ]
        return patterns_raw
