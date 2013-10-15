"""Producing type DA structures for arcslides, using local actions."""

from dastructure import SimpleDAGenerator, SimpleDAStructure
from dastructure import AddChordToDA
from localpmc import LocalStrandDiagram
from localpmc import restrictPMC, restrictStrandDiagram
from pmc import Strands
from utility import F2

class ArcslideDA:
    """Responsible for producing a type DA structure for an arcslide, using
    local actions.

    """
    @staticmethod
    def getDAStructure(slide):
        """Returns the type DA structure corresponding to slide (type Arcslide
        in arcslide.py).

        """
        dd_idems = slide.getIdems()
        da_idems = [(l_idem, r_idem.opp().comp())
                    for l_idem, r_idem in dd_idems]
        pmc1, pmc2 = slide.start_pmc, slide.end_pmc
        n = pmc1.n
        alg1, alg2 = pmc1.getAlgebra(), pmc2.getAlgebra()
        dastr = SimpleDAStructure(F2, alg1, alg2)
        for i in range(len(da_idems)):
            l_idem, r_idem = da_idems[i]
            dastr.addGenerator(
                SimpleDAGenerator(dastr, l_idem, r_idem, "%d" % i))

        alg1_gens = alg1.getGenerators()
        mod_gens = dastr.getGenerators()
        b1, c1, c2 = slide.b1, slide.c1, slide.c2

        if b1 == 1 and c1 == 2 and c2 == 0:
            local_cut1, outer_cut1, local_cut2, outer_cut2 = (
                [(0, 2)], [(3, n-1)],
                [(0, 2)], [(3, n-1)])
            short_chord = [(0, 1)]
            patterns_raw = ArcslideDA._short_underslide_up_bottom()
        elif b1 == n-2 and c1 == n-1 and c2 == n-3:
            local_cut1, outer_cut1, local_cut2, outer_cut2 = (
                [(n-3, n-1)], [(0, n-4)],
                [(n-3, n-1)], [(0, n-4)])
            short_chord = [(n-3, n-2)]
            patterns_raw = ArcslideDA._short_underslide_up_top()
        elif b1 == c1 - 1 and c2 == c1 - 2:  # *-c2-b1-c1-*
            assert c2 != 0 and c1 != n-1
            local_cut1, outer_cut1, local_cut2, outer_cut2 = (
                [(c2, c1)], [(0, c2-1), (c1+1, n-1)],
                [(c2, c1)], [(0, c2-1), (c1+1, n-1)])
            short_chord = [(c2, b1)]
            patterns_raw = ArcslideDA._short_underslide_up_middle()
        else:
            return NotImplemented

        local_pmc1, mapping1 = restrictPMC(pmc1, local_cut1)
        outer_pmc1, outer_mapping1 = restrictPMC(pmc1, outer_cut1)
        local_pmc2, mapping2 = restrictPMC(pmc2, local_cut2)
        outer_pmc2, outer_mapping2 = restrictPMC(pmc2, outer_cut2)
        # Required so the left to right transition on the outside can proceed.
        assert outer_pmc1 == outer_pmc2

        # Key is tuple of LocalStrandDiagram specifying what the A-side inputs
        # look like in the local PMC. Value is a list of possible local D-side
        # outputs.
        arrow_patterns = {}
        for pattern in patterns_raw:
            key = []
            for i in range(len(pattern)-1):
                key.append(local_pmc1.sd(pattern[i]))
            key = tuple(key)
            if key not in arrow_patterns:
                arrow_patterns[key] = []
            arrow_patterns[key].append(local_pmc2.sd(pattern[-1]))

        # Add action with zero algebra input
        AddChordToDA(dastr, Strands(pmc2, short_chord), [])

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
                    restrictStrandDiagram(pmc1, cur_a, local_pmc1, mapping1),)
                outer_a = restrictStrandDiagram(
                    pmc1, cur_a, outer_pmc1, outer_mapping1)
                # Compute product on the outside
                if cur_prod_d is None:
                    new_prod_d = 1 * outer_a
                else:
                    new_prod_d = cur_prod_d * outer_a
                if new_prod_d == 0:
                    continue
                new_prod_d = new_prod_d.getElt()
                # Local patterns match exactly one of the arrows
                if new_list_local in arrow_patterns:
                    for local_d in arrow_patterns[new_list_local]:
                        alg_d = local_d.join(new_prod_d.removeSingleHor(),
                                             pmc2, mapping2, outer_mapping2)
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
                        for pattern in arrow_patterns
                        if len(pattern) > len(new_list)]):
                    search(new_list, new_list_local, new_prod_d)

        search((), (), None)
        return dastr

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
            ([(1, 3)], [(0, 1)], [(1, 3)]),
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
            # (2, 3), (0, 2) -> (0, 3)
            ([(2, 3)], [(0, 2)], [(0, 3)]),
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
            # (0, 3)-(3, 4) -> (0, 3)-(3, 4) extension of (1, 3) -> (1, 3)
            ([(0, 3),(3, 4)], [(0, 3),(3, 4)]),
            ([2, (0, 3),(3, 4)], [2, (0, 3),(3, 4)]),
            # (0, 1)-(2, 4) -> (0, 1)-(3, 4) extension of (2, 4) -> (3, 4)
            ([(0, 1),(2, 4)], [(0, 1),(3, 4)]),
            # (0, 2)-(2, 4) -> (0, 3)-(3, 4) rather strange
            ([(0, 2),(2, 4)], [(0, 3),(3, 4)]),

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
            # (2, 3), (0, 2) -> (0, 3)
            ([(2, 3)], [(0, 2)], [(0, 3)]),
            # *** Upper extension ***
            # Extension of (2, 3), (1, 2) -> (2, 3)
            # (2, 3)-(3, 4), (1, 2) -> (2, 3)-(3, 4)
            ([(2, 3),(3, 4)], [(1, 2)], [(2, 3),(3, 4)]),
            # (2, 4), (1, 2) -> (2, 4)
            ([(2, 4)], [(1, 2)], [(2, 4)]),
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
            ([(2, 4)], [(0, 1),(1, 2)], [(0, 2),(2, 4)]),
            ([(2, 4)], [(0, 2)], [(0, 4)]),
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
