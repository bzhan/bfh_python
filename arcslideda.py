"""Producing type DA structures for arcslides, using local actions."""

from algebra import E0, TensorGenerator
from dastructure import DAStructure, SimpleDAGenerator, SimpleDAStructure
from dastructure import AddChordToDA
from hdiagram import getArcslideDiagram
from linalg import F2RowSystem
from localpmc import LocalIdempotent, LocalStrandAlgebra, LocalStrandDiagram
from localpmc import restrictPMC, restrictStrandDiagram
from pmc import Strands, StrandDiagram
from utility import memorize, subset
from utility import ACTION_LEFT, ACTION_RIGHT, F2
import itertools
from Queue import Queue

class ArcslideDA(DAStructure):
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
        self.pmc1, self.pmc2 = slide.start_pmc, slide.end_pmc

        # Initiate with usual attributes of a DAStructure
        DAStructure.__init__(self, F2, algebra1 = self.pmc1.getAlgebra(),
                             algebra2 = self.pmc2.getAlgebra(),
                             side1 = ACTION_LEFT, side2 = ACTION_RIGHT)

        n = self.pmc1.n
        b1, c1, c2 = slide.b1, slide.c1, slide.c2
        b1p, c1p, c2p = [slide.to_r[p] for p in b1, c1, c2]

        # Note intervals (start, end) with start > end are ignored.
        if b1 == c1 + 1:  # downward
            if c2 == c1 + 2:  # short underslide downward
                local_cut1, outer_cut1, local_cut2, outer_cut2 = (
                    [(c1, c2)], [(0, c1-1), (c2+1, n-1)],
                    [(c1p, c2p)], [(0, c1p-1), (c2p+1, n-1)])
                patterns_base = ArcslideDA._short_underslide_down_middle()
                if c1 == 0:
                    translator = ([-1, 0, 1, 2, 3], [-1, 0, 1, 2, 3])
                elif c2 == n - 1:
                    translator = ([0, 1, 2, 3, -1], [0, 1, 2, 3, -1])
                else:
                    translator = None
            elif c2 > c1:  # general underslide downward
                local_cut1, outer_cut1, local_cut2, outer_cut2 = (
                    [(c1, b1), (c2, c2)],
                    [(0, c1-1), (b1+1, c2-1), (c2+1, n-1)],
                    [(c1p, c1p), (b1p, c2p)],
                    [(0, c1p-1), (c1p+1, b1p-1), (c2p+1, n-1)])
                patterns_base = ArcslideDA._general_underslide_down_middle()
                if c1 == 0:
                    translator = ([-1, 0, 1, 2, 3, 4, 5],
                                  [-1, 0, 1, 2, 3, 4, 5])
                elif c2 == n - 1:
                    translator = ([0, 1, 2, 3, 4, 5, -1],
                                  [0, 1, 2, 3, 4, 5, -1])
                else:
                    translator = None
            else:  # c2 < c1, general overslide downward
                local_cut1, outer_cut1, local_cut2, outer_cut2 = (
                    [(c2, c2), (c1, b1)],
                    [(0, c2-1), (c2+1, c1-1), (b1+1, n-1)],
                    [(b1p, c2p), (c1p, c1p)],
                    [(0, b1p-1), (c2p+1, c1p-1), (c1p+1, n-1)])
                patterns_base = ArcslideDA._general_underslide_down_middle()
                if c2 == 0 and b1 == n - 1:
                    translator = ([2, 3, 4, -1, -1, 0, 1],
                                  [3, 4, -1, -1, 0, 1, 2])
                elif c2 == 0:
                    translator = ([2, 3, 4, 5, -1, 0, 1],
                                  [3, 4, 5, -1, 0, 1, 2])
                elif b1 == n - 1:
                    translator = ([3, 4, 5, -1, 0, 1, 2],
                                  [4, 5, -1, 0, 1, 2, 3])
                else:
                    translator = ([3, 4, 5, 6, 0, 1, 2], [4, 5, 6, 0, 1, 2, 3])
        elif b1 == c1 - 1:  # upward
            if c2 == c1 - 2:  # short underslide upward
                local_cut1, outer_cut1, local_cut2, outer_cut2 = (
                    [(c2, c1)], [(0, c2-1), (c1+1, n-1)],
                    [(c2p, c1p)], [(0, c2p-1), (c1p+1, n-1)])
                patterns_base = ArcslideDA._short_underslide_up_middle()
                if c2 == 0:
                    translator = ([-1, 0, 1, 2, 3], [-1, 0, 1, 2, 3])
                elif c1 == n - 1:
                    translator = ([0, 1, 2, 3, -1], [0, 1, 2, 3, -1])
                else:
                    translator = None
            elif c2 < c1:  # general underslide upward
                local_cut1, outer_cut1, local_cut2, outer_cut2 = (
                    [(c2, c2), (b1, c1)],
                    [(0, c2-1), (c2+1, b1-1), (c1+1, n-1)],
                    [(c2p, b1p), (c1p, c1p)],
                    [(0, c2p-1), (b1p+1, c1p-1), (c1p+1, n-1)])
                patterns_base = ArcslideDA._general_underslide_up_middle()
                if c2 == 0:
                    translator = ([-1, 0, 1, 2, 3, 4, 5],
                                  [-1, 0, 1, 2, 3, 4, 5])
                elif c1 == n - 1:
                    translator = ([0, 1, 2, 3, 4, 5, -1],
                                  [0, 1, 2, 3, 4, 5, -1])
                else:
                    translator = None
            else:  # c2 > c1, general overslide upward
                local_cut1, outer_cut1, local_cut2, outer_cut2 = (
                    [(b1, c1), (c2, c2)],
                    [(0, b1-1), (c1+1, c2-1), (c2+1, n-1)],
                    [(c1p, c1p), (c2p, b1p)],
                    [(0, c1p-1), (c1p+1, c2p-1), (b1p+1, n-1)])
                patterns_base = ArcslideDA._general_underslide_up_middle()
                if b1 == 0 and c2 == n - 1:
                    translator = ([3, 4, -1, -1, 0, 1, 2],
                                  [2, 3, 4, -1, -1, 0, 1])
                elif b1 == 0:
                    translator = ([3, 4, 5, -1, 0, 1, 2],
                                  [2, 3, 4, 5, -1, 0, 1])
                elif c2 == n - 1:
                    translator = ([4, 5, -1, 0, 1, 2, 3],
                                  [3, 4, 5, -1, 0, 1, 2])
                else:
                    translator = ([4, 5, 6, 0, 1, 2, 3], [3, 4, 5, 6, 0, 1, 2])
        else:
            raise NotImplementedError(
                "This slide pattern is not yet implemented.")

        if translator is None:
            patterns_raw = patterns_base
        else:
            patterns_raw = ArcslideDA._restrict_local_arrows(
                patterns_base, translator[0], translator[1])

        self.local_pmc1, self.mapping1 = restrictPMC(self.pmc1, local_cut1)
        self.outer_pmc1, self.outer_mapping1 = restrictPMC(
            self.pmc1, outer_cut1)
        self.local_pmc2, self.mapping2 = restrictPMC(self.pmc2, local_cut2)
        self.outer_pmc2, self.outer_mapping2 = restrictPMC(
            self.pmc2, outer_cut2)

        # Required so the left to right transition on the outside can proceed.
        assert self.outer_pmc1 == self.outer_pmc2

        # Compute the set of arrow patterns
        self.arrow_patterns = {}
        for pattern in patterns_raw:
            key = []
            for i in range(len(pattern)-1):
                key.append(self.local_pmc2.sd(pattern[i]))
            key = tuple(key)
            if key not in self.arrow_patterns:
                self.arrow_patterns[key] = []
            self.arrow_patterns[key].append(self.local_pmc1.sd(pattern[-1]))

        # Produce the set of possible strict prefixes
        self.strict_prefix_set = set()
        for pattern in self.arrow_patterns:
            for prefix_len in range(len(pattern)):  # excludes full pattern
                self.strict_prefix_set.add(
                    tuple([coeff.removeSingleHor()
                           for coeff in pattern[:prefix_len]]))

        # Compute the set of generators, and add to itself
        dd_idems = self.slide.getIdems()
        da_idems = [(l_idem, r_idem.opp().comp())
                    for l_idem, r_idem in dd_idems]
        self.generators = set()
        for i in range(len(da_idems)):
            l_idem, r_idem = da_idems[i]
            self.generators.add(
                SimpleDAGenerator(self, l_idem, r_idem, "%d" % i))
        # With generators set, add grading. Any generator can serve as base_gen
        for gen in self.generators:
            base_gen = gen
            break
        self.registerHDiagram(getArcslideDiagram(self.slide), base_gen)

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

    @staticmethod
    def adjustSingleHors(coeffs_a):
        """Given a tuple of A-side inputs, delete single idempotents from inputs
        in a best effort to make the idempotents match.

        """
        coeffs_a = list(coeffs_a)
        # Remove single horizontals that cause mismatch between adjacent
        # idempotents.
        while not all([coeffs_a[i].right_idem == coeffs_a[i+1].left_idem
                       for i in range(len(coeffs_a)-1)]):
            changed = False  # Must modify something on each run through
            for i in range(len(coeffs_a)-1):
                if coeffs_a[i].right_idem == coeffs_a[i+1].left_idem:
                    continue
                for idem in coeffs_a[i].right_idem:
                    if idem not in coeffs_a[i+1].left_idem and \
                       idem in coeffs_a[i].single_hor:
                        coeffs_a[i] = coeffs_a[i].removeSingleHor([idem])
                        changed = True
                for idem in coeffs_a[i+1].left_idem:
                    if idem not in coeffs_a[i].right_idem and \
                       idem in coeffs_a[i+1].single_hor:
                        coeffs_a[i+1] = coeffs_a[i+1].removeSingleHor([idem])
                        changed = True
            assert changed
        return tuple(coeffs_a)

    def getDAStructure(self):
        """Returns the simple type DA structure corresponding to slide."""
        dastr = SimpleDAStructure(F2, self.algebra1, self.algebra2)
        for gen in self.getGenerators():
            dastr.addGenerator(SimpleDAGenerator(
                dastr, gen.idem1, gen.idem2, gen.name))

        # Get gradings from Heegaard diagram
        for gen in dastr.getGenerators():
            base_gen = gen
            break
        dastr.registerHDiagram(getArcslideDiagram(self.slide), base_gen)

        alg1_gens = self.algebra1.getGenerators()
        alg2_gens = self.algebra2.getGenerators()
        mod_gens = dastr.getGenerators()

        # Add action with zero algebra input
        short_chord = [tuple(sorted([self.slide.b1, self.slide.c1]))]
        AddChordToDA(dastr, Strands(self.pmc1, short_chord), [])

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
                    self.pmc2, cur_a, self.local_pmc2, self.mapping2),)
                outer_a = restrictStrandDiagram(
                    self.pmc2, cur_a, self.outer_pmc2, self.outer_mapping2)
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
                pattern = ArcslideDA.adjustSingleHors(new_list_local)
                if pattern in self.arrow_patterns:
                    for local_d in self.arrow_patterns[pattern]:
                        alg_d = local_d.join(
                            new_prod_d.removeSingleHor(), self.pmc1,
                            self.mapping1, self.outer_mapping1)
                        if alg_d is None:
                            continue
                        for x, y in itertools.product(mod_gens, mod_gens):
                            if ArcslideDA.idemMatchDA(x, y, alg_d, new_list):
                                dastr.addDelta(x, y, alg_d, new_list, 1)
                # Local patterns match the prefix of one of the arrows
                if tuple([coeff.removeSingleHor() for coeff in new_list_local
                      ]) in self.strict_prefix_set:
                    search(new_list, new_list_local, new_prod_d)

        search((), (), None)
        return dastr

    def getGenerators(self):
        return list(self.generators)

    def delta(self, MGen, algGens):
        # Idempotent must match
        if any([algGens[i].right_idem != algGens[i+1].left_idem
                for i in range(len(algGens)-1)]):
            return E0
        if any([alg.isIdempotent() for alg in algGens]):
            return E0

        result = E0
        mod_gens = self.getGenerators()
        # Take care of case with zero algebra inputs
        if len(algGens) == 0:
            short_chord = [tuple(sorted([self.slide.b1, self.slide.c1]))]
            st = Strands(self.pmc1, short_chord)
            if st.leftCompatible(MGen.idem1):
                alg_d = StrandDiagram(self.algebra1, MGen.idem1, st)
                for y in mod_gens:
                    if ArcslideDA.idemMatchDA(MGen, y, alg_d, []):
                        result += 1 * TensorGenerator((alg_d, y), self.AtensorM)
            return result

        # Take care of remaining case
        if MGen.idem2 != algGens[0].left_idem:
            return E0
        # - Make sure the outside multiplies correctly
        prod_d = restrictStrandDiagram(self.pmc2, algGens[0],
                                       self.outer_pmc2, self.outer_mapping2)
        has_product = True
        for alg in algGens[1:]:
            prod_d = prod_d.parent.multiplyGeneral(
                prod_d, restrictStrandDiagram(
                    self.pmc2, alg, self.outer_pmc2, self.outer_mapping2),
                False)  # strict_idems = False
            if prod_d == 0:
                has_product = False
                break
            else:
                prod_d = prod_d.getElt()
        if not has_product:
            return E0

        # - Now check the inside
        alg_local = tuple([restrictStrandDiagram(
            self.pmc2, alg, self.local_pmc2, self.mapping2) for alg in algGens])
        pattern = ArcslideDA.adjustSingleHors(alg_local)
        if pattern not in self.arrow_patterns:
            return E0
        for local_d in self.arrow_patterns[pattern]:
            alg_d = local_d.join(prod_d.removeSingleHor(), self.pmc1,
                                 self.mapping1, self.outer_mapping1)
            if alg_d is None:
                continue
            for y in mod_gens:
                if ArcslideDA.idemMatchDA(MGen, y, alg_d, algGens):
                    result += 1 * TensorGenerator((alg_d, y), self.AtensorM)
        return result

    def deltaPrefix(self, MGen, algGens):
        if len(algGens) == 0:
            return True
        if MGen.idem2 != algGens[0].left_idem:
            return False
        # Note: Testing shows probably not worth it to multipliy outside.
        # Test inside.
        alg_local = [restrictStrandDiagram(
            self.pmc2, alg, self.local_pmc2, self.mapping2) for alg in algGens]
        alg_local = tuple([alg.removeSingleHor() for alg in alg_local])
        return alg_local in self.strict_prefix_set

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
                for x, y in itertools.product(mod_gens, mod_gens):
                    if ArcslideDA.idemMatchDA(x, y, coeff_d, coeffs_a):
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
            # Starting pair of the short chord
            if b1 > c1:
                start_pair = c_pair1
            else:
                start_pair = b_pair1
            alg_d = LocalStrandDiagram(alg1, [start_pair] + list(idems_to_add),
                                       local_short_chord)
            for x, y in itertools.product(mod_gens, mod_gens):
                if ArcslideDA.idemMatchDA(x, y, alg_d, []):
                    local_dastr.addDelta(x, y, alg_d, [], 1)
        return local_dastr

    @staticmethod
    def _restrict_local_arrows(patterns, point_map_d, point_map_a):
        """Given a list of patterns (in the format of patterns_raw) for more
        general local case, restrict to a more special local case using a
        mapping from points in the general local PMC to a special local PMC.
        point_map_d and point_map_a specifies the point mappings on the D-side
        and the A-side.

        This function operates entirely by translation using the given point
        map, it does not know about formats for local PMCs or local strand
        diagrams.

        """
        def translate(lst, mapping):
            """lst consists of either integers or pairs of integers. Translate
            according to mapping. If any of the translated value is -1, return
            None. Otherwise return the translated list.

            """
            result = []
            for entry in lst:
                if isinstance(entry, int):
                    result.append(mapping[entry])
                    if result[-1] == -1:
                        return None
                else:  # entry must be a pair
                    result.append((mapping[entry[0]], mapping[entry[1]]))
                    if result[-1][0] == -1 or result[-1][1] == -1:
                        return None
            return result

        new_patterns = []
        for pattern in patterns:
            new_pattern = [translate(pattern_d, point_map_a)
                           for pattern_d in pattern[0:-1]]
            new_pattern.append(translate(pattern[-1], point_map_d))
            if all([entry != None for entry in new_pattern]):
                new_patterns.append(new_pattern)

        return new_patterns

    # The next series of functions specify the local arrows. The format is as
    # follows:
    # All but the last element of the tuple is a list to be passed to the sd()
    # function of the A-side local PMC, specifying the A-side inputs. The last
    # element of the tuple is a list to be passed to the sd() function of the
    # D-side local PMC, specifying the D-side output.
    @staticmethod
    def _short_underslide_down_middle():
        """Short underslide going down, in the middle of PMC."""
        # Local PMC is 0*-1-2-3-4*, with 1 and 3 paired.
        patterns_raw = [
            #### Initial patterns
            ([(1, 2)],),
            ([], []), ([1], [1]), ([2], [2]), ([2], [1]), ([1, 2], [1, 2]),
            ([(2, 3)], [1]),
            ([(1, 2)], [(1, 3)]),
            ([(1, 3)], [(1, 3)]),
            ([2, (1, 3)], [2, (1, 3)]),
            ([(1, 2),(2, 3)], [(1, 2),(2, 3)]),
            ([(2, 3)], [(1, 2)], [(2, 3)]),
            ([(2, 3)], [(1, 3)], [(2, 3)]),

            #### Seed for top
            ([(3, 4)], [(3, 4)]),
            ([2, (3, 4)], [2, (3, 4)]),
            # Linear algebra produces
            ([(2, 3), (3, 4)], [(1, 3)], [(2, 3), (3, 4)]),
            ([(2, 3), (3, 4)], [(1, 2)], [(2, 3)], [1, (2, 4)]),
            ([2, (1, 4)], [2, (1, 4)]),
            ([(1, 4)], [(1, 4)]),
            ([1, (2, 4)], [(1, 3)], [1, (2, 4)]),
            ([(1, 2), (2, 4)], [(1, 2), (2, 4)]),
            ([2, (3, 4)], [(2, 3)], [(1, 2)], [1, (2, 4)]),
            ([(2, 4)], [(3, 4)]),
            ([(2, 3), (3, 4)], [(1, 2)], [(1, 2), (2, 4)]),
            ([(2, 3)], [(1, 4)], [(2, 4)]),

            #### Seed for bottom
            ([(0, 1)], [(0, 1)]),
            ([2, (0, 1)], [2, (0, 1)]),
            # Linear algebra produces
            ([1, 2, (0, 4)], [1, 2, (0, 4)]),
            ([1, (0, 4)], [1, (0, 4)]),
            ([(1, 2), (2, 4)], [2, (0, 1)], [1, (0, 2), (2, 4)]),
            ([(2, 3), (3, 4)], [1, (0, 2)], [1, (0, 2), (2, 4)]),
            ([2, (0, 3), (3, 4)], [1, (0, 2), (2, 4)]),
            ([2, (0, 1), (3, 4)], [2, (1, 3)], [1, 2, (0, 4)]),
            ([(0, 1), (3, 4)], [(1, 3)], [1, (0, 4)]),
            ([(1, 2)], [(0, 1), (2, 4)], [1, (0, 4)]),
            ([(2, 4)], [(0, 2)], [(0, 1), (2, 4)]),
            ([2, (0, 4)], [2, (0, 4)]),
            ([(0, 4)], [(0, 4)]),
            ([(0, 1), (1, 2)], [2, (3, 4)], [(0, 1), (1, 4)]),
            ([(2, 3), (3, 4)], [(0, 1), (1, 2)], [(0, 1), (1, 2), (2, 4)]),
            ([1, (2, 4)], [(0, 1), (1, 2)], [1, (0, 2), (2, 4)]),
            ([(2, 3)], [(0, 1), (1, 2)], [2, (3, 4)], [(0, 1), (2, 4)]),
            ([(0, 1), (1, 4)], [(0, 3), (3, 4)]),
            ([2, (0, 1), (1, 4)], [2, (0, 3), (3, 4)]),
            ([2, (1, 3)], [2, (0, 1), (3, 4)], [1, 2, (0, 4)]),
            ([(1, 3)], [(0, 1), (3, 4)], [1, (0, 4)]),
            ([2, (0, 1)], [2, (3, 4)], [(0, 1), (3, 4)]),
            ([(0, 1), (2, 4)], [(1, 3)], [(0, 1), (2, 4)]),
            ([(0, 2), (2, 3), (3, 4)], [(0, 2), (2, 3), (3, 4)]),
            ([(0, 1), (1, 2), (2, 4)], [(0, 2), (2, 3), (3, 4)]),
            ([(0, 2), (2, 4)], [(0, 2), (2, 4)]),
            ([(0, 1), (2, 4)], [(1, 2)], [1, (0, 4)]),
            ([2, (3, 4)], [(2, 3)], [(0, 1), (1, 2)], [1, (0, 2), (2, 4)]),
            ([(2, 3)], [(0, 1), (1, 4)], [1, (0, 4)]),
            ([(2, 4)], [(0, 3)], [1, (0, 4)]),
            ([2, (3, 4)], [2, (0, 3)], [1, 2, (0, 4)]),
            ([(3, 4)], [(0, 3)], [1, (0, 4)]),
            ([(0, 3)], [(0, 3)]),
            ([2, (0, 3)], [2, (0, 3)]),
            ([2, (0, 1)], [(2, 3), (3, 4)], [(0, 1), (3, 4)]),
            ([(1, 4)], [(0, 1)], [1, (0, 4)]),
            ([2, (1, 4)], [2, (0, 1)], [1, 2, (0, 4)]),
            ([2, (0, 1)], [(2, 3), (3, 4)], [(1, 2)], [(0, 1), (2, 4)]),
            ([(2, 3)], [(0, 2), (3, 4)], [(0, 1), (2, 4)]),
            ([(0, 1), (3, 4)], [(1, 2)], [1, (0, 4)]),
            ([(0, 2)], [(0, 3)]),
            ([(0, 1), (3, 4)], [(0, 1), (3, 4)]),
            ([2, (0, 1), (3, 4)], [2, (0, 1), (3, 4)]),
            ([(1, 2), (2, 3)], [2, (0, 1), (3, 4)], [1, (0, 2), (2, 4)]),
            ([(3, 4)], [(0, 2)], [1, (0, 4)]),
            ([(0, 2), (2, 3)], [(0, 2), (2, 3)]),
            ([2, (0, 1)], [1, (0, 2)]),
            ([(2, 3)], [(0, 1), (1, 2)], [(0, 2), (2, 3)]),
            ([(0, 1), (1, 2)], [(2, 3), (3, 4)], [(0, 1), (1, 4)]),
            ([(0, 2), (3, 4)], [(0, 1), (1, 4)]),
            ([(2, 4)], [(0, 2)], [1, (0, 4)]),
            ([(2, 3)], [(0, 3), (3, 4)], [(0, 1), (2, 4)]),
            ([(2, 3)], [(0, 1), (1, 4)], [(0, 1), (2, 4)]),
            ([(0, 1), (1, 2)], [1, (2, 4)], [1, (0, 4)]),
            ([(0, 3), (3, 4)], [(0, 1), (1, 4)]),
            ([2, (0, 3), (3, 4)], [2, (0, 1), (1, 4)]),
            ([(0, 1), (2, 4)], [(1, 2)], [(0, 1), (2, 4)]),
            ([2, (0, 1)], [(2, 3), (3, 4)], [(1, 2)], [1, (0, 4)]),
            ([(0, 1), (2, 4)], [(1, 3)], [1, (0, 4)]),
            ([(2, 3)], [(0, 1), (1, 2)], [(2, 3), (3, 4)], [(0, 1), (2, 4)]),
            ([(2, 4)], [(0, 3)], [(0, 1), (2, 4)]),
            ([(0, 1), (2, 4)], [(0, 1), (3, 4)]),
        ]
        return patterns_raw

    @staticmethod
    def _general_underslide_down_middle():
        """Underslide of length >= 3 going down, in the middle of PMC."""
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
            ([(3, 4), (4, 6)], [(4, 5), (5, 6)]),
            ([1, (3, 4)], [(4, 5), (5, 6)], [1, (4, 6)]),
            ([(3, 4), (5, 6)], [4, (0, 1)], [1, (0, 2), (4, 6)]),
            ([(0, 1), (3, 4)], [(0, 2), (4, 5)]),
            ([(4, 6)], [(0, 1), (3, 4)], [1, (0, 2), (4, 6)]),
            ([(3, 4), (4, 6)], [4, (0, 1)], [1, (0, 2), (4, 6)]),
            ([1, (3, 4)], [4, (5, 6)], [1, (4, 6)]),
            ([(0, 1), (5, 6)], [1, (3, 4)], [1, (0, 2), (4, 6)]),
            ([(4, 6)], [(3, 5)], [1, (4, 6)]),
            ([(0, 1), (3, 4)], [(4, 5), (5, 6)], [(0, 1), (4, 6)]),
            ([(3, 4), (5, 6)], [(4, 5), (5, 6)]),
            ([4, (0, 1), (3, 6)], [2, (0, 1), (4, 6)]),
            ([(0, 1), (3, 6)], [(0, 1), (4, 6)]),
            ([4, (3, 6)], [1, (4, 6)]),
            ([(3, 4), (4, 5)], [4, (5, 6)], [1, (4, 6)]),
            ([4, (0, 1)], [(3, 4), (4, 5), (5, 6)], [(0, 2), (4, 5), (5, 6)]),
            ([(0, 1), (3, 4), (5, 6)], [(0, 1), (1, 2), (4, 6)]),
            ([(0, 1), (4, 6)], [1, (3, 4)], [1, (0, 2), (4, 6)]),
            ([1, (3, 4)], [(1, 2), (4, 5)]),
            ([(0, 1), (3, 4)], [4, (5, 6)], [(0, 1), (4, 6)]),
            ([1, 4, (3, 6)], [1, 2, (4, 6)]),
            ([1, (3, 6)], [1, (4, 6)]),
            ([(3, 5), (5, 6)], [(4, 5), (5, 6)]),
            ([4, (3, 5), (5, 6)], [2, (4, 5), (5, 6)]),
            ([(3, 6)], [(4, 6)]),
            ([4, (3, 6)], [2, (4, 6)]),
            ([4, (0, 1)], [1, (3, 4), (4, 6)], [1, (0, 2), (4, 6)]),
            ([(5, 6)], [(0, 1), (3, 4)], [1, (0, 2), (4, 6)]),
            ([(3, 4), (4, 5)], [(1, 2), (4, 5)]),
            ([(0, 1), (3, 4), (4, 6)], [(0, 1), (1, 2), (4, 6)]),
            ([(4, 5)], [(3, 5), (5, 6)], [1, (4, 6)]),
            ([(3, 4), (4, 6)], [(4, 5)], [1, (4, 6)]),
            ([4, (3, 5)], [2, (4, 5)]),
            ([(3, 5)], [(4, 5)]),

            #### Seed for lower middle
            ([4, (1, 2)], [1, (2, 3)]),
            # Linear algebra produces
            ([4, (1, 2), (3, 6)], [2, (1, 3), (4, 6)]),
            ([(1, 2), (3, 6)], [(1, 3), (4, 6)]),
            ([4, (0, 1), (1, 2)], [4, (5, 6)], [2, (0, 3), (5, 6)]),
            ([(0, 1), (1, 2)], [(5, 6)], [(0, 3), (5, 6)]),
            ([1, (3, 4)], [(4, 5), (5, 6)], [(0, 1), (1, 2)], [1, (0, 3), (4, 6)]),
            ([(0, 1), (1, 2), (4, 6)], [(0, 1), (2, 3), (5, 6)]),
            ([(4, 5)], [(1, 2), (3, 6)], [(2, 3), (4, 6)]),
            ([4, (0, 1), (1, 2)], [4, (5, 6)], [(3, 4), (4, 5)], [1, (0, 2), (2, 3), (4, 6)]),
            ([(4, 5)], [(1, 2), (3, 4)], [(2, 3), (4, 5)]),
            ([(3, 4), (4, 6)], [(0, 2), (4, 5)], [(0, 1), (2, 3), (4, 6)]),
            ([4, (0, 2), (3, 5), (5, 6)], [2, (0, 3), (4, 5), (5, 6)]),
            ([(0, 2), (3, 5), (5, 6)], [(0, 3), (4, 5), (5, 6)]),
            ([(1, 2), (3, 4)], [(0, 1), (4, 6)], [1, (0, 3), (4, 6)]),
            ([4, (0, 2), (3, 5)], [2, (0, 3), (4, 5)]),
            ([(0, 2), (3, 5)], [(0, 3), (4, 5)]),
            ([(1, 2), (3, 4), (4, 6)], [4, (0, 1)], [1, (0, 2), (2, 3), (4, 6)]),
            ([1, (0, 2), (3, 4), (4, 6)], [1, (0, 2), (2, 3), (4, 6)]),
            ([(1, 2), (4, 6)], [(3, 5)], [1, (2, 3), (4, 6)]),
            ([4, (0, 1), (1, 2), (3, 6)], [1, (0, 2), (2, 3), (4, 6)]),
            ([(0, 2), (3, 4), (4, 5), (5, 6)], [(0, 2), (2, 3), (4, 5), (5, 6)]),
            ([(0, 1), (1, 2)], [1, (3, 4)], [4, (5, 6)], [1, (0, 3), (4, 6)]),
            ([(4, 5)], [(1, 2)], [(2, 3)]),
            ([1, (0, 2), (3, 4)], [(4, 5), (5, 6)], [1, (0, 3), (4, 6)]),
            ([(0, 2), (3, 4), (4, 5)], [(0, 2), (2, 3), (4, 5)]),
            ([(4, 5), (5, 6)], [(0, 1), (1, 2)], [1, (3, 4)], [1, (0, 2), (2, 3), (4, 6)]),
            ([1, (3, 4)], [(0, 2), (4, 5), (5, 6)], [1, (0, 3), (4, 6)]),
            ([(3, 4), (4, 5), (5, 6)], [4, (1, 2)], [(2, 3), (4, 5), (5, 6)]),
            ([4, (1, 2)], [4, (0, 1)], [(3, 4), (4, 5), (5, 6)], [1, (0, 2), (2, 3), (4, 6)]),
            ([(4, 5)], [(0, 1), (1, 2), (3, 6)], [1, (0, 3), (4, 6)]),
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
            ([(0, 1), (1, 2), (3, 6)], [(0, 3), (4, 5), (5, 6)]),
            ([4, (0, 1), (1, 2), (3, 6)], [2, (0, 3), (4, 5), (5, 6)]),
            ([4, (0, 2), (5, 6)], [2, (0, 3), (5, 6)]),
            ([(0, 2), (5, 6)], [(0, 3), (5, 6)]),
            ([(0, 1), (1, 2), (3, 4)], [4, (5, 6)], [(0, 1), (1, 3), (4, 6)]),
            ([(1, 2), (3, 4), (4, 5)], [(1, 2), (2, 3), (4, 5)]),
            ([1, (3, 4)], [(0, 1), (1, 2), (4, 6)], [1, (0, 3), (4, 6)]),
            ([(3, 4), (4, 5)], [(4, 5), (5, 6)], [(0, 1), (1, 2)], [(0, 1), (2, 3), (4, 6)]),
            ([4, (1, 2)], [2, (1, 3)]),
            ([(1, 2)], [(1, 3)]),
            ([(4, 5)], [(0, 2), (3, 5), (5, 6)], [1, (0, 3), (4, 6)]),
            ([(4, 6)], [(3, 5)], [(0, 1), (1, 2)], [(0, 1), (2, 3), (4, 6)]),
            ([(3, 4), (4, 5), (5, 6)], [4, (0, 1), (1, 2)], [(0, 2), (2, 3), (4, 5), (5, 6)]),
            ([(3, 4), (4, 5)], [4, (0, 1), (1, 2)], [(4, 5), (5, 6)], [(0, 1), (2, 3), (4, 6)]),
            ([4, (0, 2)], [2, (0, 3)]),
            ([(0, 2)], [(0, 3)]),
            ([(0, 2), (3, 4), (4, 6)], [(0, 3), (4, 5), (5, 6)]),
            ([(1, 2), (4, 5)], [1, (2, 3)]),
            ([(1, 2), (3, 4), (4, 5)], [4, (5, 6)], [1, (2, 3), (4, 6)]),
            ([(1, 2), (3, 4), (4, 6)], [(2, 3), (4, 5), (5, 6)]),
            ([(1, 2), (3, 5)], [(0, 1), (5, 6)], [1, (0, 3), (4, 6)]),
            ([4, (1, 2), (3, 5)], [4, (0, 1), (5, 6)], [1, 2, (0, 3), (4, 6)]),
            ([(1, 2), (4, 6)], [(2, 3), (5, 6)]),
            ([(3, 4), (4, 5)], [(0, 2), (4, 5), (5, 6)], [(0, 1), (2, 3), (4, 6)]),
            ([(4, 6)], [(0, 2), (3, 5)], [(0, 1), (2, 3), (4, 6)]),
            ([(0, 2), (4, 5)], [1, (0, 3)]),
            ([(4, 5)], [(0, 1), (1, 2)], [1, (0, 3)]),
            ([(0, 2), (4, 5)], [(3, 5), (5, 6)], [(0, 1), (2, 3), (4, 6)]),
            ([1, (3, 4), (4, 6)], [(1, 2), (4, 5)], [1, (2, 3), (4, 6)]),
            ([(1, 2), (3, 4)], [(1, 3), (4, 5)]),
            ([(4, 5), (5, 6)], [(1, 2), (3, 4)], [4, (0, 1)], [1, (0, 2), (2, 3), (4, 6)]),
            ([(4, 6)], [(3, 5)], [(0, 1), (1, 2)], [1, (0, 3), (4, 6)]),
            ([(3, 4), (4, 5)], [4, (0, 2), (5, 6)], [1, (0, 3), (4, 6)]),
            ([(0, 1), (3, 4), (5, 6)], [4, (1, 2)], [1, (0, 3), (4, 6)]),
            ([(4, 5)], [(3, 5), (5, 6)], [(0, 1), (1, 2)], [1, (0, 3), (4, 6)]),
            ([(3, 4), (4, 5)], [(0, 1), (1, 2), (4, 6)], [(0, 1), (2, 3), (4, 6)]),
            ([(0, 2), (4, 5)], [(0, 1), (2, 3)]),
            ([4, (0, 1), (1, 2)], [4, (5, 6)], [(0, 1), (2, 3), (5, 6)]),
            ([(3, 4), (4, 5)], [(1, 2), (4, 6)], [(2, 3), (4, 6)]),
            ([(1, 2), (4, 6)], [(0, 1), (3, 4)], [1, (0, 2), (2, 3), (4, 6)]),
            ([4, (0, 1), (1, 2)], [4, (3, 5), (5, 6)], [1, (0, 2), (2, 3), (4, 6)]),
            ([(4, 5)], [(0, 1), (1, 2), (3, 4)], [(4, 5), (5, 6)], [(0, 1), (2, 3), (4, 6)]),
            ([1, 4, (0, 2)], [(3, 4), (4, 5), (5, 6)], [1, (0, 2), (2, 3), (4, 6)]),
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
            ([(4, 5), (5, 6)], [1, (0, 2), (3, 4)], [1, (0, 2), (2, 3), (4, 6)]),
            ([(0, 1), (4, 6)], [(1, 2), (3, 5)], [(0, 1), (2, 3), (4, 6)]),
            ([(3, 4), (4, 5)], [4, (1, 2)], [(2, 3), (4, 5)]),
            ([(4, 5)], [(3, 5), (5, 6)], [(0, 1), (1, 2)], [(0, 1), (2, 3), (4, 6)]),
            ([(3, 4), (4, 6)], [(4, 5)], [(1, 2)], [(2, 3), (4, 6)]),
            ([(3, 4), (4, 6)], [(4, 5)], [(0, 1), (1, 2)], [1, (0, 3), (4, 6)]),
            ([(0, 2), (4, 6)], [(0, 3), (5, 6)]),
            ([4, (5, 6)], [(4, 5)], [(0, 1), (1, 2), (3, 4)], [1, (0, 2), (2, 3), (4, 6)]),
            ([(0, 2), (3, 4), (4, 6)], [(4, 5)], [1, (0, 3), (4, 6)]),
            ([4, (3, 5)], [4, (0, 1), (1, 2)], [2, (0, 3), (4, 5)]),
            ([(3, 5)], [(0, 1), (1, 2)], [(0, 3), (4, 5)]),
            ([4, (1, 2), (3, 6)], [4, (0, 1)], [1, 2, (0, 3), (4, 6)]),
            ([(1, 2), (3, 6)], [(0, 1)], [1, (0, 3), (4, 6)]),
            ([(3, 4), (4, 5), (5, 6)], [(1, 2), (4, 5)], [(2, 3), (4, 5), (5, 6)]),
            ([(0, 1), (1, 2), (3, 4)], [1, (4, 6)], [1, (0, 3), (4, 6)]),
            ([(0, 1), (4, 6)], [(1, 2), (3, 5)], [1, (0, 3), (4, 6)]),
            ([4, (1, 2)], [4, (0, 1), (3, 6)], [1, 2, (0, 3), (4, 6)]),
            ([(1, 2)], [(0, 1), (3, 6)], [1, (0, 3), (4, 6)]),
            ([4, (0, 1), (1, 2)], [(4, 5), (5, 6)], [1, (3, 4)], [1, (0, 2), (2, 3), (4, 6)]),
            ([(4, 5), (5, 6)], [(1, 2), (3, 4)], [(2, 3), (4, 5), (5, 6)]),
            ([(3, 4), (4, 5)], [4, (0, 1), (1, 2)], [4, (5, 6)], [(0, 1), (2, 3), (4, 6)]),
            ([(0, 1), (1, 2), (3, 4)], [(4, 5), (5, 6)], [(0, 1), (1, 3), (4, 6)]),
            ([(4, 5), (5, 6)], [1, (3, 4)], [4, (0, 1), (1, 2)], [1, (0, 2), (2, 3), (4, 6)]),
            ([(3, 4), (4, 5)], [(4, 5), (5, 6)], [(1, 2)], [(2, 3), (4, 6)]),
            ([(0, 1), (3, 4), (4, 6)], [(1, 2), (4, 5)], [(0, 1), (2, 3), (4, 6)]),
            ([(0, 2), (3, 4), (4, 6)], [(4, 5)], [(0, 1), (2, 3), (4, 6)]),
            ([4, (1, 2), (3, 5)], [2, (1, 3), (4, 5)]),
            ([(1, 2), (3, 5)], [(1, 3), (4, 5)]),
            ([(3, 4), (4, 6)], [(4, 5)], [(0, 1), (1, 2)], [(0, 1), (2, 3), (4, 6)]),
            ([(1, 2), (4, 5)], [(3, 5), (5, 6)], [1, (2, 3), (4, 6)]),
            ([(0, 2), (3, 4), (4, 5)], [4, (5, 6)], [(0, 1), (2, 3), (4, 6)]),
            ([(0, 1), (3, 4), (5, 6)], [(1, 2), (4, 5)], [1, (0, 3), (4, 6)]),
            ([4, (1, 2), (3, 6)], [1, (2, 3), (4, 6)]),
            ([4, (0, 1), (1, 2)], [(4, 5), (5, 6)], [(0, 1), (2, 3), (5, 6)]),
            ([1, (3, 4)], [4, (0, 2), (5, 6)], [1, (0, 3), (4, 6)]),
            ([(0, 2), (4, 6)], [(3, 5)], [1, (0, 3), (4, 6)]),
            ([(3, 4), (4, 5)], [4, (0, 1), (1, 2)], [(0, 2), (2, 3), (4, 5)]),
            ([1, (3, 4), (4, 6)], [4, (0, 1), (1, 2)], [1, (0, 2), (2, 3), (4, 6)]),
            ([(0, 1), (1, 2), (3, 4), (4, 6)], [(0, 1), (1, 2), (2, 3), (4, 6)]),
            ([4, (0, 2), (5, 6)], [(0, 1), (2, 3), (5, 6)]),
            ([(0, 1), (1, 2), (4, 6)], [1, (3, 4)], [1, (0, 2), (2, 3), (4, 6)]),
            ([(0, 1), (1, 2)], [(3, 5), (5, 6)], [(0, 3), (4, 5), (5, 6)]),
            ([4, (0, 1), (1, 2)], [4, (3, 5), (5, 6)], [2, (0, 3), (4, 5), (5, 6)]),
            ([1, (4, 6)], [(1, 2), (3, 5)], [1, (2, 3), (4, 6)]),
            ([(0, 1), (1, 2)], [(3, 4), (5, 6)], [(0, 3), (4, 5), (5, 6)]),
            ([(3, 4), (4, 5)], [(1, 2), (4, 5)], [(2, 3), (4, 5)]),
            ([(4, 6)], [(3, 5)], [(1, 2)], [(2, 3), (4, 6)]),
            ([4, (0, 2), (3, 6)], [1, (0, 3), (4, 6)]),
            ([4, (0, 2)], [(0, 1), (2, 3)]),
            ([(1, 2), (3, 4), (4, 6)], [(4, 5)], [1, (2, 3), (4, 6)]),
            ([(0, 2), (4, 5), (5, 6)], [(0, 1), (2, 3), (5, 6)]),
            ([(3, 4), (4, 6)], [(0, 2), (4, 5)], [1, (0, 3), (4, 6)]),
            ([(0, 2), (3, 4), (4, 5)], [4, (5, 6)], [1, (0, 3), (4, 6)]),
            ([(4, 5)], [(3, 5), (5, 6)], [(1, 2)], [(2, 3), (4, 6)]),
            ([(4, 5)], [(0, 1), (1, 2), (3, 4)], [4, (5, 6)], [(0, 1), (2, 3), (4, 6)]),
            ([(4, 5), (5, 6)], [(1, 2)], [(0, 1), (3, 4)], [1, (0, 2), (2, 3), (4, 6)]),
            ([(0, 2), (4, 6)], [(3, 5)], [(0, 1), (2, 3), (4, 6)]),
            ([(4, 5)], [(0, 1), (1, 2), (3, 6)], [(0, 1), (2, 3), (4, 6)]),
            ([(4, 5)], [(0, 1), (1, 2)], [(0, 1), (2, 3)]),
        ]
        return patterns_raw

    @staticmethod
    def _short_underslide_up_middle():
        """Short underslide going up, in the middle of PMC."""
        # Local PMC is 0*-1-2-3-4*, with 1 and 3 paired.
        patterns_raw = [
            #### Initial patterns
            ([(2, 3)],),
            ([], []), ([1], [1]), ([2], [2]), ([2], [1]), ([1, 2], [1, 2]),
            ([(1, 2)], [1]),
            ([(2, 3)], [(1, 3)]),
            ([(1, 3)], [(1, 3)]),
            ([2, (1, 3)], [2, (1, 3)]),
            ([(1, 2),(2, 3)], [(1, 2),(2, 3)]),
            ([(2, 3)], [(1, 2)], [(1, 2)]),
            ([(1, 3)], [(1, 2)], [(1, 2)]),

            #### Seed for top
            ([(3, 4)], [(3, 4)]),
            ([2, (3, 4)], [2, (3, 4)]),
            # Linear algebra produces
            ([(1, 4)], [(1, 4)]),
            ([2, (1, 4)], [2, (1, 4)]),
            ([(2, 4)], [(1, 4)]),
            ([2, (3, 4)], [1, (2, 4)]),
            ([(2, 3), (3, 4)], [(1, 2)], [(1, 2), (2, 4)]),
            ([(1, 2), (2, 4)], [(1, 2), (2, 4)]),

            #### Seed for bottom
            ([(0, 1)], [(0, 1)]),
            ([2, (0, 1)], [2, (0, 1)]),
            # Linear algebra produces
            ([(1, 2), (2, 4)], [2, (0, 1)], [1, (0, 2), (2, 4)]),
            ([2, (0, 1)], [1, (2, 4)], [1, (0, 4)]),
            ([2, (0, 1)], [(2, 3), (3, 4)], [(0, 3), (3, 4)]),
            ([(0, 1), (1, 2)], [2, (3, 4)], [(0, 1), (3, 4)]),
            ([(1, 2)], [2, (0, 1)], [(2, 3), (3, 4)], [1, (0, 4)]),
            ([(1, 3)], [(0, 1), (1, 2)], [2, (3, 4)], [1, (0, 4)]),
            ([(1, 2)], [(0, 1), (2, 4)], [1, (0, 4)]),
            ([(2, 3)], [(0, 1), (1, 2)], [(0, 1), (1, 2)]),
            ([(0, 1), (1, 2)], [2, (3, 4)], [(2, 3)], [1, (0, 4)]),
            ([(0, 4)], [(0, 4)]),
            ([2, (0, 4)], [2, (0, 4)]),
            ([(0, 1), (1, 2)], [1, (2, 4)], [1, (0, 4)]),
            ([(0, 2), (2, 3)], [2, (3, 4)], [(0, 2), (3, 4)]),
            ([(0, 1), (1, 4)], [(0, 3), (3, 4)]),
            ([2, (0, 1), (1, 4)], [2, (0, 3), (3, 4)]),
            ([2, (1, 3)], [2, (0, 1), (3, 4)], [1, 2, (0, 4)]),
            ([(1, 3)], [(0, 1), (3, 4)], [1, (0, 4)]),
            ([2, (0, 1)], [2, (3, 4)], [(0, 1), (3, 4)]),
            ([(0, 3), (3, 4)], [(1, 2)], [(0, 2), (3, 4)]),
            ([(2, 3)], [(0, 2), (3, 4)], [1, (0, 4)]),
            ([(0, 2), (2, 3)], [(0, 2), (2, 3)]),
            ([(0, 3), (3, 4)], [(0, 3), (3, 4)]),
            ([2, (0, 3), (3, 4)], [2, (0, 3), (3, 4)]),
            ([(0, 1), (2, 4)], [(0, 1), (1, 4)]),
            ([(0, 1), (1, 2), (2, 4)], [(0, 2), (2, 3), (3, 4)]),
            ([(0, 1), (2, 4)], [(1, 2)], [1, (0, 4)]),
            ([(0, 2), (3, 4)], [(0, 1), (3, 4)]),
            ([2, (3, 4)], [(2, 3)], [(0, 1), (1, 2)], [1, (0, 2), (2, 4)]),
            ([(2, 4)], [(0, 2)], [1, (0, 4)]),
            ([(1, 4)], [(0, 2)], [(0, 2), (3, 4)]),
            ([(0, 3)], [(0, 3)]),
            ([2, (0, 3)], [2, (0, 3)]),
            ([(2, 3), (3, 4)], [(1, 2)], [2, (0, 1)], [1, (0, 2), (2, 4)]),
            ([(1, 2)], [(0, 2), (2, 4)], [(0, 2), (3, 4)]),
            ([2, (0, 4)], [(0, 2), (3, 4)]),
            ([(0, 1), (1, 4)], [(1, 2)], [(0, 2), (3, 4)]),
            ([(1, 3)], [(0, 1), (1, 2)], [2, (3, 4)], [(0, 2), (3, 4)]),
            ([(1, 4)], [(0, 1)], [1, (0, 4)]),
            ([2, (1, 4)], [2, (0, 1)], [1, 2, (0, 4)]),
            ([(0, 1), (3, 4)], [(0, 1), (3, 4)]),
            ([2, (0, 1), (3, 4)], [2, (0, 1), (3, 4)]),
            ([(1, 2), (2, 3)], [2, (0, 1), (3, 4)], [1, (0, 2), (2, 4)]),
            ([(0, 2), (2, 3), (3, 4)], [(0, 2), (2, 3), (3, 4)]),
            ([(0, 1), (1, 2)], [(2, 3), (3, 4)], [(1, 2)], [(0, 2), (3, 4)]),
            ([(0, 2)], [(0, 1)]),
            ([(2, 3), (3, 4)], [(0, 1), (1, 2)], [(0, 2), (2, 3), (3, 4)]),
            ([(0, 1), (1, 2)], [(2, 3), (3, 4)], [(0, 1), (1, 4)]),
            ([1, (0, 2)], [(2, 3), (3, 4)], [1, (0, 4)]),
            ([(1, 3)], [(0, 2), (3, 4)], [(0, 2), (3, 4)]),
            ([(0, 3)], [(1, 2)], [(0, 2)]),
            ([2, (0, 1)], [(1, 2), (2, 4)], [(0, 2), (3, 4)]),
            ([(1, 4)], [(0, 2)], [1, (0, 4)]),
            ([(0, 1), (1, 2)], [(2, 3), (3, 4)], [(1, 2)], [1, (0, 4)]),
            ([(1, 3)], [(1, 2)], [2, (0, 1)], [1, (0, 2)]),
            ([(1, 3)], [(0, 2), (3, 4)], [1, (0, 4)]),
            ([(1, 2)], [(2, 3)], [(0, 1), (1, 2)], [1, (0, 2)]),
            ([(0, 2), (2, 4)], [(0, 3), (3, 4)]),
            ([(2, 3)], [1, (0, 2)], [1, (0, 2)]),
            ([2, (0, 1)], [(2, 3), (3, 4)], [(1, 2)], [(0, 2), (3, 4)]),
        ]
        return patterns_raw

    @staticmethod
    def _general_underslide_up_middle():
        """Underslide of length >= 3 going up, in the middle of PMC."""
        # Local PMC at left (D-side) is 0*-1-2*, 3*-4-5-6*, with 1 and 5 paired.
        # Local PMC at right (A-side) is 0*-1-2-3*, 4*-5-6*, with 1 and 5
        # paired.
        patterns_raw = [
            #### Initial patterns
            ([(4, 5)],),
            ([], []), ([1], [1]), ([2], [4]), ([2], [1]), ([1, 2], [1, 4]),
            ([(1, 2)], [1]),

            #### Seed for top
            ([(5, 6)], [(5, 6)]),
            ([2, (5, 6)], [4, (5, 6)]),
            # Linear algebra produces
            ([2, (5, 6)], [1, (4, 6)]),

            #### Seed for bottom
            ([(0, 1)], [(0, 1)]),
            ([2, (0, 1)], [4, (0, 1)]),
            # Linear algebra produces
            ([2, (0, 1), (5, 6)], [4, (0, 1), (5, 6)]),
            ([(0, 1), (5, 6)], [(0, 1), (5, 6)]),
            ([2, (0, 1)], [2, (5, 6)], [(0, 1), (5, 6)]),
            ([(0, 2)], [(0, 1)]),
            ([(0, 1), (1, 2)], [2, (5, 6)], [(0, 1), (5, 6)]),
            ([(0, 2), (5, 6)], [(0, 1), (5, 6)]),

            #### Seed for upper middle
            ([2, (4, 5)], [1, (3, 4)]),
            # Linear algebra produces
            ([(1, 2), (4, 5)], [1, (3, 4)]),
            ([1, 2, (4, 6)], [1, 4, (3, 6)]),
            ([1, (4, 6)], [1, (3, 6)]),
            ([2, (0, 1), (4, 6)], [(0, 1), (3, 4), (5, 6)]),
            ([(4, 5)], [(0, 1), (1, 2)], [(0, 1), (3, 4)]),
            ([2, (4, 5), (5, 6)], [1, (3, 4), (4, 6)]),
            ([(4, 5)], [(1, 2)], [(3, 4)]),
            ([2, (4, 6)], [1, (3, 6)]),
            ([2, (4, 5), (5, 6)], [4, (3, 5), (5, 6)]),
            ([(4, 5), (5, 6)], [(3, 5), (5, 6)]),
            ([(4, 6)], [(3, 6)]),
            ([2, (4, 6)], [4, (3, 6)]),
            ([2, (4, 6)], [(3, 4), (5, 6)]),
            ([(0, 1), (1, 2), (4, 6)], [(0, 1), (3, 4), (5, 6)]),
            ([(0, 2), (4, 6)], [(0, 1), (3, 6)]),
            ([(4, 5)], [(3, 5)]),
            ([2, (4, 5)], [4, (3, 5)]),
            ([(4, 5), (5, 6)], [(1, 2)], [(3, 4), (5, 6)]),
            ([2, (0, 1), (4, 6)], [4, (0, 1), (3, 6)]),
            ([(0, 1), (4, 6)], [(0, 1), (3, 6)]),
            ([(1, 2), (4, 6)], [(3, 4), (5, 6)]),
            ([(4, 5), (5, 6)], [(0, 1), (1, 2)], [(0, 1), (3, 4), (5, 6)]),
            ([(0, 1), (1, 2)], [2, (4, 5), (5, 6)], [(0, 1), (3, 4), (5, 6)]),
            ([(0, 2), (4, 5), (5, 6)], [(0, 1), (3, 4), (5, 6)]),
            ([2, (0, 1)], [2, (4, 5), (5, 6)], [(0, 1), (3, 4), (5, 6)]),
            ([(0, 2), (4, 5)], [(0, 1), (3, 4)]),
            ([(1, 2), (4, 6)], [1, (3, 6)]),

            # Seed for lower middle
            ([(2, 3)], [(1, 2)]),
            # Linear algebra produces
            ([2, (0, 3)], [4, (0, 2)]),
            ([(0, 3)], [(0, 2)]),
            ([(4, 5), (5, 6)], [(0, 1), (1, 2)], [(1, 2), (2, 3)], [(0, 2), (3, 4), (5, 6)]),
            ([(0, 1), (1, 2)], [(2, 3), (5, 6)], [(0, 2), (5, 6)]),
            ([2, (5, 6)], [2, (4, 5)], [(0, 1), (1, 2), (2, 3)], [1, (0, 2), (3, 4), (4, 6)]),
            ([(0, 3), (4, 6)], [(0, 2), (3, 6)]),
            ([2, (0, 3), (4, 6)], [4, (0, 2), (3, 6)]),
            ([(0, 1), (1, 2), (4, 6)], [(1, 2), (2, 3)], [(0, 2), (3, 4), (5, 6)]),
            ([(2, 3)], [(0, 1), (4, 6)], [1, (0, 2), (3, 6)]),
            ([(2, 3), (4, 5)], [(1, 2), (3, 5)]),
            ([(0, 1), (1, 3)], [(1, 2), (4, 5)], [1, (0, 2), (3, 4)]),
            ([2, (0, 1), (1, 3)], [2, (5, 6)], [1, (0, 2), (4, 6)]),
            ([2, (0, 1)], [(1, 2), (2, 3), (4, 6)], [(0, 2), (3, 4), (5, 6)]),
            ([(0, 1), (2, 3)], [(4, 5), (5, 6)], [(1, 2)], [(0, 2), (3, 4), (5, 6)]),
            ([(2, 3), (4, 5)], [(0, 2), (5, 6)], [1, (0, 2), (3, 6)]),
            ([(0, 2), (4, 5)], [(1, 2), (2, 3)], [(0, 2), (3, 4)]),
            ([(2, 3), (4, 6)], [(0, 2)], [1, (0, 2), (3, 6)]),
            ([(0, 1), (1, 2), (2, 3), (4, 6)], [(0, 1), (1, 2), (3, 4), (4, 6)]),
            ([1, (4, 6)], [(0, 1), (1, 3)], [1, (0, 2), (3, 6)]),
            ([1, 2, (4, 6)], [2, (0, 1), (1, 3)], [1, 4, (0, 2), (3, 6)]),
            ([1, (2, 3)], [(0, 2), (5, 6)], [1, (0, 2), (4, 6)]),
            ([(1, 3), (4, 5)], [(1, 2), (3, 5)]),
            ([2, (1, 3), (4, 5)], [4, (1, 2), (3, 5)]),
            ([(0, 3), (4, 5), (5, 6)], [(0, 2), (3, 5), (5, 6)]),
            ([2, (0, 3), (4, 5), (5, 6)], [4, (0, 2), (3, 5), (5, 6)]),
            ([2, (0, 3), (4, 6)], [(0, 2), (3, 4), (5, 6)]),
            ([2, (0, 1), (1, 3), (4, 6)], [4, (0, 2), (3, 5), (5, 6)]),
            ([(0, 1), (1, 3), (4, 6)], [(0, 2), (3, 5), (5, 6)]),
            ([1, (2, 3)], [(0, 1), (1, 2)], [2, (5, 6)], [1, (0, 2), (4, 6)]),
            ([2, (0, 1)], [1, (2, 3)], [(4, 5), (5, 6)], [1, (0, 2), (3, 6)]),
            ([2, (0, 3), (4, 5)], [4, (0, 2), (3, 5)]),
            ([(0, 3), (4, 5)], [(0, 2), (3, 5)]),
            ([(2, 3)], [(0, 1)], [(4, 5), (5, 6)], [1, (0, 2), (3, 6)]),
            ([(0, 1), (2, 3), (5, 6)], [(0, 2), (4, 5), (5, 6)]),
            ([2, (0, 1)], [(2, 3), (4, 5), (5, 6)], [(0, 2), (3, 5), (5, 6)]),
            ([2, (0, 3), (5, 6)], [4, (0, 2), (5, 6)]),
            ([(0, 3), (5, 6)], [(0, 2), (5, 6)]),
            ([(0, 1), (2, 3)], [1, (4, 6)], [1, (0, 2), (3, 6)]),
            ([(1, 2)], [(2, 3), (4, 5)], [(0, 1), (1, 2)], [1, (0, 2), (3, 4)]),
            ([(0, 1), (1, 3), (4, 6)], [(1, 2)], [(0, 2), (3, 4), (5, 6)]),
            ([(0, 1), (1, 2), (2, 3)], [2, (4, 5), (5, 6)], [(0, 2), (3, 4), (4, 5), (5, 6)]),
            ([(1, 2), (4, 5)], [(0, 2), (2, 3), (5, 6)], [(0, 2), (3, 4), (5, 6)]),
            ([2, (4, 5)], [(1, 2), (2, 3)], [(1, 2), (3, 4)]),
            ([(1, 2), (2, 3)], [2, (0, 1)], [2, (4, 5), (5, 6)], [1, (0, 2), (3, 4), (4, 6)]),
            ([(0, 2), (2, 3)], [2, (4, 6)], [(0, 2), (3, 4), (5, 6)]),
            ([(1, 2), (2, 3), (4, 6)], [(1, 2), (3, 4), (4, 6)]),
            ([(1, 2), (4, 6)], [(0, 2), (2, 3)], [(0, 2), (3, 4), (5, 6)]),
            ([(0, 2), (2, 3), (4, 6)], [(0, 2), (3, 5), (5, 6)]),
            ([(0, 1), (1, 2)], [2, (4, 5), (5, 6)], [(1, 2), (2, 3)], [(0, 2), (3, 4), (5, 6)]),
            ([2, (1, 3), (4, 6)], [2, (0, 1)], [1, 4, (0, 2), (3, 6)]),
            ([(1, 3), (4, 6)], [(0, 1)], [1, (0, 2), (3, 6)]),
            ([(1, 3), (4, 5)], [(1, 2)], [2, (0, 1)], [1, (0, 2), (3, 4)]),
            ([1, (2, 3)], [(1, 2), (4, 5)]),
            ([2, (0, 1), (1, 3)], [4, (0, 1), (1, 2)]),
            ([(0, 1), (1, 3)], [(0, 1), (1, 2)]),
            ([(1, 3)], [(4, 5)], [(0, 1), (1, 2)], [1, (0, 2), (3, 4)]),
            ([(1, 2), (2, 3), (4, 5)], [(1, 2), (3, 4), (4, 5)]),
            ([2, (0, 1), (1, 3), (4, 6)], [1, (0, 2), (3, 4), (4, 6)]),
            ([(1, 3)], [(0, 2), (4, 6)], [(0, 2), (3, 4), (5, 6)]),
            ([(1, 3), (4, 5)], [(1, 2)], [(1, 2), (3, 4)]),
            ([2, (0, 1)], [2, (4, 5), (5, 6)], [(1, 2), (2, 3)], [(0, 2), (3, 4), (5, 6)]),
            ([(0, 1), (1, 2), (2, 3)], [2, (5, 6)], [(0, 2), (4, 5), (5, 6)]),
            ([(1, 2), (4, 5)], [(0, 1), (1, 2), (2, 3)], [2, (5, 6)], [(0, 2), (3, 4), (5, 6)]),
            ([(4, 5)], [(1, 2)], [(0, 2), (2, 3)], [(0, 2), (3, 4)]),
            ([(1, 2), (4, 5)], [(1, 2), (2, 3)], [2, (0, 1)], [1, (0, 2), (3, 4)]),
            ([(1, 3), (4, 5)], [(0, 1), (5, 6)], [1, (0, 2), (3, 6)]),
            ([2, (1, 3), (4, 5)], [2, (0, 1), (5, 6)], [1, 4, (0, 2), (3, 6)]),
            ([(1, 3), (4, 6)], [(0, 2)], [1, (0, 2), (3, 6)]),
            ([2, (4, 5), (5, 6)], [(1, 2), (2, 3)], [(1, 2), (3, 4), (4, 6)]),
            ([2, (4, 5)], [1, (0, 2), (2, 3)], [1, (0, 2), (3, 4)]),
            ([(0, 1), (1, 3)], [1, (4, 6)], [1, (0, 2), (3, 6)]),
            ([2, (0, 1), (1, 3)], [1, 2, (4, 6)], [1, 4, (0, 2), (3, 6)]),
            ([2, (4, 5), (5, 6)], [1, (2, 3)], [(0, 1), (1, 2)], [1, (0, 2), (3, 4), (4, 6)]),
            ([2, (4, 5), (5, 6)], [2, (0, 1), (1, 3)], [1, (0, 2), (3, 4), (4, 6)]),
            ([(1, 2)], [2, (4, 5)], [(0, 1), (1, 2), (2, 3)], [1, (0, 2), (3, 4)]),
            ([(0, 2), (2, 3), (4, 5)], [2, (5, 6)], [(0, 2), (3, 4), (5, 6)]),
            ([(2, 3), (4, 5), (5, 6)], [1, (0, 2)], [1, (0, 2), (3, 4), (4, 6)]),
            ([(1, 2)], [(0, 2), (2, 3), (4, 6)], [(0, 2), (3, 4), (5, 6)]),
            ([1, (0, 3)], [1, (0, 2)]),
            ([1, 2, (0, 3)], [1, 4, (0, 2)]),
            ([(1, 3), (4, 6)], [(0, 2)], [(0, 2), (3, 4), (5, 6)]),
            ([1, (0, 2), (2, 3), (4, 6)], [1, (0, 2), (3, 4), (4, 6)]),
            ([(4, 5), (5, 6)], [(1, 2)], [(0, 2), (2, 3)], [(0, 2), (3, 4), (5, 6)]),
            ([(1, 3), (4, 6)], [(1, 2), (3, 6)]),
            ([2, (1, 3), (4, 6)], [4, (1, 2), (3, 6)]),
            ([(4, 5)], [(0, 1), (1, 2)], [(1, 2), (2, 3)], [(0, 2), (3, 4)]),
            ([2, (1, 3)], [2, (0, 1), (4, 6)], [1, 4, (0, 2), (3, 6)]),
            ([(1, 3)], [(0, 1), (4, 6)], [1, (0, 2), (3, 6)]),
            ([1, 2, (4, 6)], [(0, 1), (1, 2), (2, 3)], [1, (0, 2), (3, 4), (4, 6)]),
            ([(0, 1), (1, 2)], [1, (2, 3)], [1, (0, 2)]),
            ([2, (0, 1)], [(2, 3), (5, 6)], [(0, 2), (5, 6)]),
            ([(2, 3), (5, 6)], [(0, 2), (4, 5)], [1, (0, 2), (3, 4), (4, 6)]),
            ([(0, 2), (2, 3), (4, 5), (5, 6)], [(0, 2), (3, 4), (4, 5), (5, 6)]),
            ([(1, 3), (4, 5)], [(0, 1), (1, 2)], [2, (5, 6)], [1, (0, 2), (3, 6)]),
            ([(1, 2)], [(0, 1), (2, 3)], [(4, 5), (5, 6)], [1, (0, 2), (3, 6)]),
            ([(4, 5)], [(1, 3)], [(0, 2)], [(0, 2), (3, 4)]),
            ([(0, 1), (1, 3)], [(1, 2), (4, 6)], [(0, 2), (3, 4), (5, 6)]),
            ([(1, 2), (2, 3)], [(1, 2), (4, 5)]),
            ([(2, 3), (4, 5), (5, 6)], [(0, 1), (1, 2)], [(0, 1), (1, 2), (3, 4), (4, 6)]),
            ([(1, 3), (4, 5)], [(0, 1), (1, 2)], [2, (5, 6)], [(0, 2), (3, 4), (5, 6)]),
            ([1, (2, 3)], [(0, 1), (5, 6)], [1, (0, 2), (4, 6)]),
            ([(2, 3), (4, 6)], [(1, 2), (3, 6)]),
            ([(2, 3), (4, 6)], [(0, 1)], [1, (0, 2), (3, 6)]),
            ([2, (4, 5), (5, 6)], [(0, 1), (1, 2), (2, 3)], [(0, 2), (3, 4), (4, 5), (5, 6)]),
            ([(0, 1), (1, 2)], [(2, 3), (4, 5), (5, 6)], [(1, 2)], [(0, 2), (3, 4), (5, 6)]),
            ([(0, 1), (2, 3), (4, 6)], [(0, 2), (3, 5), (5, 6)]),
            ([(0, 3), (4, 5)], [(1, 2)], [(0, 2), (3, 4)]),
            ([(1, 2), (4, 5)], [2, (0, 1), (1, 3)], [1, (0, 2), (3, 4)]),
            ([(1, 3)], [(0, 2), (4, 5)], [1, (0, 2), (3, 4)]),
            ([(0, 1), (2, 3)], [(5, 6)], [(0, 2), (5, 6)]),
            ([(1, 3)], [(1, 2)]),
            ([2, (1, 3)], [4, (1, 2)]),
            ([(2, 3)], [(0, 1)], [1, (0, 2)]),
            ([(2, 3), (4, 5)], [1, (0, 2)], [1, (0, 2), (3, 4)]),
            ([2, (0, 1)], [(2, 3), (4, 5), (5, 6)], [(1, 2)], [(0, 2), (3, 4), (5, 6)]),
            ([2, (4, 5)], [(0, 1), (1, 2), (2, 3)], [(0, 1), (1, 2), (3, 4)]),
            ([(0, 2), (2, 3)], [(0, 2), (4, 5)]),
            ([1, (2, 3)], [(0, 2), (4, 5), (5, 6)], [1, (0, 2), (3, 4), (4, 6)]),
            ([2, (0, 1), (1, 3)], [2, (4, 5), (5, 6)], [1, (0, 2), (3, 4), (4, 6)]),
            ([(1, 2)], [(0, 1), (2, 3)], [1, (0, 2)]),
            ([(0, 2), (4, 5), (5, 6)], [(1, 2), (2, 3)], [(0, 2), (3, 4), (5, 6)]),
            ([(0, 2), (2, 3), (4, 5)], [(0, 2), (3, 4), (4, 5)]),
            ([(4, 5)], [(0, 1), (1, 3)], [(1, 2)], [(0, 2), (3, 4)]),
            ([(0, 1), (2, 3)], [(4, 5), (5, 6)], [(0, 1), (1, 2), (3, 6)]),
            ([(0, 2), (2, 3), (5, 6)], [(0, 2), (4, 5), (5, 6)]),
            ([(0, 1), (1, 3)], [(1, 2), (4, 6)], [1, (0, 2), (3, 6)]),
            ([(2, 3), (5, 6)], [(0, 1)], [1, (0, 2), (4, 6)]),
            ([(0, 1), (2, 3)], [(4, 5), (5, 6)], [(1, 2)], [1, (0, 2), (3, 6)]),
            ([(2, 3), (4, 5), (5, 6)], [(1, 2)], [(1, 2), (3, 4), (4, 6)]),
            ([(2, 3), (5, 6)], [(1, 2), (4, 6)]),
            ([(0, 3), (4, 5), (5, 6)], [(1, 2)], [(0, 2), (3, 4), (5, 6)]),
            ([(2, 3)], [(0, 2), (4, 5)], [1, (0, 2), (3, 4)]),
            ([(1, 3), (4, 5)], [(0, 2), (5, 6)], [(0, 2), (3, 4), (5, 6)]),
            ([(2, 3), (5, 6)], [(0, 2)], [1, (0, 2), (4, 6)]),
            ([(1, 3)], [(0, 2), (4, 6)], [1, (0, 2), (3, 6)]),
            ([(2, 3), (4, 5)], [(0, 1), (5, 6)], [1, (0, 2), (3, 6)]),
            ([(2, 3), (4, 5)], [(0, 1), (1, 2)], [(0, 1), (1, 2), (3, 4)]),
            ([(1, 3), (4, 5)], [(0, 2), (5, 6)], [1, (0, 2), (3, 6)]),
            ([(0, 1), (2, 3)], [(0, 2), (4, 5)]),
            ([(1, 2), (4, 5)], [(1, 2), (2, 3)], [(1, 2), (3, 4)]),
            ([(0, 1), (1, 3)], [(1, 2)], [1, (0, 2)]),
            ([(2, 3)], [(0, 2), (4, 6)], [1, (0, 2), (3, 6)]),
            ([(2, 3)], [(0, 2)], [1, (0, 2)]),
            ([2, (4, 5)], [(0, 1), (1, 2), (2, 3)], [2, (5, 6)], [1, (0, 2), (3, 6)]),
            ([(0, 1), (2, 3)], [(5, 6)], [(4, 5)], [1, (0, 2), (3, 6)]),
            ([(1, 2), (2, 3)], [2, (0, 1)], [2, (5, 6)], [1, (0, 2), (4, 6)]),
            ([2, (4, 5), (5, 6)], [(1, 2), (2, 3)], [2, (0, 1)], [1, (0, 2), (3, 4), (4, 6)]),
            ([(1, 3)], [(0, 2)], [1, (0, 2)]),
            ([1, (2, 3)], [(0, 1), (1, 2)], [2, (4, 5), (5, 6)], [1, (0, 2), (3, 4), (4, 6)]),
            ([(4, 5), (5, 6)], [(0, 1), (1, 3)], [(1, 2)], [(0, 2), (3, 4), (5, 6)]),
            ([(0, 1), (2, 3)], [(1, 2), (4, 6)], [(0, 2), (3, 4), (5, 6)]),
            ([(0, 1), (1, 2), (2, 3)], [1, 2, (4, 6)], [1, (0, 2), (3, 4), (4, 6)]),
            ([(2, 3), (4, 5)], [(1, 2)], [(1, 2), (3, 4)]),
            ([(0, 1), (1, 2)], [(2, 3), (4, 5), (5, 6)], [(0, 2), (3, 5), (5, 6)]),
            ([(4, 5), (5, 6)], [(1, 3)], [(0, 2)], [(0, 2), (3, 4), (5, 6)]),
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

        # Combine some of the one-step arrows
        combined_one_step_arrows = []  # list of lists of arrows
        while len(one_step_arrows) != 0:
            arrow = one_step_arrows.pop()
            coeff_d, coeffs_a = arrow
            valid = True
            for idem_d, idem_a in single_idems:
                if idem_a in coeff_d.single_hor and \
                   not all(
                       idem_d in coeff_a.single_hor for coeff_a in coeffs_a):
                    valid = False
                    break
            if not valid:
                continue

            alt_idems = ArcslideDA.getAltIdempotents([arrow], single_idems)
            for alt_idem in alt_idems:
                one_step_arrows.remove(alt_idem)
            combined_one_step_arrows.append([arrow] + alt_idems)

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

        lin_sys = F2RowSystem(matrix_map)
        comb = lin_sys.getComb(target_vec)
        assert comb is not None

        result = []
        for i in range(len(combined_one_step_arrows)):
            if comb[i] != 0 and \
               combined_one_step_arrows[i][0] not in arrows_new:
                result.extend(combined_one_step_arrows[i])
        return result
