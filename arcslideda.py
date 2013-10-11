"""Producing type DA structures for arcslides, using local actions."""

from arcslide import Arcslide
from dastructure import SimpleDAGenerator, SimpleDAStructure
from dastructure import AddChordToDA
from localpmc import LocalStrandDiagram
from localpmc import restrictPMC, restrictStrandDiagram
from pmc import Strands
from utility import F2

# debug only
from pmc import splitPMC

class ArcslideDA:
    """Responsible for producing a type DA structure for an arcslide, using
    local actions.

    """
    @staticmethod
    def getDAStructure(slide):
        """Returns the type DA structure corresponding to slide (type Arcslide).

        """
        genus = slide.start_pmc.genus
        assert slide == Arcslide(splitPMC(genus), 1, 2)
        dd_idems = slide.getIdems()
        da_idems = [(l_idem, r_idem.opp().comp())
                    for l_idem, r_idem in dd_idems]
        pmc = splitPMC(genus)
        alg = pmc.getAlgebra()
        dastr = SimpleDAStructure(F2, alg, alg)
        for i in range(len(da_idems)):
            l_idem, r_idem = da_idems[i]
            dastr.addGenerator(
                SimpleDAGenerator(dastr, l_idem, r_idem, "%d" % i))

        alg_gens = alg.getGenerators()
        mod_gens = dastr.getGenerators()
        local_pmc, mapping = restrictPMC(pmc, [(0, 2)])
        outer_pmc, outer_mapping = restrictPMC(pmc, [(3, 4*genus-1)])

        single_patterns = {}
        # Format: A-side idem, A-side strands, D-side idem, D-side strands
        single_patterns_raw = [
            ([], [], [], []),
            ([0], [], [0], []),
            ([0], [], [1], []),
            ([1], [], [1], []),
            ([1], [], [0], []),
            ([1], [(1, 2)], [0], []),
            ([0], [(0, 1)], [0], [(0, 2)]),
            ([0], [(2, 3)], [0], [(2, 3)]),
            ([0, 1], [(2, 3)], [0, 1], [(2, 3)]),
            ([0], [(0, 2)], [0], [(0, 2)]),
            ([0, 1], [(0, 2)], [0, 1], [(0, 2)]),
            ([0, 1], [(0, 1),(1, 2)], [0, 1], [(0, 1),(1, 2)]),
            ([1], [(1, 3)], [0], [(2, 3)]),
            ([0], [(0, 3)], [0], [(0, 3)]),
            ([0, 1], [(0, 3)], [0, 1], [(0, 3)]),
            ([0, 1], [(0, 1),(1, 3)], [0, 1], [(0, 1),(1, 3)]),
        ]
        for a_idem, a_strands, d_idem, d_strands in single_patterns_raw:
            key = LocalStrandDiagram(local_pmc, a_idem, a_strands)
            if key not in single_patterns:
                single_patterns[key] = []
            single_patterns[key].append(
                LocalStrandDiagram(local_pmc, d_idem, d_strands))

        # Format: A-side1 idem, A-side1 strands, A-side2 idem, A-side2 strands,
        # D-side idem, D-side strands
        double_patterns_raw = [
            ([1], [(1, 2)], [0], [(0, 1)], [1], [(1, 2)]),
            ([1], [(1, 2)], [0], [(0, 2)], [1], [(1, 2)]),
            ([1], [(1, 2)], [0], [(0, 3)], [1], [(1, 3)]),
            ([0, 1], [(1, 2),(2, 3)], [0], [(0, 1)], [0, 1], [(1, 2),(2, 3)]),
            ([0, 1], [(1, 2),(2, 3)], [0], [(0, 2)], [0, 1], [(1, 2),(2, 3)]),
            ([0, 1], [(1, 2),(2, 3)], [0, 1], [(0, 2)],
             [0, 1], [(1, 2),(2, 3)]),
            ([1], [(1, 3)], [0], [(0, 1)], [1], [(1, 3)]),
            ([0, 1], [(1, 3)], [0], [(0, 1)], [0, 1], [(1, 3)]),
            ([1], [(1, 3)], [0], [(0, 2)], [1], [(1, 3)]),
            ([0, 1], [(1, 3)], [0, 1], [(0, 2)], [0, 1], [(1, 3)]),
            ([0, 1], [(1, 3)], [0], [(0, 2)], [0, 1], [(1, 3)]),
        ]
        double_patterns = {}
        for a1_idem, a1_strands, a2_idem, a2_strands, d_idem, d_strands in \
            double_patterns_raw:
            key = (LocalStrandDiagram(local_pmc, a1_idem, a1_strands),
                   LocalStrandDiagram(local_pmc, a2_idem, a2_strands))
            double_patterns[key] = LocalStrandDiagram(
                local_pmc, d_idem, d_strands)

        # Add action with zero algebra input
        AddChordToDA(dastr, Strands(pmc, [(0, 1)]), [])

        for a1 in alg_gens:
            # process single patterns
            if a1.isIdempotent():
                continue
            local_a1 = restrictStrandDiagram(pmc, a1, local_pmc, mapping)
            outer_a1 = restrictStrandDiagram(pmc, a1, outer_pmc, outer_mapping)
            outer_a1 = outer_a1.removeSingleHor()
            if local_a1 in single_patterns:
                for local_d in single_patterns[local_a1]:
                    alg_d = local_d.join(outer_a1, pmc, mapping, outer_mapping)
                    if alg_d is not None:
                        for x in mod_gens:
                            for y in mod_gens:
                                if x.idem1 == alg_d.left_idem and \
                                   x.idem2 == a1.left_idem and \
                                   y.idem1 == alg_d.right_idem and \
                                   y.idem2 == a1.right_idem:
                                    dastr.addDelta(x, y, alg_d, [a1], 1)
            for a2 in alg_gens:
                # process double patterns
                if a2.isIdempotent():
                    continue
                if a1.right_idem != a2.left_idem:
                    continue
                local_a2 = restrictStrandDiagram(pmc, a2, local_pmc, mapping)
                outer_a2 = restrictStrandDiagram(
                    pmc, a2, outer_pmc, outer_mapping)
                # print a1, a2
                # print local_a1, local_a2
                if (local_a1, local_a2) in double_patterns:
                    local_d = double_patterns[(local_a1, local_a2)]
                    outer_prod = outer_a1.multiply(outer_a2)
                    if outer_prod is None:
                        continue
                    alg_d = local_d.join(
                        outer_prod, pmc, mapping, outer_mapping)
                    if alg_d is not None:
                        for x in mod_gens:
                            for y in mod_gens:
                                if x.idem1 == alg_d.left_idem and \
                                   x.idem2 == a1.left_idem and \
                                   y.idem1 == alg_d.right_idem and \
                                   y.idem2 == a2.right_idem:
                                    dastr.addDelta(x, y, alg_d, [a1, a2], 1)
        return dastr
