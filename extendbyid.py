"""Extension by identity of type DA structures."""

from algebra import TensorGenerator
from algebra import E0
from dastructure import DAStructure, SimpleDAGenerator
from localpmc import restrictIdempotent, restrictStrandDiagram
from utility import ACTION_LEFT, ACTION_RIGHT, F2

class ExtendedDAStructure(DAStructure):
    """Type DA structure obtained by extension by identity from a local type DA
    structure.

    """
    def __init__(self, local_da, outer_pmc, pmc1, pmc2,
                 mapping1, outer_mapping1, mapping2, outer_mapping2):
        """Specifies the local type DA structure (local_da, of type
        DAStructure), and the following data needed for an extension:
         - outer_pmc: of type LocalPMC, the local PMC on the outside.
         - pmc1, mapping1, outer_mapping1: full PMC on the D-side, and mappings
           specifying the decomposition into local PMCs.
         - pmc2, mapping2, outer_mapping2: full PMC on the A-side, and mappings
           specifying the decomposition into local PMCs.
        
        """
        self.local_da = local_da
        self.outer_pmc = outer_pmc
        self.pmc1, self.pmc2 = pmc1, pmc2
        self.mapping1, self.outer_mapping1 = mapping1, outer_mapping1
        self.mapping2, self.outer_mapping2 = mapping2, outer_mapping2

        self.local_pmc1 = local_da.algebra1.local_pmc
        self.local_pmc2 = local_da.algebra2.local_pmc

        # Initiate the DA structure
        DAStructure.__init__(self, F2, algebra1 = self.pmc1.getAlgebra(),
                             algebra2 = self.pmc2.getAlgebra(),
                             side1 = ACTION_LEFT, side2 = ACTION_RIGHT)

        local_gens = self.local_da.getGenerators()
        # Verify that each generator in local_da has a unique idempotent.
        # Build a dictionary mapping local idempotents to local generators.
        self.local_idem_to_gen = dict()
        for local_gen in local_gens:
            local_idem = (local_gen.idem1, local_gen.idem2)
            assert local_idem not in self.local_idem_to_gen
            self.local_idem_to_gen[local_idem] = local_gen

        # Obtain the set of extended generators
        self.generators = []
        idem1_list = self.algebra1.getIdempotents()
        idem2_list = self.algebra2.getIdempotents()
        for idem1 in idem1_list:
            for idem2 in idem2_list:
                # Restrict idem1 and idem2 to local idempotents
                local_idem = (
                    restrictIdempotent(pmc1, idem1, self.local_pmc1, mapping1),
                    restrictIdempotent(pmc2, idem2, self.local_pmc2, mapping2))
                outside_idem = (
                    restrictIdempotent(
                        pmc1, idem1, self.outer_pmc, outer_mapping1),
                    restrictIdempotent(
                        pmc2, idem2, self.outer_pmc, outer_mapping2))
                outside_idem = [idem.removeSingleHor() for idem in outside_idem]
                if outside_idem[0] == outside_idem[1] and \
                   local_idem in self.local_idem_to_gen:
                    cur_gen = SimpleDAGenerator(
                        self, idem1, idem2, "%d" % len(self.generators))
                    self.generators.append(cur_gen)

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
    def adjustSingleHors(init_idems, coeffs_a):
        """Given the idempotent (both sides) of the initial (local) generator
        and a tuple of (local) A-side inputs, delete single idempotents from
        inputs in a best effort to make the idempotents match.

        If failed, return (None, None)

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
            if not changed:
                return (None, None)
        init_l_idem, init_r_idem = init_idems
        if len(coeffs_a) > 0 and init_r_idem != coeffs_a[0].left_idem:
            for idem in init_r_idem:
                if idem not in coeffs_a[0].left_idem:
                    init_r_idem = init_r_idem.removeSingleHor([idem])
                    # Warning: very bad when there are multiple single
                    # horizontals. Need to fix.
                    init_l_idem = init_l_idem.removeSingleHor()
        return ((init_l_idem, init_r_idem), tuple(coeffs_a))

    def getGenerators(self):
        return self.generators

    def delta(self, MGen, algGens):
        # Idempotent must match
        if any([algGens[i].right_idem != algGens[i+1].left_idem
                for i in range(len(algGens)-1)]):
            return E0
        if any([alg.isIdempotent() for alg in algGens]):
            return E0

        result = E0
        mod_gens = self.getGenerators()

        # Multiply the outside
        has_product = True
        prod_d = restrictIdempotent(self.pmc1, MGen.idem1, self.outer_pmc,
                                    self.outer_mapping1).toAlgElt()
        for alg in algGens:
            prod_d = prod_d.parent.multiplyGeneral(
                prod_d, restrictStrandDiagram(
                    self.pmc2, alg, self.outer_pmc, self.outer_mapping2),
                False)  # strict_idems = False
            if prod_d == 0:
                has_product = False
                break
            else:
                prod_d = prod_d.getElt()
        if not has_product:
            return E0

        # Test the inside
        idem_local = (
            restrictIdempotent(
                self.pmc1, MGen.idem1, self.local_pmc1, self.mapping1),
            restrictIdempotent(
                self.pmc2, MGen.idem2, self.local_pmc2, self.mapping2))
        alg_local = tuple([restrictStrandDiagram(
            self.pmc2, alg, self.local_pmc2, self.mapping2) for alg in algGens])
        (idem_local, alg_local) = ExtendedDAStructure.adjustSingleHors(
            idem_local, alg_local)
        local_MGen = self.local_idem_to_gen[idem_local]
        local_delta = self.local_da.delta(local_MGen, alg_local)
        for (local_d, local_y), ring_coeff in local_delta.items():
            alg_d = local_d.join(prod_d.removeSingleHor(), self.pmc1,
                                 self.mapping1, self.outer_mapping1)
            if alg_d is None:
                continue
            for y in mod_gens:
                if ExtendedDAStructure.idemMatchDA(MGen, y, alg_d, algGens):
                    result += 1 * TensorGenerator((alg_d, y), self.AtensorM)
        return result

    def deltaPrefix(self, MGen, algGens):
        if len(algGens) == 0:
            return True
        if MGen.idem2 != algGens[0].left_idem:
            return False

        # Note: Testing shows probably not worth it to multipliy outside.
        # Just test inside.
        idem_local = (
            restrictIdempotent(
                self.pmc1, MGen.idem1, self.local_pmc1, self.mapping1),
            restrictIdempotent(
                self.pmc2, MGen.idem2, self.local_pmc2, self.mapping2))
        alg_local = tuple([restrictStrandDiagram(
            self.pmc2, alg, self.local_pmc2, self.mapping2) for alg in algGens])
        (idem_local, alg_local) = ExtendedDAStructure.adjustSingleHors(
            idem_local, alg_local)
        if idem_local not in self.local_idem_to_gen:
            return False
        local_MGen = self.local_idem_to_gen[idem_local]
        if self.local_da.deltaPrefix(local_MGen, alg_local):
            return True

        # One more situation: if idem_local and alg_local has single idempotent
        # throughout, try the case without the single idempotent.
        # Note: fix when there are more than one single idempotents.
        if alg_local[0].removeSingleHor() != alg_local[0]:
            alg_local = (alg_local[0].removeSingleHor(),) + alg_local[1:]
            (idem_local, alg_local) = ExtendedDAStructure.adjustSingleHors(
                idem_local, alg_local)
            if idem_local not in self.local_idem_to_gen:
                return False
            local_MGen = self.local_idem_to_gen[idem_local]
            if self.local_da.deltaPrefix(local_MGen, alg_local):
                return True
        return False
