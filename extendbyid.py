"""Extension by identity of type DA structures."""

from algebra import TensorGenerator
from algebra import E0
from dastructure import DAGenerator, DAStructure, SimpleDAGenerator
from localpmc import PMCSplitting
from utility import ACTION_LEFT, ACTION_RIGHT, F2

class ExtendedDAGenerator(SimpleDAGenerator):
    """Represents a generator of the extended DA structure. Stores the generator
    of the local DA structure this comes from, and the outer idempotent.

    """
    def __init__(self, parent, local_gen, outer_idem, name):
        assert local_gen.parent == parent.local_da
        assert outer_idem.local_pmc == parent.outer_pmc
        self.local_gen = local_gen
        self.outer_idem = outer_idem
        idem1 = parent.splitting1.joinIdempotent(local_gen.idem1, outer_idem)
        idem2 = parent.splitting2.joinIdempotent(local_gen.idem2, outer_idem)
        SimpleDAGenerator.__init__(self, parent, idem1, idem2, name)

class ExtendedDAStructure(DAStructure):
    """Type DA structure obtained by extension by identity from a local type DA
    structure.

    """
    def __init__(self, local_da, splitting1, splitting2):
        """Specifies the local type DA structure (local_da, of type
        DAStructure), and the following data needed for an extension:
         - outer_pmc: of type LocalPMC, the local PMC on the outside.
         - pmc1, mapping1, outer_mapping1: full PMC on the D-side, and mappings
           specifying the decomposition into local PMCs.
         - pmc2, mapping2, outer_mapping2: full PMC on the A-side, and mappings
           specifying the decomposition into local PMCs.
        
        """
        self.local_da = local_da
        self.splitting1 = splitting1
        self.splitting2 = splitting2

        self.pmc1, self.pmc2 = splitting1.pmc, splitting2.pmc

        self.outer_pmc = splitting1.outer_pmc
        assert self.outer_pmc == splitting2.outer_pmc

        self.local_pmc1 = local_da.algebra1.local_pmc
        self.local_pmc2 = local_da.algebra2.local_pmc
        assert self.local_pmc1 == splitting1.local_pmc
        assert self.local_pmc2 == splitting2.local_pmc

        mapping1 = splitting1.local_mapping
        mapping2 = splitting2.local_mapping

        self.idem_size1 = self.pmc1.genus
        self.idem_size2 = self.pmc2.genus

        # Possible values of single assignments
        self.NONE, self.LOCAL, self.OUTER = 0, 1, 2

        def checkSingle(local_pmc):
            # Checks the local PMC has a single unpaired point. Return the pair
            # index of that point.
            pairs = local_pmc.getSingleIdems()
            assert len(pairs) == 1
            return pairs[0]

        # Index of the pair / single point in local_pmcs and outer_pmc that has
        # a single point.
        self.single1 = checkSingle(self.local_pmc1)
        self.single2 = checkSingle(self.local_pmc2)
        self.single_outer = checkSingle(self.outer_pmc)
        # Index of pair in each of the two full PMCs that is smeared.
        for p in range(self.pmc1.n):
            if p in mapping1 and \
                    self.local_pmc1.pairid[mapping1[p]] == self.single1:
                self.smeared1 = self.pmc1.pairid[p]
        for p in range(self.pmc2.n):
            if p in mapping2 and \
                    self.local_pmc2.pairid[mapping2[p]] == self.single2:
                self.smeared2 = self.pmc2.pairid[p]

        # Initiate the DA structure
        DAStructure.__init__(self, F2, algebra1 = self.pmc1.getAlgebra(),
                             algebra2 = self.pmc2.getAlgebra(),
                             side1 = ACTION_LEFT, side2 = ACTION_RIGHT)

        # Obtain the set of extended generators
        self.generators = []
        local_gens = self.local_da.getGenerators()
        outer_idems = [idem for idem in self.outer_pmc.getIdempotents()
                       if self.single_outer not in idem]
        # Maps (local_gen, outer_idem) to generator
        self.gen_index = dict()
        for local_gen in local_gens:
            for outer_idem in outer_idems:
                if len(local_gen.idem1) + len(outer_idem) != self.idem_size1:
                    continue
                assert len(local_gen.idem2) + len(outer_idem) == self.idem_size2
                cur_gen = ExtendedDAGenerator(self, local_gen, outer_idem,
                                              "%d" % len(self.generators))
                self.generators.append(cur_gen)
                self.gen_index[(local_gen, outer_idem)] = cur_gen

        # Obtain the u map on local generators. Currently assuming that
        # generators have unique idempotents.
        self.u_map = dict()
        self.uinv_map = dict()
        for local_gen in local_gens:
            idem1, idem2 = local_gen.idem1, local_gen.idem2
            if self.single1 in idem1 and self.single2 in idem2:
                target_idem1 = idem1.removeSingleHor()
                target_idem2 = idem2.removeSingleHor()
                for target_gen in local_gens:
                    if target_gen.idem1 == target_idem1 and \
                            target_gen.idem2 == target_idem2:
                        self.u_map[local_gen] = target_gen
                        self.uinv_map[target_gen] = local_gen

    def getGenerators(self):
        return self.generators

    def getSingleAssignments(self, algGens):
        """Compute how to assign the smeared pair to local or outside. Return
        None if there is no valid assignment.

        """
        assignments = [self.NONE] * len(algGens)
        # Enforce as much as possible constraints from the previous algebra
        # input.
        for i in range(1, len(algGens)):
            alg_prev, alg = algGens[i-1], algGens[i]
            if self.smeared2 in alg.double_hor:
                # Note assignments[i] will stay NONE if there is not enough
                # information.
                if assignments[i-1] != self.NONE:
                    assignments[i] = assignments[i-1]
                elif self.smeared2 not in alg_prev.double_hor:
                    for q in [q for p, q in alg_prev.strands
                              if self.pmc2.pairid[q] == self.smeared2]:
                        if q in self.mapping2:
                            assignments[i] = self.LOCAL
                        else:
                            assignments[i] = self.OUTER
        # Now enforce constraints from next algebra input. Return None if there
        # is a contradiction.
        for i in reversed(range(len(algGens)-1)):  # need test case
            alg_next, alg = algGens[i+1], algGens[i]
            if self.smeared2 in alg.double_hor:
                to_assign = self.NONE
                if assignments[i+1] != self.NONE:
                    to_assign = assignments[i+1]
                elif self.smeared2 not in alg_next.double_hor:
                    for p in [p for p, q in alg_next.strands
                              if self.pmc2.pairid[p] == self.smeared2]:
                        if p in self.mapping2:
                            to_assign = self.LOCAL
                        else:
                            to_assign = self.OUTER
                if assignments[i] != self.NONE and to_assign != self.NONE and \
                        assignments[i] != to_assign:
                    return None  # need test case
                if to_assign != self.NONE:
                    assignments[i] = to_assign
        # There is one final case where everything is double horizontal. Assign
        # all to local (this agrees with how generators are created for the
        # extended bimodule, where the smeared pair always goes to local).
        if all([self.smeared2 in alg.double_hor for alg in algGens]):
            assignments = [self.LOCAL] * len(algGens)
        return assignments

    def delta(self, MGen, algGens):
        # Idempotent must match
        if len(algGens) > 0 and algGens[0].left_idem != MGen.idem2:
            return E0
        if any([algGens[i].right_idem != algGens[i+1].left_idem
                for i in range(len(algGens)-1)]):
            return E0
        if any([alg.isIdempotent() for alg in algGens]):
            return E0

        result = E0
        mod_gens = self.getGenerators()

        assignments = self.getSingleAssignments(algGens)
        if assignments is None:
            return E0

        # Restrict and multiply on the outside
        has_product = True
        prod_d = None
        for i in range(len(algGens)):
            alg = algGens[i]
            outer_sd = self.splitting2.restrictStrandDiagramOuter(alg)
            if assignments[i] == self.LOCAL:
                outer_sd = outer_sd.removeSingleHor()
            if prod_d is None:
                prod_d = 1*outer_sd
            else:
                prod_d = prod_d * outer_sd
            if prod_d == 0:
                has_product = False
                break
            else:
                prod_d = prod_d.getElt()
        if not has_product:
            return E0
        if prod_d is None:  # No algebra input
            prod_d = \
                self.splitting2.restrictIdempotentOuter(MGen.idem2).toAlgElt()
            prod_d = prod_d.removeSingleHor()

        # Test the inside
        local_MGen = MGen.local_gen
        alg_local = []
        for i in range(len(algGens)):
            cur_alg_local = \
                self.splitting2.restrictStrandDiagramLocal(algGens[i])
            if assignments[i] == self.OUTER:
                cur_alg_local = cur_alg_local.removeSingleHor()
            alg_local.append(cur_alg_local)
        alg_local = tuple(alg_local)

        if self.single2 in local_MGen.idem2 and len(alg_local) > 0 and \
                self.single2 not in alg_local[0].left_idem:
            if local_MGen not in self.u_map:
                return E0  # need test case
            local_delta = self.local_da.delta(self.u_map[local_MGen], alg_local)
        else:
            local_delta = self.local_da.delta(local_MGen, alg_local)

        for (local_d, local_y), ring_coeff in local_delta.items():
            alg_d = self.splitting1.joinStrandDiagram(local_d, prod_d)
            if alg_d is None:
                continue
            if self.single_outer in prod_d.right_idem:
                if self.single1 not in local_y.idem1:
                    local_y = self.uinv_map[local_y]
                y = self.gen_index[(local_y,
                                    prod_d.right_idem.removeSingleHor())]
            else:
                y = self.gen_index[(local_y, prod_d.right_idem)]
            result += 1 * TensorGenerator((alg_d, y), self.AtensorM)
        return result

    def deltaPrefix(self, MGen, algGens):
        if len(algGens) == 0:
            return True
        if algGens[0].left_idem != MGen.idem2:
            return False
        if any([algGens[i].right_idem != algGens[i+1].left_idem
                for i in range(len(algGens)-1)]):
            return False

        # Note: Testing shows probably not worth it to multiply outside.
        # Just test inside.
        def testWithAssignment(assignments):
            # Test with a specific single assignments
            local_MGen = MGen.local_gen
            alg_local = []
            for i in range(len(algGens)):
                cur_alg_local = \
                    self.splitting2.restrictStrandDiagramLocal(algGens[i])
                if assignments[i] == self.OUTER:
                    cur_alg_local = cur_alg_local.removeSingleHor()
                alg_local.append(cur_alg_local)
            alg_local = tuple(alg_local)
            if self.single2 in local_MGen.idem2 and len(alg_local) > 0 and \
                    self.single2 not in alg_local[0].left_idem:
                if local_MGen not in self.u_map:
                    return False  # need test case
                return self.local_da.deltaPrefix(self.u_map[local_MGen],
                                                 alg_local)
            else:
                return self.local_da.deltaPrefix(local_MGen, alg_local)

        assignments = self.getSingleAssignments(algGens)
        if assignments is None:
            return False

        if all([a == self.LOCAL for a in assignments]):
            return testWithAssignment(assignments) or \
                testWithAssignment([self.OUTER] * len(algGens)) # need test case
        else:
            return testWithAssignment(assignments)
