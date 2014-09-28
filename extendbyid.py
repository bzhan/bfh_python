"""Extension by identity of type DA structures."""

from algebra import TensorGenerator
from algebra import E0
from dastructure import DAGenerator, DAStructure, MorDAtoDAComplex, \
    SimpleDAGenerator, SimpleDAStructure
from localpmc import LocalStrandAlgebra, PMCSplitting
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

class LocalDAStructure(SimpleDAStructure):
    """Represents a local type DA structure. So far we always assume that a
    local type DA structure is simple (delta map is explicitly given). The extra
    data is the map between single idempotents on the two sides, and u_maps
    between generators.

    """
    def __init__(self, ring, algebra1, algebra2,
                 side1 = ACTION_LEFT, side2 = ACTION_RIGHT,
                 single_idems1 = None, single_idems2 = None):
        """single_idems1 and single_idems2 are two lists that order the unpaired
        idempotents on the two sides. Idempotents on the two sides that appear
        in the same position correspond to each other. If there are 0 or 1
        unpaired idempotents, they can be omitted by specifying None. Otherwise
        they must be provided.

        """
        assert isinstance(algebra1, LocalStrandAlgebra)
        assert isinstance(algebra2, LocalStrandAlgebra)
        if single_idems1 is None or single_idems2 is None:
            self.single_idems1 = algebra1.local_pmc.getSingleIdems()
            self.single_idems2 = algebra2.local_pmc.getSingleIdems()
            assert len(self.single_idems1) < 2, \
                "There are more than two unpaired idempotents."
        else:
            assert tuple(sorted(algebra1.local_pmc.getSingleIdems())) == \
                tuple(sorted(single_idems1))
            assert tuple(sorted(algebra2.local_pmc.getSingleIdems())) == \
                tuple(sorted(single_idems2))
            self.single_idems1 = single_idems1
            self.single_idems2 = single_idems2

        self.num_single_idems = len(self.single_idems1)
        assert self.num_single_idems == len(self.single_idems2)

        SimpleDAStructure.__init__(self, ring, algebra1, algebra2, side1, side2)
        self.u_maps = [dict() for i in range(self.num_single_idems)]
        self.uinv_maps = [dict() for i in range(self.num_single_idems)]

    def add_u_map(self, idem_id, source, target):
        """Add to the u-map a mapping from source to target. idem_id refers to
        the position in single_idems1 and single_idems2 given in the
        constructor.

        All entries in u-map remove the corresponding idempotents from idem1 and
        idem2.

        """
        self.u_maps[idem_id][source] = target
        self.uinv_maps[idem_id][target] = source

    def auto_u_map(self):
        """Autocompletes the u-maps. To call this function, one of the following
        must hold for each generator and each u-map for which it is eligible:
        1. The generator already appears as a key in the u_map.
        2. There is unique way to choose the target for this generator.

        """
        for i in range(self.num_single_idems):
            single1, single2 = self.single_idems1[i], self.single_idems2[i]
            for local_gen in self.generators:
                idem1, idem2 = local_gen.idem1, local_gen.idem2
                if single1 in idem1 and single2 in idem2:
                    # local_gen is eligible for u_maps[i]
                    if local_gen in self.u_maps[i]:
                        continue
                    # Otherwise, check there is a unique target and map there
                    target_idem1 = idem1.removeSingleHor([single1])
                    target_idem2 = idem2.removeSingleHor([single2])
                    target_gen = [gen for gen in self.generators
                                  if gen.idem1 == target_idem1
                                  and gen.idem2 == target_idem2]
                    assert len(target_gen) == 1, "Cannot autocomplete u-map"
                    self.add_u_map(i, local_gen, target_gen[0])

class LocalMorDAtoDAComplex(MorDAtoDAComplex):
    """Represents the complex of type DA morphisms between two local type DA
    structures.

    """
    def __init__(self, ring, source, target):
        assert isinstance(source, LocalDAStructure)
        assert isinstance(target, LocalDAStructure)
        MorDAtoDAComplex.__init__(self, ring, source, target)

    def getMappingCone(self, morphism):
        """In addition to what is done in the parent class, need to set up the
        u_map.

        """
        result = LocalDAStructure(
            F2, self.source.algebra1, self.source.algebra2,
            self.source.side1, self.source.side2,
            self.source.single_idems1, self.source.single_idems2)
        gen_map = dict()
        for gen in self.source.getGenerators():
            gen_map[gen] = SimpleDAGenerator(
                result, gen.idem1, gen.idem2, "S_%s" % gen.name)
            gen_map[gen].filtration = [0]
            if hasattr(gen, "filtration"):
                gen_map[gen] += gen.filtration
            result.addGenerator(gen_map[gen])
        for gen in self.target.getGenerators():
            gen_map[gen] = SimpleDAGenerator(
                result, gen.idem1, gen.idem2, "T_%s" % gen.name)
            gen_map[gen].filtration = [1]
            if hasattr(gen, "filtration"):
                gen_map[gen] += gen.filtration
            result.addGenerator(gen_map[gen])

        for (x1, coeffs_a), target in self.source.da_action.items():
            for (coeff_d, x2), ring_coeff in target.items():
                result.addDelta(
                    gen_map[x1], gen_map[x2], coeff_d, coeffs_a, ring_coeff)
        for (y1, coeffs_a), target in self.target.da_action.items():
            for (coeff_d, y2), ring_coeff in target.items():
                result.addDelta(
                    gen_map[y1], gen_map[y2], coeff_d, coeffs_a, ring_coeff)
        for gen, ring_coeff in morphism.items():
            # coeffs_a is a tuple of A-side inputs
            coeff_d, coeffs_a = gen.coeff
            result.addDelta(gen_map[gen.source], gen_map[gen.target],
                            coeff_d, tuple(coeffs_a), ring_coeff)

        # Set up u_map
        num_single_idems = len(self.source.single_idems1)
        for idem_id in range(num_single_idems):
            for x, u_x in self.source.u_maps[idem_id].items():
                result.add_u_map(idem_id, gen_map[x], gen_map[u_x])
            for y, u_y in self.target.u_maps[idem_id].items():
                result.add_u_map(idem_id, gen_map[y], gen_map[u_y])

        return result

class ExtendedDAStructure(DAStructure):
    """Type DA structure obtained by extension by identity from a local type DA
    structure.

    """
    def __init__(self, local_da, splitting1, splitting2):
        """Specifies the local type DA structure (local_da, of type
        DAStructure), and splittings of the two full PMCs on the two sides (of
        type PMCSplitting). The parameters should be consistent in the following
        way:

        self.local_pmc1 = splitting1.local_pmc = local_da.algebra1.local_pmc
        self.local_pmc2 = splitting2.local_pmc = local_da.algebra2.local_pmc
        self.outer_pmc = splitting1.outer_pmc = splitting2.outer_pmc
        
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

        self.mapping1 = splitting1.local_mapping
        self.mapping2 = splitting2.local_mapping
        self.outer_mapping1 = splitting1.outer_mapping
        self.outer_mapping2 = splitting2.outer_mapping

        self.idem_size1 = self.pmc1.genus
        self.idem_size2 = self.pmc2.genus

        # Possible values of single assignments, for use in delta and
        # deltaPrefix (through the function getSingleAssignments).
        self.NONE, self.LOCAL, self.OUTER = 0, 1, 2

        # Record the local and outer single idempotents.
        # Everything is indexed by 0 ... self.num_single-1
        self.single_idems1 = self.local_da.single_idems1  # idems in local_pmc1
        self.single_idems2 = self.local_da.single_idems2  # idems in local_pmc2
        self.num_singles = len(self.single_idems1)
        assert self.num_singles == len(self.single_idems2)
        self.smeared_idems1 = []  # idems in pmc1
        self.smeared_idems2 = []  # idems in pmc2
        self.single_idems_outer = []  # idems in outer_pmc
        self.single_pts_outer = []  # pts in outer_pmc

        for i in range(self.num_singles):
            # Fill in data, and verify that the correspondence of idempotents on
            # the two sides is consistent on the outer PMC.
            single_idems1 = self.single_idems1[i]
            single_idems2 = self.single_idems2[i]
            single_pt1 = self.local_pmc1.pairs[single_idems1][0]
            single_pt2 = self.local_pmc2.pairs[single_idems2][0]
            for p in range(self.pmc1.n):
                if p in self.mapping1 and self.mapping1[p] == single_pt1:
                    self.smeared_idems1.append(self.pmc1.pairid[p])
                    q = self.pmc1.otherp[p]
                    assert q in self.outer_mapping1
                    q_outer = self.outer_mapping1[q]
                    self.single_pts_outer.append(q_outer)
                    self.single_idems_outer.append(
                        self.outer_pmc.pairid[q_outer])
            for p in range(self.pmc2.n):
                if p in self.mapping2 and self.mapping2[p] == single_pt2:
                    self.smeared_idems2.append(self.pmc2.pairid[p])
                    q = self.pmc2.otherp[p]
                    assert q in self.outer_mapping2
                    assert self.single_pts_outer[-1] == self.outer_mapping2[q]
                    
        # Initiate the DA structure
        DAStructure.__init__(self, F2, algebra1 = self.pmc1.getAlgebra(),
                             algebra2 = self.pmc2.getAlgebra(),
                             side1 = ACTION_LEFT, side2 = ACTION_RIGHT)

        # Obtain the set of extended generators, and create a map self.gen_index
        # from (local_gen, outer_idem) to the extended generators.
        self.generators = []
        local_gens = self.local_da.getGenerators()
        outer_idems = [idem for idem in self.outer_pmc.getIdempotents()
                       if all(single_idem_outer not in idem for
                              single_idem_outer in self.single_idems_outer)]
        self.gen_index = dict()
        for local_gen in local_gens:
            cur_count = 0  # number of generators so far with local_gen
            for outer_idem in outer_idems:
                if len(local_gen.idem1) + len(outer_idem) != self.idem_size1:
                    continue
                assert len(local_gen.idem2) + len(outer_idem) == self.idem_size2
                cur_gen = ExtendedDAGenerator(
                    self, local_gen, outer_idem,
                    "%s%%%d" % (local_gen.name, cur_count))
                cur_count += 1
                if hasattr(local_gen, "filtration"):
                    cur_gen.filtration = local_gen.filtration
                self.generators.append(cur_gen)
                self.gen_index[(local_gen, outer_idem)] = cur_gen

    def __len__(self):
        return len(self.generators)

    def getGenerators(self):
        return self.generators

    def getSingleAssignments(self, algGens):
        """Compute how to assign the smeared pair to local or outside. Return
        None if there is no valid assignment.

        """
        # The final result. This will be a list with self.num_singles elements,
        # each element will be a list of length len(algGens), with each entry
        # one of self.LOCAL and self.OUTER (referring to the assignment of the
        # single horizontal line to either local or outer PMC.
        all_assignments = []

        for idem_id in self.smeared_idems2:
            # Current list of assignments.
            assignments = [self.NONE] * len(algGens)

            # Enforce as much as possible constraints from the previous algebra
            # input.
            for i in range(1, len(algGens)):
                alg_prev, alg = algGens[i-1], algGens[i]
                if idem_id in alg.double_hor:
                    # Note cur_assign[i] will stay NONE if there is not enough
                    # information.
                    if assignments[i-1] != self.NONE:
                        assignments[i] = assignments[i-1]
                    elif idem_id not in alg_prev.double_hor:
                        # There should be a moving strand in alg_prev ending at
                        # one of the points in the idempotent. Find whether the
                        # point is local or on the outside.
                        for t in [t for s, t in alg_prev.strands
                                  if self.pmc2.pairid[t] == idem_id]:
                            if t in self.mapping2:
                                assignments[i] = self.LOCAL
                            else:
                                assignments[i] = self.OUTER

            # Now enforce constraints from next algebra input. Return None if
            # there is a contradiction.
            for i in reversed(range(len(algGens)-1)):  # need test case
                alg_next, alg = algGens[i+1], algGens[i]
                if idem_id in alg.double_hor:
                    to_assign = self.NONE
                    if assignments[i+1] != self.NONE:
                        to_assign = assignments[i+1]
                    elif idem_id not in alg_next.double_hor:
                        # Now find the moving strand in alg_next starting at one
                        # of the points in the idempotent. Find whether the
                        # point is local or on the outside.
                        for s in [s for s, t in alg_next.strands
                                  if self.pmc2.pairid[s] == idem_id]:
                            if s in self.mapping2:
                                to_assign = self.LOCAL
                            else:
                                to_assign = self.OUTER
                    # Return no assignment if there is a contradiction
                    if assignments[i] != self.NONE and to_assign != self.NONE \
                            and assignments[i] != to_assign:
                        return None  # need test case
                    if to_assign != self.NONE:
                        assignments[i] = to_assign
            # There is one final case where everything is double horizontal.
            # Assign all to local (this agrees with how generators are created
            # for the extended bimodule, where the smeared pair always goes to
            # local).
            if all([idem_id in alg.double_hor for alg in algGens]):
                assignments = [self.LOCAL] * len(algGens)

            all_assignments.append(assignments)

        return all_assignments

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
            outer_sd = outer_sd.removeSingleHor(tuple(
                [self.single_idems_outer[single_id]
                 for single_id in range(self.num_singles)
                 if assignments[single_id][i] == self.LOCAL]))
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
        if prod_d is None:  # Case len(algGens) == 0.
            prod_d = \
                self.splitting2.restrictIdempotentOuter(MGen.idem2).toAlgElt()
            prod_d = prod_d.removeSingleHor()  # always goes to LOCAL

        # Test the inside
        local_MGen = MGen.local_gen
        alg_local = []
        for i in range(len(algGens)):
            cur_alg_local = \
                self.splitting2.restrictStrandDiagramLocal(algGens[i])
            cur_alg_local = cur_alg_local.removeSingleHor(tuple(
                [self.single_idems2[single_id]
                 for single_id in range(self.num_singles)
                 if assignments[single_id][i] == self.OUTER]))
            alg_local.append(cur_alg_local)
        alg_local = tuple(alg_local)

        # Assigning the smeared idempotents for the starting generator,
        # according to the rule that local_MGen.idem2 (A-side idempotent) must
        # match the left idempotent of the first algebra input (if there is
        # any).
        for i in range(self.num_singles):
            single2 = self.single_idems2[i]
            if single2 in local_MGen.idem2 and len(alg_local) > 0 and \
                    single2 not in alg_local[0].left_idem:
                if local_MGen not in self.local_da.u_maps[i]:
                    return E0  # need test case
                local_MGen = self.local_da.u_maps[i][local_MGen]
        local_delta = self.local_da.delta(local_MGen, alg_local)

        for (local_d, local_y), ring_coeff in local_delta.items():
            alg_d = self.splitting1.joinStrandDiagram(local_d, prod_d)
            if alg_d is None:
                continue
            outer_idem = prod_d.right_idem
            for i in range(self.num_singles):
                single1 = self.single_idems1[i]
                single_outer = self.single_idems_outer[i]
                if single_outer in prod_d.right_idem:
                    # If the split idempotent ended up on the outside, switch it
                    # to the inside.
                    if single1 not in local_y.idem1:
                        local_y = self.local_da.uinv_maps[i][local_y]
                    outer_idem = outer_idem.removeSingleHor([single_outer])
            y = self.gen_index[(local_y, outer_idem)]
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

        def testWithAssignment(assignments):
            # Test with a specific single assignments
            local_MGen = MGen.local_gen
            alg_local = []
            for i in range(len(algGens)):
                cur_alg_local = \
                    self.splitting2.restrictStrandDiagramLocal(algGens[i])
                cur_alg_local = cur_alg_local.removeSingleHor(tuple(
                    [self.single_idems2[single_id]
                     for single_id in range(self.num_singles)
                     if assignments[single_id][i] == self.OUTER]))
                alg_local.append(cur_alg_local)
            alg_local = tuple(alg_local)

            for i in range(self.num_singles):
                single2 = self.single_idems2[i]
                if single2 in local_MGen.idem2 and len(alg_local) > 0 and \
                        single2 not in alg_local[0].left_idem:
                    if local_MGen not in self.local_da.u_maps[i]:
                        return False  # need test case
                    local_MGen = self.local_da.u_maps[i][local_MGen]
            if not self.local_da.deltaPrefix(local_MGen, alg_local):
                return False

            # Test that the product on the outside is nonzero.
            has_product = True
            prod_d = None
            for i in range(len(algGens)):
                alg = algGens[i]
                outer_sd = self.splitting2.restrictStrandDiagramOuter(alg)
                outer_sd = outer_sd.removeSingleHor(tuple(
                    [self.single_idems_outer[single_id]
                     for single_id in range(self.num_singles)
                     if assignments[single_id][i] == self.LOCAL]))
                if prod_d is None:
                    prod_d = 1*outer_sd
                else:
                    prod_d = prod_d * outer_sd
                if prod_d == 0:
                    has_product = False
                    break
                else:
                    prod_d = prod_d.getElt()
            return has_product

        assignments = self.getSingleAssignments(algGens)
        if assignments is None:
            return False

        def fix_double_horizontal(assignments, single_id):
            """If assignments[single_id] is all self.LOCAL, return a list with
            two entries: the original assignment and a new assignment with all
            self.OUTER at single_id (both will need to be tested in
            deltaPrefix). Otherwise return a list with only the original
            assignment.

            """
            if all(a == self.LOCAL for a in assignments[single_id]):
                outer_assign = \
                    assignments[:single_id] + [[self.OUTER] * len(algGens)] + \
                    assignments[single_id+1:]
                return [assignments, outer_assign]
            else:
                return [assignments]

        all_assignments = [assignments]
        for i in range(self.num_singles):
            tmp = []
            for assignment in all_assignments:
                tmp.extend(fix_double_horizontal(assignment, i))
            all_assignments = tmp

        return any(testWithAssignment(assignment)
                   for assignment in all_assignments)

def identityDALocal(local_pmc):
    """Returns the identity type DA structure for a given local PMC.

    Actually the same as the non-local case, except we don't have Heegaard
    diagrams.

    """
    alg = local_pmc.getAlgebra()
    single_idems = local_pmc.getSingleIdems()
    dastr = LocalDAStructure(F2, alg, alg, single_idems1 = single_idems,
                             single_idems2 = single_idems)
    idems = local_pmc.getIdempotents()
    idem_to_gen_map = {}
    for i in range(len(idems)):
        cur_gen = SimpleDAGenerator(dastr, idems[i], idems[i], i)
        idem_to_gen_map[idems[i]] = cur_gen
        dastr.addGenerator(cur_gen)
    alg_gen = alg.getGenerators()
    for gen in alg_gen:
        if not gen.isIdempotent():
            gen_from = idem_to_gen_map[gen.getLeftIdem()]
            gen_to = idem_to_gen_map[gen.getRightIdem()]
            dastr.addDelta(gen_from, gen_to, gen, (gen,), 1)
    dastr.auto_u_map()
    return dastr
