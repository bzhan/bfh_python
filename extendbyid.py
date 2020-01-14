"""Extension by identity of type DA structures."""

from algebra import TensorGenerator
from algebra import E0
from dastructure import DAGenerator, DAStructure, DATensorDGenerator, \
    MorDAtoDAComplex, SimpleDAGenerator, SimpleDAStructure
from dstructure import SimpleDStructure
from grading import GeneralGradingSet, GeneralGradingSetElement
from localpmc import LocalStrandAlgebra, PMCSplitting
from utility import subset
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

    def delta(self, MGen, algGens):
        if len(algGens) == 1 and algGens[0].isIdempotent() and \
           algGens[0].left_idem == MGen.idem2:
               return MGen.idem1.toAlgElt() * MGen
        elif (MGen, algGens) not in self.da_action:
            return E0
        else:
            return self.da_action[(MGen, algGens)]

                    
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

        for (x1, coeffs_a), target in list(self.source.da_action.items()):
            for (coeff_d, x2), ring_coeff in list(target.items()):
                result.addDelta(
                    gen_map[x1], gen_map[x2], coeff_d, coeffs_a, ring_coeff)
        for (y1, coeffs_a), target in list(self.target.da_action.items()):
            for (coeff_d, y2), ring_coeff in list(target.items()):
                result.addDelta(
                    gen_map[y1], gen_map[y2], coeff_d, coeffs_a, ring_coeff)
        for gen, ring_coeff in list(morphism.items()):
            # coeffs_a is a tuple of A-side inputs
            coeff_d, coeffs_a = gen.coeff
            result.addDelta(gen_map[gen.source], gen_map[gen.target],
                            coeff_d, tuple(coeffs_a), ring_coeff)

        # Set up u_map
        num_single_idems = len(self.source.single_idems1)
        for idem_id in range(num_single_idems):
            for x, u_x in list(self.source.u_maps[idem_id].items()):
                result.add_u_map(idem_id, gen_map[x], gen_map[u_x])
            for y, u_y in list(self.target.u_maps[idem_id].items()):
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

        # Possible values of single assignments, for use in tensorD, delta and
        # deltaPrefix (through the function getSingleAssignments).
        self.NONE, self.LOCAL, self.OUTER, self.DOUBLE = 0, 1, 2, 3

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

    def adjustLocalMGen(self, local_MGen, alg_local0):
        """Assigning the smeared idempotents for the starting generator,
        according to the rule that local_MGen.idem2 (A-side idempotent) must
        match the left idempotent of the first algebra input (if there is any).

        """
        for i in range(self.num_singles):
            single2 = self.single_idems2[i]
            if single2 in local_MGen.idem2 and \
               single2 not in alg_local0.left_idem:
                if local_MGen not in self.local_da.u_maps[i]:
                    return None  # need test case
                local_MGen = self.local_da.u_maps[i][local_MGen]
        return local_MGen

    def testPrefix(self, local_MGen, algs_local):
        """Query deltaPrefix for the given set of local algebra inputs. Perform
        the adjustment on local_MGen if necessary.

        """
        if len(algs_local) > 0:
            local_MGen = self.adjustLocalMGen(local_MGen, algs_local[0])
            if local_MGen is None:
                return False
        return self.local_da.deltaPrefix(local_MGen, tuple(algs_local))

    def extendRestrictions(self, last_assign, algs_local, prod_d, new_alg):
        """Update the idempotent assignments when a new algebra input, new_alg,
        is considered. Apply possible changes to previous idempotent assignments
        to both algs_local and prod_d. Adds the local restriction of new_alg to
        algs_local, but does NOT multiply outer restriction of new_alg to
        prod_d (this is done in a separate function getNewProdD for efficiency
        considerations.

        """
        # Update single assignments
        new_assign = []
        for i in range(self.num_singles):
            idem_id = self.smeared_idems2[i]
            if idem_id in new_alg.double_hor:
                # Double horizontal in the new algebra element. Just continue
                # the previous assignment.
                assert last_assign[i] != self.NONE
                new_assign.append(last_assign[i])
            else:
                # First determine new assignment.
                if idem_id in new_alg.right_idem:
                    end_pt = [t for s, t in new_alg.strands
                              if self.pmc2.pairid[t] == idem_id]
                    assert len(end_pt) == 1
                    end_pt = end_pt[0]
                    if end_pt in self.mapping2:
                        new_assign.append(self.LOCAL)
                    else:
                        new_assign.append(self.OUTER)
                else:
                    new_assign.append(self.NONE)

                # Now correct previous assignment if necessary.
                if idem_id in new_alg.left_idem:
                    assert last_assign[i] != self.NONE
                    start_pt = [s for s, t in new_alg.strands
                                if self.pmc2.pairid[s] == idem_id]
                    assert len(start_pt) == 1
                    start_pt = start_pt[0]
                    if start_pt in self.mapping2:
                        if last_assign[i] == self.OUTER:
                            return (None, None, None)  # conflict
                    else:
                        if last_assign[i] == self.LOCAL:
                            return (None, None, None)  # conflict
                        elif last_assign[i] == self.DOUBLE:
                            # Previous assignment changes from DOUBLE to local.
                            # Need to update algs_local and prod_d.
                            to_remove = (self.single_idems2[i],)
                            algs_local = [alg.removeSingleHor(to_remove)
                                          for alg in algs_local]
                            to_add = (self.single_idems_outer[i],)
                            prod_d = prod_d.addSingleHor(to_add)

        # Restrict current algebra element to local and form new_local.
        new_local = [alg for alg in algs_local]
        cur_alg_local = self.splitting2.restrictStrandDiagramLocal(new_alg)
        idems_to_remove = [self.single_idems2[single_id]
                           for single_id in range(self.num_singles)
                           if new_assign[single_id] == self.OUTER]
        cur_alg_local = cur_alg_local.removeSingleHor(tuple(idems_to_remove))
        if len(new_local) != 0:
            assert new_local[-1].right_idem == cur_alg_local.left_idem
        new_local.append(cur_alg_local)

        return (new_assign, new_local, prod_d)

    def getNewProdD(self, new_assign, new_alg, last_prod_d):
        """Multiplies the outer restriction of new_alg onto last_prod_d."""
        outer_sd = self.splitting2.restrictStrandDiagramOuter(new_alg)
        outer_sd = outer_sd.removeSingleHor(tuple(
            [self.single_idems_outer[single_id]
             for single_id in range(self.num_singles)
             if new_assign[single_id] in (self.LOCAL, self.DOUBLE)]))

        assert last_prod_d.right_idem == outer_sd.left_idem
        new_prod_d = last_prod_d * outer_sd
        if new_prod_d == 0:
            return None
        else:
            return new_prod_d.getElt()

    def getAssignments(self, MGen, algs):
        """Returns the triple (assignment, algs_local, prod_d)."""
        assignment = [self.DOUBLE] * self.num_singles
        algs_local = []
        prod_d = self.splitting2.restrictIdempotentOuter(MGen.idem2).toAlgElt()
        prod_d = prod_d.removeSingleHor()

        for alg in algs:
            assignment, algs_local, prod_d = self.extendRestrictions(
                assignment, algs_local, prod_d, alg)
            if assignment is None:
                return (None, None, None)
            prod_d = self.getNewProdD(assignment, alg, prod_d)
            if prod_d is None:
                return (None, None, None)
        return (assignment, algs_local, prod_d)

    def joinOutput(self, local_d, local_y, outer_d):
        """Joins local_d and outer_d. Adjust idempotents if necessary."""
        alg_d = self.splitting1.joinStrandDiagram(local_d, outer_d)
        if alg_d is None:
            return (None, None)
        outer_idem = outer_d.right_idem
        for i in range(self.num_singles):
            single1 = self.single_idems1[i]
            single_outer = self.single_idems_outer[i]
            if single_outer in outer_idem:
                # If the split idempotent ended up on the outside, switch it to
                # the inside.
                if single1 not in local_y.idem1:
                    local_y = self.local_da.uinv_maps[i][local_y]
                outer_idem = outer_idem.removeSingleHor([single_outer])
        y = self.gen_index[(local_y, outer_idem)]
        return (alg_d, y)

    def tensorD(self, dstr):
        """Compute the box tensor product DA * D of this bimodule with the given
        type D structure. Returns the resulting type D structure. Uses delta()
        and deltaPrefix() functions of this type DA structure.

        """
        dstr_result = SimpleDStructure(F2, self.algebra1)
        # Compute list of generators in the box tensor product
        for gen_left in self.getGenerators():
            for gen_right in dstr.getGenerators():
                if gen_left.idem2 == gen_right.idem:
                    dstr_result.addGenerator(DATensorDGenerator(
                        dstr_result, gen_left, gen_right))

        def search(start_gen, cur_dgen, algs, last_assign, algs_local,
                   last_prod_d):
            """Searching for an arrow in the box tensor product.
            - start_gen: starting generator in the box tensor product. The
              resulting arrow will start from here.
            - cur_dgen: current location in the type D structure.
            - algs: current list of A-side inputs to the type DA structure (or
              alternatively, list of algebra outputs produced by the existing
              path through the type D structure).
            - algs_local: current list of local restrictions of algs.
            - last_assign: a list of length self.num_singles. For each split
              idempotent, specify the single assignments at the last algebra
              input.
            - prod_d: product of the outer restrictions, except for the last
              algebra input.

            """
            start_dagen, start_dgen = start_gen
            local_MGen = start_dagen.local_gen

            # Preliminary tests
            if len(algs) > 0:
                assert algs[0].left_idem == start_dagen.idem2
            for i in range(len(algs)-1):
                assert algs[i].right_idem == algs[i+1].left_idem
            if any(alg.isIdempotent() for alg in algs):
                return

            # First, adjust local module generator, and check for delta.
            if len(algs_local) > 0:
                local_MGen = self.adjustLocalMGen(local_MGen, algs_local[0])
                if local_MGen is None:
                    return
            local_delta = self.local_da.delta(local_MGen, tuple(algs_local))
            has_delta = (local_delta != E0)

            # Second, check for delta prefix.
            has_delta_prefix = False
            if len(algs) == 0:
                has_delta_prefix = True
            else:
                dbls = [self.single_idems2[i] for i in range(self.num_singles)
                        if last_assign[i] == self.DOUBLE]
                for to_remove in subset(dbls):
                    if len(to_remove) != 0:
                        cur_algs_local = tuple([alg.removeSingleHor(to_remove)
                                                for alg in algs_local])
                    else:
                        cur_algs_local = algs_local
                    if self.testPrefix(local_MGen, cur_algs_local):
                        has_delta_prefix = True
                        break

            if (not has_delta) and (not has_delta_prefix):
                return

            # Now, compute new prod_d.
            if len(algs) > 0:
                prod_d = self.getNewProdD(last_assign, algs[-1], last_prod_d)
            else:
                prod_d = last_prod_d
            if prod_d is None:
                return

            # If has_delta is True, add to delta
            for (local_d, local_y), ring_coeff in list(local_delta.items()):
                alg_d, y = self.joinOutput(local_d, local_y, prod_d)
                if alg_d is not None:
                    dstr_result.addDelta(start_gen, DATensorDGenerator(
                        dstr_result, y, cur_dgen), alg_d, 1)

            if not has_delta_prefix:
                return
            for (new_alg, dgen_to), ring_coeff in list(dstr.delta(cur_dgen).items()):
                new_assign, new_local, last_prod_d = self.extendRestrictions(
                    last_assign, algs_local, prod_d, new_alg)
                if new_assign is not None:
                    search(start_gen, dgen_to, algs + [new_alg],
                           new_assign, new_local, last_prod_d)

        # Perform search for each generator in dstr_result.
        for x in dstr_result.getGenerators():
            dagen, dgen = x
            prod_d = \
                self.splitting2.restrictIdempotentOuter(dagen.idem2).toAlgElt()
            prod_d = prod_d.removeSingleHor()  # always goes to LOCAL
            search(x, dgen, [], [self.DOUBLE] * self.num_singles, [], prod_d)
            # Add arrows coming from idempotent output on the D-side
            for (coeff_out, dgen_to), ring_coeff in list(dstr.delta(dgen).items()):
                if coeff_out.isIdempotent():
                    dstr_result.addDelta(
                        x, DATensorDGenerator(dstr_result, dagen, dgen_to),
                        dagen.idem1.toAlgElt(self.algebra1), 1)

        # Find grading set if available on both components
        def tensorGradingSet():
            """Find the grading set of the new type D structure."""
            return GeneralGradingSet([self.gr_set, dstr.gr_set])

        def tensorGrading(gr_set, dagen, dgen):
            """Find the grading of the generator (x, y) in the tensor type D
            structure. The grading set need to be provided as gr_set.

            """
            return GeneralGradingSetElement(
                gr_set, [self.grading[dagen], dstr.grading[dgen]])

        if hasattr(self, "gr_set") and hasattr(dstr, "gr_set"):
            dstr_result.gr_set = tensorGradingSet()
            dstr_result.grading = dict()
            for x in dstr_result.getGenerators():
                dagen, dgen = x
                dstr_result.grading[x] = tensorGrading(
                    dstr_result.gr_set, dagen, dgen)

        return dstr_result

    def delta(self, MGen, algGens):
        # Preliminary tests
        if len(algGens) > 0 and algGens[0].left_idem != MGen.idem2:
            return E0
        if any([algGens[i].right_idem != algGens[i+1].left_idem
                for i in range(len(algGens)-1)]):
            return E0
        if any([alg.isIdempotent() for alg in algGens]):
            return E0

        assignment, algs_local, prod_d = self.getAssignments(MGen, algGens)
        if assignment is None:
            return E0

        local_MGen = MGen.local_gen
        if len(algs_local) > 0:
            local_MGen = self.adjustLocalMGen(local_MGen, algs_local[0])
            if local_MGen is None:
                return E0

        local_delta = self.local_da.delta(local_MGen, tuple(algs_local))
        if local_delta == 0:
            return E0

        result = E0
        for (local_d, local_y), ring_coeff in list(local_delta.items()):
            alg_d, y = self.joinOutput(local_d, local_y, prod_d)
            result += 1 * TensorGenerator((alg_d, y), self.AtensorM)
        return result

    def deltaPrefix(self, MGen, algGens):
        # Preliminary tests
        if len(algGens) == 0:
            return True
        if algGens[0].left_idem != MGen.idem2:
            return False
        if any([algGens[i].right_idem != algGens[i+1].left_idem
                for i in range(len(algGens)-1)]):
            return False

        assignment, algs_local, prod_d = self.getAssignments(MGen, algGens)
        if assignment is None:
            return E0

        local_MGen = MGen.local_gen
        dbls = [self.single_idems2[i] for i in range(self.num_singles)
                if assignment[i] == self.DOUBLE]
        for to_remove in subset(dbls):
            if len(to_remove) != 0:
                cur_algs_local = tuple([alg.removeSingleHor(to_remove)
                                        for alg in algs_local])
            else:
                cur_algs_local = algs_local
            if self.testPrefix(local_MGen, cur_algs_local):
                return True
        return False

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
