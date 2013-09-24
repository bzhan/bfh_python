"""Defines type DD structures."""

from dstructure import *

class DDGenerator(Generator):
    """Represents a generator of type DD structure. Distinguished by (python)
    identity.

    """
    def __init__(self, parent, idem1, idem2):
        """Every generator has two idempotents (for the two type D actions)."""
        Generator.__init__(self, parent)
        self.idem1, self.idem2 = idem1, idem2

    def toSimpleDDGenerator(self, name):
        """Convert to a SimpleDDGenerator with the given name. All fields are
        preserved, except ``name`` which is overwritten, and _hash_val which is
        removed, if present.

        """
        new_obj = SimpleDDGenerator(self.parent, self.idem1, self.idem2, name)
        new_obj.__dict__.update(self.__dict__)
        new_obj.name = name # to make sure original name is overwritten
        if hasattr(new_obj, '_hash_val'):
            del new_obj._hash_val # reset hash value
        return new_obj

class SimpleDDGenerator(DDGenerator, NamedObject):
    """Represents a generator of type DD structure, distinguished by name."""
    def __init__(self, parent, idem1, idem2, name):
        """Specifies name in addition."""
        DDGenerator.__init__(self, parent, idem1, idem2)
        NamedObject.__init__(self, name)

    def __str__(self):
        return "%s:%s,%s" % (self.name, str(self.idem1), str(self.idem2))
        
    def __repr__(self):
        return str(self)

    def toDGenerator(self, new_parent):
        """Convert to a generator of a type D structure over the bialgebra."""
        new_idem = TensorIdempotent((self.idem1, self.idem2))
        new_gen = SimpleDGenerator(new_parent, new_idem, self.name)
        new_gen.__dict__.update(self.__dict__)
        new_gen.parent, new_gen.idem, new_gen.name = \
            new_parent, new_idem, self.name
        return new_gen

class MorDDtoDGenerator(DGenerator, MorObject):
    """Represents a generator of the type D structure of morphisms from a type
    DD structure to a type D structure.

    """
    def __init__(self, parent, source, coeff, target):
        """Specifies the morphism source -> coeff * target."""
        DGenerator.__init__(self, parent, source.idem2.opp())
        MorObject.__init__(self, source, coeff, target)

class DDStructure(FreeModule):
    """Represents a type DD structure. Note delta() returns an element in the
    tensor module Tensor((A,A,M)).

    """
    def __init__(self, ring, algebra1, algebra2, side1, side2):
        """Specifies the algebras and sides of the type DD action."""
        FreeModule.__init__(self, ring)
        assert isinstance(algebra1, DGAlgebra)
        assert isinstance(algebra2, DGAlgebra)
        self.algebra1 = algebra1
        self.side1 = side1
        self.algebra2 = algebra2
        self.side2 = side2

        # Construct A tensor A tensor M. Add diff and the left action on this
        # tensor product.
        self.AAtensorM = Tensor((algebra1, algebra2, self))
        def _mul_AA_AAtensorM((AGen1, AGen2, MGen), (ACoeff1, ACoeff2)):
            """To be used as rmultiply() in AAtensorM. Multiply ACoeff1 with
            AGen1 and ACoeff2 with AGen2.

            """
            return expandTensor((ACoeff1*AGen1, ACoeff2*AGen2, MGen),
                                 self.AAtensorM)

        def _diff_AAtensorM((AGen1, AGen2, MGen)):
            """To be used as diff() in AAtensorM."""
            return expandTensor((AGen1.diff(), AGen2, MGen), self.AAtensorM) \
                + expandTensor((AGen1, AGen2.diff(), MGen), self.AAtensorM) \
                + (AGen1, AGen2) * (MGen.delta())

        self.AAtensorM.rmultiply = _mul_AA_AAtensorM
        self.AAtensorM.diff = _diff_AAtensorM

    def delta(self, generator):
        """Returns delta^1 of the generator."""
        raise NotImplementedError("Differential not implemented.")

class SimpleDDStructure(DDStructure):
    """Represents a type DD structure with a finite number of generators, and
    explicitly stored generating set and delta operation.

    """
    def __init__(self, ring, algebra1, algebra2,
                 side1 = ACTION_LEFT, side2 = ACTION_LEFT):
        """Specifies the algebras and sides of the type D action."""
        assert side1 == ACTION_LEFT and side2 == ACTION_LEFT, \
            "Right action not implemented."
        DDStructure.__init__(self, ring, algebra1, algebra2, side1, side2)
        self.generators = set()
        self.delta_map = dict()

    def __len__(self):
        return len(self.generators)

    def delta(self, generator):
        return self.delta_map[generator]

    def getGenerators(self):
        return list(self.generators)

    def addGenerator(self, generator):
        """Add a generator. No effect if the generator already exists."""
        assert generator.parent == self
        assert isinstance(generator, DDGenerator)
        self.generators.add(generator)
        if not self.delta_map.has_key(generator):
            self.delta_map[generator] = E0

    def addDelta(self, gen_from, gen_to, alg_coeff1, alg_coeff2, ring_coeff):
        """Add ring_coeff * (alg_coeff1, alg_coeff2) * gen_to to the delta of
        gen_from. The first four arguments should be generators.

        """
        assert gen_from.parent == self and gen_to.parent == self
        assert alg_coeff1.getLeftIdem() == gen_from.idem1
        assert alg_coeff1.getRightIdem() == gen_to.idem1
        assert alg_coeff2.getLeftIdem() == gen_from.idem2
        assert alg_coeff2.getRightIdem() == gen_to.idem2
        target_gen = TensorGenerator((alg_coeff1, alg_coeff2, gen_to),
                                     self.AAtensorM)
        self.delta_map[gen_from] += target_gen.elt(ring_coeff)

    def deltaCoeff(self, gen_from, gen_to):
        """Return the coefficient (as bialgebra element) of gen_to in delta of
        gen_from.

        """
        if self.delta_map[gen_from] == 0:
            return E0
        else:
            # No need for algebra structure on the tensor product
            return self.delta_map[gen_from].fixLast(gen_to)

    def reindex(self):
        """Replace the generators by simple generators indexed by integers."""
        gen_list = list(self.generators)
        new_gen_list = []
        translate_dict = dict()
        for i in range(len(gen_list)):
            new_gen = gen_list[i].toSimpleDDGenerator("g%d"%(i+1))
            new_gen_list.append(new_gen)
            translate_dict[gen_list[i]] = new_gen
        self.generators = set(new_gen_list)
        new_delta = dict()
        for k, v in self.delta_map.items():
            new_v = E0
            for (AGen1, AGen2, MGen), coeff in v.items():
                target_gen = TensorGenerator(
                    (AGen1, AGen2, translate_dict[MGen]), self.AAtensorM)
                new_v += target_gen.elt(coeff)
            new_delta[translate_dict[k]] = new_v
        self.delta_map = new_delta
        if hasattr(self, "grading"):
            new_grading = dict()
            for gen, gr in self.grading.items():
                if gen in translate_dict: # gen is still in ddstr
                    new_grading[translate_dict[gen]] = gr
            self.grading = new_grading

    def testDelta(self):
        """Verify d^2 = 0 for this structure."""
        for gen in self.generators:
            if gen.delta().diff() != 0:
                return False
        return True

    def __str__(self):
        result = "Type DD Structure.\n"
        for k, v in self.delta_map.items():
            result += "d(%s) = %s\n" % (k, v)
        return result

    def morToD(self, other):
        """Compute the type D structure of morphisms from self to other. Note
        ``other`` must be a type D structure.

        """
        assert self.algebra1 == other.algebra
        alg_gens = self.algebra1.getGenerators()
        xlist = self.getGenerators()
        ylist = other.getGenerators()
        gens = list()
        dstr = SimpleDStructure(F2, self.algebra2.opp())
        genType = MorDDtoDGenerator

        def morGradingSet():
            """Find the grading set of the new type D structure."""
            lr_domains = [(d1, d2.opp())
                          for d1, d2 in self.gr_set.periodic_domains]
            self.lr_set = SimpleDbGradingSet(
                self.gr_set.gr_group1, ACTION_LEFT,
                self.gr_set.gr_group2.opp(), ACTION_RIGHT, lr_domains)
            return GeneralGradingSet([self.lr_set.inverse(), other.gr_set])

        def morGrading(gr_set, x, a, y):
            """Find the grading of the generator x -> ay in the morphism
            type D structure. The grading set need to be provided as gr_set.

            """
            gr_x1, gr_x2 = self.grading[x].data
            gr_x_lr = SimpleDbGradingSetElement(self.lr_set,
                                                (gr_x1, gr_x2.opp()))
            gr = [gr_x_lr.inverse(), other.grading[y] * a.getGrading()]
            return GeneralGradingSetElement(gr_set, gr)

        # Prepare rev_delta for the last step in computing differentials
        rev_delta = dict()
        for x in xlist:
            rev_delta[x] = []
        for p in xlist:
            for (b1, b2, q), coeff in p.delta().items():
                rev_delta[q].append(((b1, b2, p), coeff))

        # Get the list of generators
        for x in xlist:
            for a in alg_gens:
                for y in ylist:
                    if x.idem1 == a.getLeftIdem() and \
                            y.idem == a.getRightIdem():
                        gens.append(genType(dstr, x, a, y))
        for gen in gens:
            dstr.addGenerator(gen)

        # Get the type D structure maps
        for gen in gens:
            # Differential of y in (x -> ay)
            x, a, y = gen.source, gen.coeff, gen.target
            ady = a * y.delta()
            for (b, q), coeff in ady.items():
                dstr.addDelta(gen, genType(dstr, x, b, q), None, coeff)
            # Differential of a
            for da_gen, coeff in a.diff().items():
                dstr.addDelta(gen, genType(dstr, x, da_gen, y), None, coeff)
            # For each p such that (b1,b2)*x is in dp, add opp(b2)*(p->(b1*a)y)
            for (b1, b2, p), coeff1 in rev_delta[x]:
                for b1a_gen, coeff2 in (b1*a).items():
                    dstr.addDelta(gen, genType(dstr, p, b1a_gen, y),
                                  b2.opp(), coeff1*coeff2)

        # Find grading set and grading of elements
        if hasattr(self, "gr_set") and hasattr(other, "gr_set"):
            dstr.gr_set = morGradingSet()
            dstr.grading = dict()
            for gen in gens:
                dstr.grading[gen] = morGrading(
                    dstr.gr_set, gen.source, gen.coeff, gen.target)
        return dstr

    def simplify(self):
        """Simplify a type DD structure using cancellation lemma."""
        # Simplification is best done in terms of coefficients
        # Build dictionary of coefficients
        arrows = dict()
        for gen in self.generators:
            arrows[gen] = dict()
        bialgebra = TensorDGAlgebra((self.algebra1, self.algebra2))
        for gen in self.generators:
            for (AGen1, AGen2, MGen), coeff in self.delta_map[gen].items():
                if MGen not in arrows[gen]:
                    arrows[gen][MGen] = E0
                arrows[gen][MGen] += TensorGenerator(
                    (AGen1, AGen2), bialgebra) * coeff

        arrows = simplifyComplex(arrows, E0)

        # Now rebuild the type DD structure
        self.generators = set()
        self.delta_map = dict()
        for x in arrows:
            self.generators.add(x)
            self.delta_map[x] = E0
            for y, coeff in arrows[x].items():
                for (a1, a2), ring_coeff in coeff.items():
                    target_gen = TensorGenerator((a1, a2, y), self.AAtensorM)
                    self.delta_map[x] += ring_coeff * target_gen

    def toDStructure(self):
        """Convert this type DD structure into a type D structure over the
        tensor product of two algebras.

        """
        bialgebra = TensorDGAlgebra((self.algebra1, self.algebra2))
        dstr = SimpleDStructure(self.ring, bialgebra, ACTION_LEFT)
        gen_map = dict()
        for gen in self.generators:
            new_gen = gen.toDGenerator(dstr)
            gen_map[gen] = new_gen
            dstr.addGenerator(new_gen)
        for gen_from in self.generators:
            for (a1, a2, gen_to), coeff in self.delta_map[gen_from].items():
                dstr.addDelta(gen_map[gen_from], gen_map[gen_to],
                              TensorGenerator((a1, a2), bialgebra), coeff)
        return dstr

    def registerHDiagram(self, diagram, base_gen, base_gr = None):
        """Associate the given diagram as the Heegaard diagram from which this
        type DD structure can be derived. We will attempt to match generators
        of the type DD structure to generators of the Heegaard diagram.
        Currently this is possible only if no two generators have the same
        idempotents (so the match can be made by comparing idempotents).

        As a result, computes grading of each generator from the Heegaard
        diagram and checks it against type DD operations. Attributes added are:
        *. self.hdiagram - the Heegaard diagram.
        *. self.hdiagram_gen_map - dictionary mapping generators in the type DD
        structure to generators in Heegaard diagram.
        *. self.gr_set - the grading set (of type SimpleDbGradingSet).
        *. self.grading - dictionary mapping generators in the type DD
        structure to their gradings.

        """
        self.hdiagram = diagram
        # Get PMC's and check that they make sense
        hd_pmc1, hd_pmc2 = self.hdiagram.pmc_list
        dds_pmc1, dds_pmc2 = self.algebra1.pmc, self.algebra2.pmc
        assert hd_pmc1.opp() == dds_pmc1
        assert hd_pmc2.opp() == dds_pmc2
        # Now attempt to match generators
        self.hdiagram_gen_map = dict()
        gens, dgens = self.generators, diagram.getHFGenerators()
        for gen in gens:
            for dgen in dgens:
                dgen_idem1, dgen_idem2 = dgen.getDIdem()
                if (gen.idem1, gen.idem2) == (dgen_idem1, dgen_idem2):
                    self.hdiagram_gen_map[gen] = dgen
                    break
            assert gen in self.hdiagram_gen_map
        # Compute grading and check consistency with algebra actions
        base_hgen = self.hdiagram_gen_map[base_gen]
        self.gr_set, gr = self.hdiagram.computeDDGrading(base_hgen, base_gr)
        self.grading = dict()
        for gen in gens:
            self.grading[gen] = gr[self.hdiagram_gen_map[gen]]
        if ASSERT_LEVEL > 0:
            self.checkGrading()

    @memorize
    def dual(self):
        """Returns the dual of this type DD structure, which is the type DD
        invariant of the orientation reversed bordered 3-manifold. If self has
        left action over A1 and A2, then the new DD structure has left action
        over A2.opp() and A1.opp() (in that order).

        """
        dual_str = SimpleDDStructure(
            self.ring, self.algebra2.opp(), self.algebra1.opp(),
            self.side2, self.side1)
        # Map from generators in self to generators in dual_str
        gen_map = dict()
        for x in self.generators:
            # As in the type D case, simple generators only
            assert isinstance(x, SimpleDDGenerator)
            new_x = SimpleDDGenerator(dual_str, x.idem2.opp(), x.idem1.opp(),
                                      x.name)
            dual_str.addGenerator(new_x)
            gen_map[x] = new_x
        for x in self.generators:
            for (a, b, y), coeff in x.delta().items():
                dual_str.addDelta(gen_map[y], gen_map[x], b.opp(), a.opp(),
                                  coeff)
        return dual_str

    def checkGrading(self):
        for x in self.generators:
            for (a1, a2, y), coeff in x.delta().items():
                gr_x = self.grading[x]
                gr_y = self.grading[y]
                assert gr_x - 1 == gr_y * [a1.getGrading(), a2.getGrading()]

    def compareDDStructures(self, other):
        """Compare two type DD structures, print out any differences."""
        # Some basic tests:
        if len(self) != len(other):
            print "Different number of generators."""
            return False
        if (self.algebra1, self.algebra2) != (other.algebra1, other.algebra2):
            print "Different algebra action."
            return False
        gen_map = dict()
        for gen1 in self.generators:
            for gen2 in other.generators:
                if (gen1.idem1, gen1.idem2) == (gen2.idem1, gen2.idem2):
                    gen_map[gen1] = gen2
                    break
        for gen1 in self.generators:
            for gen2 in self.generators:
                coeff1 = self.deltaCoeff(gen1, gen2)
                coeff2 = other.deltaCoeff(gen_map[gen1], gen_map[gen2])
                if coeff1 != coeff2:
                    print "Different coefficient at %s->%s" % (gen1, gen2)
                    print "%s vs %s" % (coeff1, coeff2)
                    return False
        return True

def DDStrFromChords(alg1, alg2, idem_pairs, chord_pairs):
    """Construct type DD structure from list of idempotant pairs and chord
    pairs.

    """
    ddstr = SimpleDDStructure(F2, alg1, alg2)
    for i in range(len(idem_pairs)):
        ddstr.addGenerator(SimpleDDGenerator(
                ddstr, idem_pairs[i][0], idem_pairs[i][1], i))
    gen_set = ddstr.getGenerators()
    for x in gen_set:
        for y in gen_set:
            for l_chord, r_chord in chord_pairs:
                if alg1.mult_one and not l_chord.isMultOne():
                    continue
                if alg2.mult_one and not r_chord.isMultOne():
                    continue
                if l_chord.idemCompatible(x.idem1, y.idem1) and \
                        r_chord.idemCompatible(x.idem2, y.idem2):
                    ddstr.addDelta(x, y, StrandDiagram(alg1, x.idem1, l_chord),
                                   StrandDiagram(alg2, x.idem2, r_chord), 1)
    return ddstr

def identityDD(pmc):
    """Returns the identity type DD structure for a given PMC."""
    n = pmc.n
    pmcopp = pmc.opp()
    alg1 = pmc.getAlgebra()
    alg2 = pmcopp.getAlgebra()
    ddstr = SimpleDDStructure(F2, alg1, alg2)
    idems = pmc.getIdempotents()
    idem_pairs = [(idem, idem.opp().comp()) for idem in idems]
    chord_pairs = [(Strands(pmc, [(i, j)]), Strands(pmcopp, [(n-1-j, n-1-i)]))
                   for i in range(n) for j in range(i+1, n)]
    ddstr = DDStrFromChords(alg1, alg2, idem_pairs, chord_pairs)
    # Any generator can serve as base_gen
    for gen in ddstr.getGenerators():
        base_gen = gen
        break
    ddstr.registerHDiagram(getIdentityDiagram(pmc), base_gen)
    return ddstr

def DDStrFromDStr(dstr, genus1):
    """Obtain the type DD structure from a type D structure that is related by
    drilling. See section 7.3 of paper 'Computing HF by Factoring Mapping
    Classes'.

    If dstr has left type D action by algebra A(Z_1 # Z_2), where genus1
    specifies the genus of Z_1, then ddstr will have left type DD action by
    A(Z_1) and A(Z_2).

    """
    assert dstr.side == ACTION_LEFT
    pmc_all = dstr.algebra.pmc
    assert dstr.algebra.idem_size == pmc_all.genus
    pmc1, pmc2 = unconnectSumPMC(pmc_all, genus1)
    mult_one = dstr.algebra.mult_one
    
    ddstr = SimpleDDStructure(F2, pmc1.getAlgebra(None, mult_one),
                              pmc2.getAlgebra(None, mult_one))
    gen_map = {}
    for x in dstr.getGenerators():
        # Split idempotent of x into two parts
        xidem = x.idem
        x1_idem = Idempotent(pmc1,
                             [pairid for pairid in xidem if pairid < 2*genus1])
        x2_idem = Idempotent(pmc2,
                             [pairid-2*genus1 for pairid in xidem
                              if pairid >= 2*genus1])
        if len(x1_idem) != genus1:
            continue
        gen_map[x] = SimpleDDGenerator(ddstr, x1_idem, x2_idem, x.name)
        ddstr.addGenerator(gen_map[x])

    cut_point = 4 * genus1
    for x in dstr.getGenerators():
        for (a, y), coeff in x.delta().items():
            if a.multiplicity[cut_point-1] == 0:
                # The interval (cut_point-1, cut_point) is unoccupied
                a1, a2 = unconnectSumStrandDiagram(a, genus1)
                ddstr.addDelta(gen_map[x], gen_map[y], a1, a2, coeff)

    return ddstr
