"""Defines type D structures."""

from fractions import Fraction
from algebra import DGAlgebra, FreeModule, Generator, SimpleChainComplex, \
    Tensor, TensorGenerator
from algebra import simplifyComplex
from algebra import E0
from grading import GeneralGradingSet, GeneralGradingSetElement
from hdiagram import getZeroFrameDiagram, getInfFrameDiagram, getPlatDiagram
from pmc import Idempotent, Strands, StrandDiagram
from pmc import connectSumPMC, splitPMC, linearPMC
from utility import MorObject, NamedObject
from utility import memorize
from utility import ACTION_LEFT, DEFAULT_GRADING, F2, SMALL_GRADING

class DGenerator(Generator):
    """Represents a generator of type D structure. Distinguished by (python)
    identity.

    """
    def __init__(self, parent, idem):
        """Every generator must have an idempotent."""
        Generator.__init__(self, parent)
        self.idem = idem

    def toSimpleDGenerator(self, name):
        """Convert to a SimpleDGenerator with the given name. All fields are
        preserved, except ``name`` which is overwritten, and _hash_val which is
        removed, if present.

        """
        new_obj = SimpleDGenerator(self.parent, self.idem, name)
        new_obj.__dict__.update(self.__dict__)
        new_obj.name = name # to make sure original name is overwritten
        if hasattr(new_obj, '_hash_val'):
            del new_obj._hash_val # reset hash value
        return new_obj

class SimpleDGenerator(DGenerator, NamedObject):
    """Represents a generator of type D structure, distinguished by name."""
    def __init__(self, parent, idem, name):
        """Specifies name in addition."""
        DGenerator.__init__(self, parent, idem)
        NamedObject.__init__(self, name)

class MorDtoDGenerator(Generator, MorObject):
    """Represents a generator of the morphism complex from a type D structure
    to another type D structure.

    """
    def __init__(self, parent, source, coeff, target):
        """Specifies the morphism source -> coeff * target."""
        Generator.__init__(self, parent)
        MorObject.__init__(self, source, coeff, target)

class DStructure(FreeModule):
    """Represents a type D structure. Note delta() returns an element in the
    tensor module Tensor((A,M)).

    """
    def __init__(self, ring, algebra, side):
        """Specifies the algebra and side of the type D action."""
        FreeModule.__init__(self, ring)
        assert isinstance(algebra, DGAlgebra)
        self.algebra = algebra
        self.side = side

        # Construct A tensor M. Add diff and the left action of A on this
        # tensor product.
        self.AtensorM = Tensor((algebra, self))
        def _mul_A_AtensorM((AGen, MGen), ACoeff):
            """To be used as rmultiply() in AtensorM. Multiply ACoeff with
            AGen.

            """
            return (ACoeff * AGen) * MGen

        def _diff_AtensorM((AGen, MGen)):
            """To be used as diff() in AtensorM."""
            return (AGen.diff() * MGen) + (AGen * MGen.delta())

        self.AtensorM.rmultiply = _mul_A_AtensorM
        self.AtensorM.diff = _diff_AtensorM

    def delta(self, generator):
        """Returns delta^1 of the generator."""
        raise NotImplementedError("Differential not implemented.")

    def rmultiply(self, MGen, AGen):
        """Multiply a generator of the DStructure with an algebra generator
        means forming the tensor.

        """
        return 1*TensorGenerator((AGen, MGen), self.AtensorM)

class SimpleDStructure(DStructure):
    """Represents a type D structure with a finite number of generators, and
    explicitly stored generating set and delta operation.

    """
    def __init__(self, ring, algebra, side = ACTION_LEFT):
        """Initializes an empty type D structure."""
        assert side == ACTION_LEFT, "Right action not implemented."
        DStructure.__init__(self, ring, algebra, side)
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
        assert isinstance(generator, DGenerator)
        self.generators.add(generator)
        if not self.delta_map.has_key(generator):
            self.delta_map[generator] = E0

    def addDelta(self, gen_from, gen_to, alg_coeff, ring_coeff):
        """Add ring_coeff * alg_coeff * gen_to to the delta of gen_from. Both
        arguments should be generators.

        """
        assert gen_from.parent == self and gen_to.parent == self
        if alg_coeff is None:
            alg_coeff = gen_to.idem.toAlgElt(self.algebra)
        assert alg_coeff.getLeftIdem() == gen_from.idem
        assert alg_coeff.getRightIdem() == gen_to.idem
        self.delta_map[gen_from] += (alg_coeff * gen_to) * ring_coeff

    def reindex(self):
        """Replace the generators by simple generators indexed by integers."""
        gen_list = list(self.generators)
        new_gen_list = []
        translate_dict = dict()
        for i in range(len(gen_list)):
            new_gen = gen_list[i].toSimpleDGenerator("g%d"%(i+1))
            new_gen_list.append(new_gen)
            translate_dict[gen_list[i]] = new_gen
        self.generators = set(new_gen_list)
        new_delta = dict()
        for k, v in self.delta_map.items():
            new_v = E0
            for (AGen, MGen), coeff in v.items():
                new_v += (AGen * translate_dict[MGen]) * coeff
            new_delta[translate_dict[k]] = new_v
        self.delta_map = new_delta
        if hasattr(self, "grading"):
            new_grading = dict()
            for gen, gr in self.grading.items():
                if gen in translate_dict: # gen is still in dstr
                    new_grading[translate_dict[gen]] = gr
            self.grading = new_grading

    def deltaCoeff(self, gen_from, gen_to):
        """Return the coefficient (as algebra element) of gen_to in delta of
        gen_from.

        """
        if self.delta_map[gen_from] == 0:
            return E0
        else:
            return self.delta_map[gen_from].fixLast(gen_to)

    def testDelta(self):
        """Verify d^2 = 0 for this structure."""
        for gen in self.generators:
            if gen.delta().diff() != 0:
                # Print the offending terms in d^2 for one generator.
                print gen, "==>"
                for k, v in gen.delta().diff().items():
                    print v, "*", k
                return False
        return True

    def __str__(self):
        result = "Type D Structure.\n"
        for k, v in self.delta_map.items():
            result += "d(%s) = %s\n" % (k, v)
        return result

    def morToD(self, other):
        """Compute the chain complex of morphisms from self to other."""
        assert self.algebra == other.algebra
        alg_gens = self.algebra.getGenerators()
        xlist = self.getGenerators()
        ylist = other.getGenerators()
        gens = list()
        cx = SimpleChainComplex(F2)
        genType = MorDtoDGenerator

        def morGradingSet():
            """Find the grading set of the new chain complex."""
            return GeneralGradingSet([self.gr_set.inverse(), other.gr_set])

        def morGrading(gr_set, x, a, y):
            """Find the grading of the generator x -> ay in the morphism
            complex. The grading set need to be provided as gr_set.

            """
            gr = [self.grading[x].inverse(), other.grading[y] * a.getGrading()]
            return GeneralGradingSetElement(gr_set, gr)

        # Prepare rev_delta for the last step in computing differentials
        rev_delta = dict()
        for x in xlist:
            rev_delta[x] = []
        for p in xlist:
            for (b, q), coeff in p.delta().items():
                rev_delta[q].append(((b, p), coeff))

        # Get the list of generators
        for x in xlist:
            for a in alg_gens:
                for y in ylist:
                    if x.idem == a.getLeftIdem() and \
                            y.idem == a.getRightIdem():
                        gens.append(genType(cx, x, a, y))
        for gen in gens:
            cx.addGenerator(gen)

        # Get differentials
        for gen in gens:
            # Differential of ay in (x -> ay)
            x, a, y = gen.source, gen.coeff, gen.target
            day = a * y.delta() + a.diff() * y
            for (b, q), coeff in day.items():
                cx.addDifferential(gen, genType(cx, x, b, q), coeff)
            # For each p such that b*x is in dp, add p->(ba)y
            for (b, p), coeff1 in rev_delta[x]:
                for ba_gen, coeff2 in (b*a).items():
                    cx.addDifferential(
                        gen, genType(cx, p, ba_gen, y), coeff1*coeff2)

        # Find grading set and grading of elements
        if hasattr(self, "gr_set") and hasattr(other, "gr_set"):
            cx.gr_set = morGradingSet()
            cx.grading = dict()
            for gen in gens:
                cx.grading[gen] = morGrading(cx.gr_set,
                                             gen.source, gen.coeff, gen.target)
        return cx

    def simplify(self):
        """Simplify a type D structure using cancellation lemma."""
        # Simplification is best done in terms of coefficients
        # Build dictionary of coefficients
        arrows = dict()
        for gen in self.generators:
            arrows[gen] = dict()
        for gen in self.generators:
            for (AGen, MGen), coeff in self.delta_map[gen].items():
                if MGen not in arrows[gen]:
                    arrows[gen][MGen] = E0
                arrows[gen][MGen] += AGen * coeff

        arrows = simplifyComplex(arrows, E0)

        # Now rebuild the type D structure
        self.generators = set()
        self.delta_map = dict()
        for x in arrows:
            self.generators.add(x)
            self.delta_map[x] = E0
            for y, coeff in arrows[x].items():
                self.delta_map[x] += coeff * y

        # This is a good place to simplify gradings
        if hasattr(self, "gr_set"):
            new_gr_set = self.gr_set.simplifiedSet()
            for gen in self.generators:
                self.grading[gen] = self.gr_set.simplifiedElt(self.grading[gen])
            self.gr_set = new_gr_set

    def registerHDiagram(self, diagram, base_gen, base_gr = None):
        """Associate the given diagram as the Heegaard diagram from which this
        type D structure can be derived. Broadly similar (and somewhat simpler)
        than the type DD case. See the corresponding method for ddstructure for
        details.

        """
        self.hdiagram = diagram
        # Match PMC's and check that they make sense
        hd_pmc = self.hdiagram.pmc_list[0]
        dds_pmc = self.algebra.pmc
        assert hd_pmc.opp() == dds_pmc
        # Now attempt to match generators
        self.hdiagram_gen_map = dict()
        gens, dgens = self.generators, diagram.getHFGenerators()
        for gen in gens:
            for dgen in dgens:
                dgen_idem = dgen.getDIdem()[0]
                if gen.idem == dgen_idem:
                    self.hdiagram_gen_map[gen] = dgen
                    break
            assert gen in self.hdiagram_gen_map
        # Compute grading and check consistency with algebra actions
        base_hgen = self.hdiagram_gen_map[base_gen]
        self.gr_set, gr = self.hdiagram.computeDGrading(base_hgen, base_gr)
        self.grading = dict()
        for gen in gens:
            self.grading[gen] = gr[self.hdiagram_gen_map[gen]]
        self.checkGrading()

    @memorize
    def dual(self):
        """Returns the dual of this type D structure, which is the type D
        invariant of the orientation reversed bordered 3-manifold. Reverse all
        arrows and take the opp() of all coefficients. The result is a type D
        structure over the opposite algebra (acting from the same side).

        """
        dual_str = SimpleDStructure(self.ring, self.algebra.opp(), self.side)
        # Map from generators in self to generators in dual_str:
        gen_map = dict()
        for x in self.generators:
            # Don't want to deal with the case where x represents more
            # complicated information. Use reindex() to reduce to this case.
            assert isinstance(x, SimpleDGenerator)
            new_x = SimpleDGenerator(dual_str, x.idem.opp(), x.name)
            dual_str.addGenerator(new_x)
            gen_map[x] = new_x
        for x in self.generators:
            for (a, y), coeff in x.delta().items():
                dual_str.addDelta(gen_map[y], gen_map[x], a.opp(), coeff)
        if hasattr(self, "gr_set"):
            dual_str.gr_set = self.gr_set.inverse().opp()
            dual_str.grading = dict()
            for x in self.generators:
                dual_str.grading[gen_map[x]] = self.grading[x].inverse().opp()
        return dual_str

    def checkGrading(self):
        """Check grading is consistent with the type D operations."""
        for x in self.generators:
            for (a, y), coeff in x.delta().items():
                gr_x = self.grading[x]
                gr_y = self.grading[y]
                assert gr_x - 1 == gr_y * [a.getGrading()]

def connectSumTypeD(dstr1, dstr2):
    """Form the connect sum of two type D structures."""
    algebra1, algebra2 = dstr1.algebra, dstr2.algebra
    assert algebra1.mult_one == algebra2.mult_one
    pmc1, pmc2 = algebra1.pmc, algebra2.pmc
    pmc = connectSumPMC(pmc1, pmc2)
    algebra = pmc.getAlgebra(mult_one = algebra1.mult_one)
    dstr = SimpleDStructure(F2, algebra)
    # Maps pairs of generators in dstr1 and dstr2 to a generator in dstr, and
    # vice versa.
    pair_map = dict()
    rev_pair_map = dict()
    for gen1 in dstr1.getGenerators():
        for gen2 in dstr2.getGenerators():
            assert all([isinstance(x, SimpleDGenerator) for x in (gen1, gen2)])
            idem = list(gen1.idem) + [p+pmc1.num_pair for p in gen2.idem]
            idem = Idempotent(pmc, idem)
            gen = SimpleDGenerator(dstr, idem, gen1.name + gen2.name)
            dstr.addGenerator(gen)
            pair_map[(gen1, gen2)] = gen
            rev_pair_map[gen] = (gen1, gen2)
    for gen in dstr.getGenerators():
        gen1, gen2 = rev_pair_map[gen]
        for (a, y), coeff in gen1.delta().items():
            new_strands = Strands(pmc, a.strands)
            new_a = StrandDiagram(algebra, gen.idem, new_strands)
            dstr.addDelta(gen, pair_map[(y, gen2)], new_a, coeff)
        for (a, y), coeff in gen2.delta().items():
            new_strands = Strands(
                pmc, [(p+pmc1.n, q+pmc1.n) for p,q in a.strands])
            new_a = StrandDiagram(algebra, gen.idem, new_strands)
            dstr.addDelta(gen, pair_map[(gen1, y)], new_a, coeff)
    return dstr

typeDGrs1 = {"zeroDual" : (0, [0,0]),
             "zeroReg" : (Fraction(1,2), [0,Fraction(1,2)]),
             "infDual" : (Fraction(1,4), [0,0]),
             "infReg" : (-Fraction(1,4), [-Fraction(1,2),0])}

typeDGrs2 = {"zeroDual" : (0, [0,0]),
             "zeroReg" : (0, [0,-Fraction(1,2)]),
             "infDual" : (Fraction(1,4), [0,0]),
             "infReg" : (Fraction(1,4), [Fraction(1,2),0])}

typeDGrs3 = {"zeroDual" : (0, [0,-1]),
             "zeroReg" : (0, [0,Fraction(1,2)]),
             "infDual" : (-Fraction(3,4), [1,0]),
             "infReg" : (-Fraction(3,4), [-Fraction(1,2),0])}

typeDGrs4 = {"zeroDual" : (Fraction(1,2), [0,-1]),
             "zeroReg" : (0, [0,-Fraction(1,2)]),
             "infDual" : (-Fraction(1,4), [1,0]),
             "infReg" : (Fraction(1,4), [Fraction(1,2),0])}

typeDGrs5 = {"zeroDual" : (0, [0,Fraction(-1,2)]),
             "zeroReg" : (Fraction(1,2), [0,-1]),
             "infDual" : (Fraction(1,4), [Fraction(1,2),0]),
             "infReg" : (-Fraction(1,4), [1,0])}

typeDGrs6 = {"zeroDual" : (0, [0,Fraction(-1,2)]),
             "zeroReg" : (0, [0,0]),
             "infDual" : (Fraction(1,4), [Fraction(1,2),0]),
             "infReg" : (Fraction(1,4), [0,0])}

typeDGrs7 = {"zeroDual" : (0, [0,Fraction(1,2)]),
             "zeroReg" : (0, [0,-1]),
             "infDual" : (-Fraction(3,4), [Fraction(-1,2),0]),
             "infReg" : (-Fraction(3,4), [1,0])}

typeDGrs8 = {"zeroDual" : (Fraction(1,2), [0,Fraction(1,2)]),
             "zeroReg" : (0, [0,0]),
             "infDual" : (-Fraction(1,4), [Fraction(-1,2),0]),
             "infReg" : (Fraction(1,4), [0,0])}

typeDGrs = [typeDGrs1, typeDGrs2, typeDGrs3, typeDGrs4,
            typeDGrs5, typeDGrs6, typeDGrs7, typeDGrs8]

def getDGrs(abs_gr_info, code_str):
    """Returns the maslov and spinc components of the absolute grading using
    the given grading info and code string.

    """
    maslov, spinc = 0, []
    for info in abs_gr_info:
        cur_maslov, cur_spinc = typeDGrs[info][code_str]
        maslov += cur_maslov
        spinc += cur_spinc
    return maslov, spinc

def zeroTypeD(genus, is_dual = False, abs_gr_info = None):
    """Returns the type D structure for the 0-framed handlebody of a given
    genus.

    """
    pmc = splitPMC(genus)
    algebra = pmc.getAlgebra()
    dstr = SimpleDStructure(F2, algebra)
    idem = pmc.idem([4*i for i in range(genus)])
    genx = SimpleDGenerator(dstr, idem, "x")
    dstr.addGenerator(genx)
    for i in range(genus):
        sd = StrandDiagram(algebra, idem, [(4*i,4*i+2)])
        dstr.addDelta(genx, genx, sd, 1)
    if abs_gr_info is None:
        genx_gr = None
    else:
        assert DEFAULT_GRADING == SMALL_GRADING
        if is_dual:
            maslov, spinc = getDGrs(reversed(abs_gr_info), "zeroDual")
        else:
            maslov, spinc = getDGrs(abs_gr_info, "zeroReg")
        genx_gr = pmc.small_gr(maslov, spinc) # really pmc_opp
    dstr.registerHDiagram(getZeroFrameDiagram(genus), genx, genx_gr)
    if is_dual:
        dstr = dstr.dual()
    return dstr

def zeroTypeDAdm(genus):
    """Returns a larger type D structure for the 0-framed handlebody of a given
    genus. The diagram for this is obtained by isotopying the beta circles to
    create more intersections, so it is more likely to create admissible
    diagrams when tensored with another bordered diagram.

    """
    if genus > 1:
        return connectSumTypeD(zeroTypeDAdm(genus-1), zeroTypeDAdm(1))
    # genus == 1 case
    pmc = splitPMC(1)
    algebra = pmc.getAlgebra()
    dstr = SimpleDStructure(F2, algebra)
    idem_x = pmc.idem([0])
    idem_o = pmc.idem([1]) # idem for the other two generators
    genx = SimpleDGenerator(dstr, idem_x, "x")
    geny = SimpleDGenerator(dstr, idem_o, "y")
    genz = SimpleDGenerator(dstr, idem_o, "z")
    [dstr.addGenerator(gen) for gen in [genx, geny, genz]]
    dstr.addDelta(genz, geny, StrandDiagram(algebra, idem_o, []), 1)
    dstr.addDelta(genz, genx, StrandDiagram(algebra, idem_o, [(1,2)]), 1)
    dstr.addDelta(genx, geny, StrandDiagram(algebra, idem_x, [(0,1)]), 1)
    return dstr

def infTypeD(genus, is_dual = False, abs_gr_info = None):
    """Returns the type D structure for the inf-framed handlebody of a given
    genus.

    """
    pmc = splitPMC(genus)
    algebra = pmc.getAlgebra()
    dstr = SimpleDStructure(F2, algebra)
    idem = pmc.idem([4*i+1 for i in range(genus)])
    geny = SimpleDGenerator(dstr, idem, "y")
    dstr.addGenerator(geny)
    for i in range(genus):
        sd = StrandDiagram(algebra, idem, [(4*i+1, 4*i+3)])
        dstr.addDelta(geny, geny, sd, 1)
    if abs_gr_info is None:
        geny_gr = None
    else:
        assert DEFAULT_GRADING == SMALL_GRADING
        if is_dual:
            maslov, spinc = getDGrs(reversed(abs_gr_info), "infDual")
        else:
            maslov, spinc = getDGrs(abs_gr_info, "infReg")
        geny_gr = pmc.small_gr(maslov, spinc) # really pmc_opp
    dstr.registerHDiagram(getInfFrameDiagram(genus), geny, geny_gr)
    if is_dual:
        dstr = dstr.dual()
    return dstr

def platTypeD(genus):
    """Returns the type D structure for the plat handlebody of a given
    genus.

    """
    pmc = linearPMC(genus)
    algebra = pmc.getAlgebra()
    dstr = SimpleDStructure(F2, algebra)
    idem = pmc.idem([4*i+1 for i in range(genus-1)]+[4*genus-3])
    genx = SimpleDGenerator(dstr, idem, "x")
    dstr.addGenerator(genx)
    strands = [(4*i+1,4*i+4) for i in range(genus-1)]+[(4*genus-3, 4*genus-1)]
    for st in strands:
        sd = StrandDiagram(algebra, idem, [st])
        dstr.addDelta(genx, genx, sd, 1)
    dstr.registerHDiagram(getPlatDiagram(genus), genx)
    return dstr
