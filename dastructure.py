"""Defines type DA structures."""

from algebra import CobarAlgebra, DGAlgebra, FreeModule, Generator, Tensor, \
    TensorGenerator, TensorStarGenerator
from algebra import ChainComplex, E0, TensorDGAlgebra
from dstructure import DGenerator, SimpleDStructure
from ddstructure import DDGenerator, SimpleDDGenerator, SimpleDDStructure
from grading import GeneralGradingSet, GeneralGradingSetElement, \
    SimpleDbGradingSetElement
from hdiagram import getIdentityDiagram
from pmc import StrandDiagram
from utility import MorObject, NamedObject
from utility import sumColumns
from utility import ACTION_LEFT, ACTION_RIGHT, F2

class DAGenerator(Generator):
    """Represents a generator of type DA structure. Distinguished by (python)
    identity.

    """
    def __init__(self, parent, idem1, idem2):
        """Every generator has two idempotents. idem1 is the type D idempotent
        on the left. idem2 is the type A idempotent on the right.

        """
        Generator.__init__(self, parent)
        self.idem1, self.idem2 = idem1, idem2

class SimpleDAGenerator(DAGenerator, NamedObject):
    """Represents a generator of type DA structure, distinguished by name."""
    def __init__(self, parent, idem1, idem2, name):
        """Specifies name in addition."""
        DAGenerator.__init__(self, parent, idem1, idem2)
        NamedObject.__init__(self, name)

    def __str__(self):
        return "%s:%s,%s" % (self.name, self.idem1, self.idem2)

    def __repr__(self):
        return str(self)

class DATensorDGenerator(DGenerator, tuple):
    """Generator of a type D structure formed by tensoring a type DA structure
    and a type D structure. Also serves as the generator of the type D structure
    formed by tensoring DD * CFAA(Id) * D, with the generator of CFAA(Id)
    implicit from the idempotents.

    gen_left is a generator of either a type DA structure (DA * D
    interpretation) or a type DD structure (DD * CFAA(Id) * D interpretation).

    gen_right is a generator of a type D structure.

    """
    def __new__(cls, parent, gen_left, gen_right):
        return tuple.__new__(cls, (gen_left, gen_right))

    def __init__(self, parent, gen_left, gen_right):
        """Specify generators on two sides of the tensor (DA/DD and D
        generators).

        """
        # Note tuple initialization is automatic
        DGenerator.__init__(self, parent, gen_left.idem1)
        self.gen_left = gen_left
        self.gen_right = gen_right
        filt = []
        if hasattr(gen_left, "filtration"):
            filt += gen_left.filtration
        if hasattr(gen_right, "filtration"):
            filt += gen_right.filtration
        if filt != []:
            self.filtration = filt

class DATensorDDGenerator(DDGenerator, tuple):
    """Generator of a type DD structure formed by tensoring a type DA structure
    and a type DD structure. Also serves as the generator of the type DD
    structure formed by tensoring DD * CFAA(Id) * DD, with the generator of
    CFAA(Id) implicit from the idempotents.

    gen_left is a generator of either a type DA structure (DA * DD
    interpretation) or a type DD structure (DA * CFAA(Id) * DD interpretation).

    gen_right is a generator of a type DD structure.

    """
    def __new__(cls, parent, gen_left, gen_right):
        return tuple.__new__(cls, (gen_left, gen_right))

    def __init__(self, parent, gen_left, gen_right):
        """Specify generators on two sides of the tensor (DA/DD and DD
        generators).

        """
        # Note tuple initialization is automatic
        DDGenerator.__init__(self, parent, gen_left.idem1, gen_right.idem2)
        self.gen_left = gen_left
        self.gen_right = gen_right
        filt = []
        if hasattr(gen_left, "filtration"):
            filt += gen_left.filtration
        if hasattr(gen_right, "filtration"):
            filt += gen_right.filtration
        if filt != []:
            self.filtration = filt

class MorDAtoDAGenerator(Generator, MorObject):
    """Represents a generator of the chain complex of bimodule morphisms from a
    type DA structure to a type DA structure.

    """
    def __init__(self, parent, coeff_d, coeffs_a, source, target):
        """Specifies the morphism m(source, coeffs_a) -> coeff_d * target.
        source and target are generators in two type DA bimodules with same
        algebra actions. If the bimodules have left type D action by algebra1
        and right type A action by algebra2, then as a MorObject coeff is of
        type TensorDGAlgebra(algebra1, CobarAlgebra(algebra2)).

        """
        Generator.__init__(self, parent)
        self.coeff_d, self.coeffs_a = coeff_d, coeffs_a
        cobar_alg = CobarAlgebra(source.parent.algebra2)
        tensor_alg = TensorDGAlgebra((source.parent.algebra1, cobar_alg))
        coeff = TensorGenerator(
            (coeff_d, TensorStarGenerator(coeffs_a, cobar_alg, source.idem2)),
            tensor_alg)
        MorObject.__init__(self, source, coeff, target)

    def __str__(self):
        coeff_d, coeffs_a = self.coeff
        return "m(%s:%s; %s) = %s*%s:%s" % \
            (self.source.name, self.source, coeffs_a, coeff_d,
             self.target.name, self.target)

class MorDAtoDAComplex(ChainComplex):
    """Represents the complex of type DA morphisms between two type DA
    structures.

    """
    def __init__(self, ring, source, target):
        """Specifies the source and target DA structures."""
        ChainComplex.__init__(self, ring)
        assert source.algebra1 == target.algebra1 and \
            source.algebra2 == target.algebra2
        assert source.side1 == target.side1 and source.side2 == target.side2
        self.source = source
        self.target = target

    def __eq__(self, other):
        # Unlike other structures, MorDDtoDDComplex is distinguished by its
        # source and target
        return self.source == other.source and self.target == other.target

    def __ne__(self, other):
        return not (self == other)

    def __hash__(self, other):
        return hash(tuple((self.source, self.target)))

    def diff(self, gen):
        result = E0
        x, c, y = gen.source, gen.coeff, gen.target
        c_d, cs_a = gen.coeff  # D-side output and list of A-side inputs

        # Differential of coefficient
        for dc, ring_coeff in c.diff().items():
            coeff_d, coeffs_a = dc
            result += ring_coeff * MorDAtoDAGenerator(
                self, coeff_d, coeffs_a, x, y)

        # Pre-compose with differential in source
        for (x1, coeffs_a), target in self.source.da_action.items():
            for (coeff_d, x2), ring_coeff in target.items():
                if x == x2 and coeff_d * c_d != E0:
                    result += 1*MorDAtoDAGenerator(
                        self, (coeff_d * c_d).getElt(), coeffs_a + cs_a, x1, y)

        # Post-compose with differential in target
        for (y1, coeffs_a), target in self.target.da_action.items():
            for (coeff_d, y2), ring_coeff in target.items():
                if y == y1 and c_d * coeff_d != E0:
                    result += 1*MorDAtoDAGenerator(
                        self, (c_d * coeff_d).getElt(), cs_a + coeffs_a, x, y2)
        return result

    def getMappingCone(self, morphism):
        """Returns the mapping cone of a morphism. This is broadly similar to
        that for DDStructures.

        """
        result = SimpleDAStructure(
            F2, self.source.algebra1, self.source.algebra2,
            self.source.side1, self.source.side2)
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
        return result

class DAStructure(FreeModule):
    """Represents a type DA structure. delta() takes a generator and a sequence
    of algebra generators (as a generator of the tensor algebra), and returns
    an element in the tensor module Tensor((A,M)), where A is algebra1 (D-side
    algebra).

    """
    def __init__(self, ring, algebra1, algebra2, side1, side2):
        """Specifies the algebras and sides of the type DA action."""
        FreeModule.__init__(self, ring)
        assert isinstance(algebra1, DGAlgebra)
        assert isinstance(algebra2, DGAlgebra)
        self.algebra1 = algebra1
        self.side1 = side1
        self.algebra2 = algebra2
        self.side2 = side2

        # Construct A tensor M. Add diff and the left action of A on this
        # tensor product.
        self.AtensorM = Tensor((algebra1, self))
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
            for i in range(len(coeffs_a)-1):
                if coeffs_a[i].right_idem != coeffs_a[i+1].left_idem:
                    return False
            return x.idem2 == coeffs_a[0].left_idem and \
                y.idem2 == coeffs_a[-1].right_idem

    def delta(self, MGen, algGens):
        """algGens = (a_1, ..., a_n) is an element of TensorAlgebra. Evaluates
        the type DA operation delta^1(MGen; a_1, ..., a_n).

        """
        raise NotImplementedError("Differential not implemented.")

    def deltaPrefix(self, MGen, algGens):
        """algGens = (a_1, ..., a_n) is an element of TensorAlgebra. Returns a
        boolean value indicating whether there exists an arrow in the type DA
        action with starting generator MGen, and whose list of algebra
        generators has a_1, ..., a_n as a *strict* prefix.

        """
        raise NotImplementedError("Prefix differential not implemented.")

    def rmultiply(self, MGen, AGen):
        """Multiply a generator of type DAStructure with an algebra generator
        means forming the tensor.

        """
        return 1*TensorGenerator((AGen, MGen), self.AtensorM)

    def toSimpleDAStructure(self):
        """Using delta and deltaPrefix, product a simple DA structure (with
        explicit DA action). Does not work when there are infinitely many
        actions.

        """
        assert self.side1 == ACTION_LEFT and self.side2 == ACTION_RIGHT
        dastr = SimpleDAStructure(F2, self.algebra1, self.algebra2,
                                  side1 = ACTION_LEFT, side2 = ACTION_RIGHT)
        gen_map = dict()
        for gen in self.getGenerators():
            gen_map[gen] = SimpleDAGenerator(
                dastr, gen.idem1, gen.idem2, gen.name)
            dastr.addGenerator(gen_map[gen])

        alg2_gens = [alg_gen for alg_gen in self.algebra2.getGenerators()
                     if not alg_gen.isIdempotent()]

        def search(start_gen, cur_coeffs_a):
            """Search for terms in the action, starting from the generator
            start_gen, and with the current list of A-side inputs cur_coeffs_a
            (to be possibly extended).

            """
            cur_delta = self.delta(start_gen, cur_coeffs_a)
            for (coeff_d, gen_to), ring_coeff in cur_delta.items():
                dastr.addDelta(gen_map[start_gen], gen_map[gen_to], coeff_d,
                               cur_coeffs_a, ring_coeff)
            if self.deltaPrefix(start_gen, cur_coeffs_a):
                for coeff_a in alg2_gens:
                    search(start_gen, cur_coeffs_a + (coeff_a,))

        for gen in self.getGenerators():
            search(gen, ())

        # Copy over Heegaard diagram and grading information
        if hasattr(self, "hdiagram"):
            dastr.hdiagram = self.hdiagram
            dastr.hdiagram_gen_map = dict()
            for dagen, hgen in self.hdiagram_gen_map.items():
                dastr.hdiagram_gen_map[gen_map[dagen]] = hgen
        if hasattr(self, "gr_set"):
            dastr.gr_set = self.gr_set
            dastr.grading = dict()
            for gen, gr in self.grading.items():
                dastr.grading[gen_map[gen]] = gr

        return dastr

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

        def search(start_gen, cur_dgen, cur_coeffs_a):
            """Searching for an arrow in the box tensor product.
            - start_gen: starting generator in the box tensor product. The
              resulting arrow will start from here.
            - cur_dgen: current location in the type D structure.
            - cur_coeffs_a: current list of A-side inputs to the type DA
              structure (or alternatively, list of algebra outputs produced by
              the existing path through the type D structure).

            """
            start_dagen, start_dgen = start_gen
            cur_delta = self.delta(start_dagen, cur_coeffs_a)
            for (coeff_d, gen_to), ring_coeff in cur_delta.items():
                dstr_result.addDelta(start_gen, DATensorDGenerator(
                    dstr_result, gen_to, cur_dgen), coeff_d, 1)
            if self.deltaPrefix(start_dagen, cur_coeffs_a):
                for (coeff_out, dgen_to), ring_coeff in \
                    dstr.delta(cur_dgen).items():
                    search(start_gen, dgen_to, cur_coeffs_a + (coeff_out,))

        for x in dstr_result.getGenerators():
            dagen, dgen = x
            search(x, dgen, ())
            # Add arrows coming from idempotent output on the D-side
            for (coeff_out, dgen_to), ring_coeff in dstr.delta(dgen).items():
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

    def tensorDD(self, ddstr):
        """Compute the box tensor product DA * DD of this bimodule with the
        given type DD structure. Returns the resulting type DD structure. Uses
        delta() and deltaPrefix() functions of this type DA structure.

        """
        ddstr_result = SimpleDDStructure(F2, self.algebra1, ddstr.algebra2)
        # Compute list of generators in the box tensor product
        for gen_left in self.getGenerators():
            for gen_right in ddstr.getGenerators():
                if gen_left.idem2 == gen_right.idem1:
                    ddstr_result.addGenerator(DATensorDDGenerator(
                        ddstr_result, gen_left, gen_right))

        def search(start_gen, cur_ddgen, cur_algd, cur_coeffs_a):
            """Searching for an arrow in the box tensor product.
            - start_gen: starting generator in the box tensor product. The
              resulting arrow will start from here.
            - cur_ddgen: current location in the type DD structure.
            - cur_algd: current product algebra outputs on the right side of the
              DD structure.
            - cur_coeffs_a: current list of A-side inputs to the type DA
              structure (or alternatively, list of algebra outputs on the left
              side of the DD structure).

            """
            start_dagen, start_dgen = start_gen
            cur_delta = self.delta(start_dagen, cur_coeffs_a)
            for (coeff_d, gen_to), ring_coeff in cur_delta.items():
                ddstr_result.addDelta(start_gen, DATensorDDGenerator(
                    ddstr_result, gen_to, cur_ddgen), coeff_d, cur_algd, 1)
            if self.deltaPrefix(start_dagen, cur_coeffs_a):
                for (coeff_out1, coeff_out2, dgen_to), ring_coeff in \
                    ddstr.delta(cur_ddgen).items():
                    new_algd = cur_algd * coeff_out2
                    if new_algd != E0:
                        search(start_gen, dgen_to, new_algd.getElt(),
                               cur_coeffs_a + (coeff_out1,))

        for x in ddstr_result.getGenerators():
            dagen, ddgen = x
            search(x, ddgen, ddgen.idem2.toAlgElt(ddstr.algebra2), ())
            # Add arrows coming from idempotent output on the left DD-side
            for (coeff_out1, coeff_out2, dgen_to), ring_coeff in \
                ddstr.delta(ddgen).items():
                if coeff_out1.isIdempotent():
                    ddstr_result.addDelta(
                        x, DATensorDDGenerator(ddstr_result, dagen, dgen_to),
                        dagen.idem1.toAlgElt(self.algebra1), coeff_out2, 1)

        # Grading is omitted.
        return ddstr_result

    def registerHDiagram(self, diagram, base_gen, base_gr = None):
        """Associate the given diagram as the Heegaard diagram from which this
        type DA structure can be derived. We will attempt to match generators
        of the type DA structure to generators of the Heegaard diagram.
        Currently this is possible only if no two generators have the same
        idempotents (so the match can be made by comparing idempotents).

        As a result, computes grading of each generator from the Heegaard
        diagram and checks it against type DA operations. Attributes added are:
        *. self.hdiagram - the Heegaard diagram.
        *. self.hdiagram_gen_map - dictionary mapping generators in the type DA
        structure to generators in Heegaard diagram.
        *. self.gr_set - the grading set (of type SimpleDbGradingSet).
        *. self.grading - dictionary mapping generators in the type DA
        structure to their gradings.

        Requires self.getGenerators() to be implemented.

        """
        self.hdiagram = diagram
        # Get PMC's and check that they make sense
        hd_pmc1, hd_pmc2 = self.hdiagram.pmc_list
        das_pmc1, das_pmc2 = self.algebra1.pmc, self.algebra2.pmc
        assert hd_pmc1.opp() == das_pmc1
        assert hd_pmc2 == das_pmc2
        # Now attempt to match generators
        self.hdiagram_gen_map = dict()
        idem_size = 2 * das_pmc1.genus - len(base_gen.idem1)
        gens = self.getGenerators()
        hgens = diagram.getHFGenerators(idem_size = idem_size)
        for gen in gens:
            for hgen in hgens:
                hgen_idem1, hgen_idem2 = hgen.getDIdem()[0], hgen.getIdem()[1]
                if (gen.idem1, gen.idem2) == (hgen_idem1, hgen_idem2):
                    self.hdiagram_gen_map[gen] = hgen
                    break
            assert gen in self.hdiagram_gen_map
        # Compute grading and check consistency with algebra actions
        base_hgen = self.hdiagram_gen_map[base_gen]
        self.gr_set, gr = self.hdiagram.computeDAGrading(base_hgen, base_gr)
        self.grading = dict()
        for gen in gens:
            self.grading[gen] = gr[self.hdiagram_gen_map[gen]]
                
class SimpleDAStructure(DAStructure):
    """Represents a type DA structure with a finite number of generators and a
    finite number of type DA operations.

    """
    def __init__(self, ring, algebra1, algebra2,
                 side1 = ACTION_LEFT, side2 = ACTION_RIGHT):
        """Specifies the algebras and sides of the type DA action. algebra1 and
        side1 are for the type D action. algebra2 and side2 are for the type A
        action.

        """
        assert side1 == ACTION_LEFT and side2 == ACTION_RIGHT, \
            "Actions other than left/right are not implemented for DA."
        DAStructure.__init__(self, ring, algebra1, algebra2, side1, side2)
        self.generators = set()
        self.da_action = dict()

    def __len__(self):
        return len(self.generators)

    def getGenerators(self):
        return list(self.generators)

    def addGenerator(self, generator):
        """Add a generator. No effect if the generator already exists."""
        assert generator.parent == self
        assert isinstance(generator, DAGenerator)
        self.generators.add(generator)

    def addDelta(self, gen_from, gen_to, coeff_d, coeffs_a, ring_coeff):
        """Add ring_coeff * (coeff_d * gen_to) to the delta of gen_from, with
        coeffs_a as A-side inputs. The arguments gen_form, gen_to, and coeff_d
        should be generators, and coeff_a should be a list of generators.

        """
        assert gen_from.parent == self and gen_to.parent == self
        assert DAStructure.idemMatchDA(gen_from, gen_to, coeff_d, coeffs_a)
        coeffs_a = tuple(coeffs_a)
        target_gen = TensorGenerator((coeff_d, gen_to), self.AtensorM)
        if (gen_from, coeffs_a) not in self.da_action:
            self.da_action[(gen_from, coeffs_a)] = E0
        self.da_action[(gen_from, coeffs_a)] += target_gen.elt(ring_coeff)
        # Clean out the zero maps
        if self.da_action[(gen_from, coeffs_a)] == 0:
            del self.da_action[(gen_from, coeffs_a)]

    def delta(self, MGen, algGens):
        if len(algGens) == 1 and algGens[0].isIdempotent() and \
           algGens[0].left_idem == MGen.idem2:
            return MGen.idem1.toAlgElt() * MGen
        elif (MGen, algGens) not in self.da_action:
            return E0
        else:
            return self.da_action[(MGen, algGens)]

    def deltaPrefix(self, MGen, algGens):
        # Should be called only after all addDelta has completed
        if not hasattr(self, "strict_prefix"):
            self.strict_prefix = set()
            for (gen_from, coeffs_a), target in self.da_action.items():
                for i in range(len(coeffs_a)):
                    self.strict_prefix.add((gen_from, tuple(coeffs_a[0:i])))

        return (MGen, algGens) in self.strict_prefix

    def toDDStructure(self):
        """Convert this to a type DD structure over algebra1 and cobar of
        algebra2.

        """
        cobar2 = CobarAlgebra(self.algebra2)
        ddstr = SimpleDDStructure(self.ring, self.algebra1, cobar2)
        dagen_to_ddgen_map = dict()
        for gen in self.generators:
            ddgen = SimpleDDGenerator(ddstr, gen.idem1, gen.idem2, gen.name)
            dagen_to_ddgen_map[gen] = ddgen
            ddstr.addGenerator(ddgen)
        for (gen_from, coeffs_a), target in self.da_action.items():
            for (coeff_d, gen_to), ring_coeff in target.items():
                idem = None
                if len(coeffs_a) == 0:
                    idem = gen_from.idem2
                    assert idem == gen_to.idem2
                cobar_gen = TensorStarGenerator(coeffs_a, cobar2, idem)
                ddstr.addDelta(dagen_to_ddgen_map[gen_from],
                               dagen_to_ddgen_map[gen_to],
                               coeff_d, cobar_gen, ring_coeff)
        return ddstr

    def testDelta(self):
        """Verify the type DA structure equations."""
        return self.toDDStructure().testDelta()

    def __str__(self):
        result = "Type DA Structure.\n"
        for (gen_from, coeffs_a), target in self.da_action.items():
            result += "m(%s; %s) = %s\n" % (gen_from, coeffs_a, target)
        return result

    def toStrWithMultA(self, mult_a):
        """Print all arrows with the given multiplicities on the D side."""
        result = "Type DA Structure.\n"
        for (gen_from, coeffs_a), target in self.da_action.items():
            total_mult = sumColumns([coeff.multiplicity for coeff in coeffs_a],
                                    len(mult_a))
            if mult_a == total_mult:
                result += "m(%s; %s) = %s\n" % (gen_from, coeffs_a, target)
        return result

    def restrictToMultA(self, start, end):
        """Restrict actions to those with multiplicity in the interval
        (start, end). Here start and end specify points on the PMC. For example,
        start, end = 0, pmc.n-1 will not change the action. Returns the new type
        DA structure without changing the original structure.

        """
        translate_dict = dict()
        dastr = SimpleDAStructure(F2, self.algebra1, self.algebra2,
                                  self.side1, self.side2)
        for gen in self.generators:
            translate_dict[gen] = SimpleDAGenerator(
                dastr, gen.idem1, gen.idem2, gen.name)
            dastr.addGenerator(translate_dict[gen])

        mult_len = self.algebra1.pmc.n - 1
        for (gen_from, coeffs_a), target in self.da_action.items():
            total_mult = sumColumns([coeff.multiplicity for coeff in coeffs_a],
                                    mult_len)
            if all([total_mult[i] <= 0
                    for i in range(0, start) + range(end, mult_len)]):
                for (coeff_d, gen_to), ring_coeff in target.items():
                    dastr.addDelta(translate_dict[gen_from],
                                   translate_dict[gen_to],
                                   coeff_d, coeffs_a, ring_coeff)
        return dastr

    def checkGrading(self):
        for (x, coeffs_a), target in self.da_action.items():
            for (coeff_d, y), ring_coeff in target.items():
                gr_x1, gr_x2 = self.grading[x].data
                gr_y1, gr_y2 = self.grading[y].data
                for coeff_a in coeffs_a:
                    gr_x2 = gr_x2 * coeff_a.getGrading()
                gr_y1 = coeff_d.getGrading() * gr_y1
                new_gr_x = SimpleDbGradingSetElement(
                    self.gr_set, [gr_x1, gr_x2])
                new_gr_y = SimpleDbGradingSetElement(
                    self.gr_set, [gr_y1, gr_y2])
                assert new_gr_x - (1-len(coeffs_a)) == new_gr_y

def identityDA(pmc):
    """Returns the identity type DA structure for a given PMC."""
    alg = pmc.getAlgebra()
    dastr = SimpleDAStructure(F2, alg, alg)
    idems = pmc.getIdempotents()
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
    # Now add grading. Any generator can serve as base_gen
    for gen in dastr.getGenerators():
        base_gen = gen
        break
    dastr.registerHDiagram(getIdentityDiagram(pmc), base_gen)
    return dastr
