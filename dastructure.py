"""Defines type DA structures."""

from algebra import CobarAlgebra, DGAlgebra, FreeModule, Generator, Tensor, \
    TensorGenerator, TensorStarGenerator
from algebra import E0
from dstructure import DGenerator, SimpleDStructure
from ddstructure import SimpleDDGenerator, SimpleDDStructure
from grading import GeneralGradingSet, GeneralGradingSetElement, \
    SimpleDbGradingSetElement
from hdiagram import getIdentityDiagram
from pmc import StrandDiagram
from utility import NamedObject
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
        return "%s,%s" % (str(self.idem1), str(self.idem2))

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
        """Specify generators on two sides of the tensor (DD and D generators).

        """
        # Note tuple initialization is automatic
        DGenerator.__init__(self, parent, gen_left.idem1)

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
        return NotImplementedError("toSimpleDAStructure not yet implemented.")

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
        assert coeff_d.getRightIdem() == gen_to.idem1
        assert coeff_d.getLeftIdem() == gen_from.idem1
        if len(coeffs_a) == 0:
            assert gen_from.idem2 == gen_to.idem2
        else:
            assert gen_from.idem2 == coeffs_a[0].getLeftIdem()
            for i in range(len(coeffs_a)-1):
                assert coeffs_a[i].getRightIdem() == \
                    coeffs_a[i+1].getLeftIdem()
            assert coeffs_a[-1].getRightIdem() == gen_to.idem2
        coeffs_a = tuple(coeffs_a)
        target_gen = TensorGenerator((coeff_d, gen_to), self.AtensorM)
        if (gen_from, coeffs_a) not in self.da_action:
            self.da_action[(gen_from, coeffs_a)] = E0
        self.da_action[(gen_from, coeffs_a)] += target_gen.elt(ring_coeff)
        # Clean out the zero maps
        for source, target in self.da_action.items():
            if target == 0:
                del self.da_action[source]

    def delta(self, MGen, algGens):
        if len(algGens) == 1 and algGens[0].isIdempotent() and \
           algGens[0].left_idem == MGen.idem2:
            return MGen.idem1.toAlgElt() * MGen
        elif (MGen, algGens) not in self.da_action:
            return E0
        else:
            return self.da_action[(MGen, algGens)]

    def deltaPrefix(self, MGen, algGens):
        # Very inefficient implementation.
        for (gen_from, coeffs_a), target in self.da_action.items():
            if MGen == gen_from and len(algGens) < len(coeffs_a) and \
               algGens == coeffs_a[:len(algGens)]:
                return True
        return False

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

def AddChordToDA(dastr, coeff_d, coeffs_a):
    alg1 = dastr.algebra1
    alg2 = dastr.algebra2
    gen_set = dastr.getGenerators()
    for x in gen_set:
        for y in gen_set:
            if alg1.mult_one and not coeff_d.isMultOne():
                continue
            if alg2.mult_one and \
               not all([coeff_a.isMultOne() for coeff_a in coeffs_a]):
                continue
            if not coeff_d.idemCompatible(x.idem1, y.idem1):
                continue
            alg_d = StrandDiagram(alg1, x.idem1, coeff_d)
            # Now test that the sequence of A inputs match the idempotents.
            # Construct the algebra inputs at the same time.
            alg_a = []
            cur_idem = x.idem2
            idem_ok = True
            for coeff_a in coeffs_a:
                next_idem = coeff_a.propagateRight(cur_idem)
                if next_idem is None:
                    idem_ok = False
                    break
                alg_a.append(StrandDiagram(alg2, cur_idem, coeff_a))
                cur_idem = next_idem
            if idem_ok and cur_idem == y.idem2:
                dastr.addDelta(x, y, alg_d, alg_a, 1)

def DAStrFromChords(alg1, alg2, idem_pairs, chord_pairs):
    """Construct type DA structure from list of idempotent pairs and chord
    pairs. This is usually used to construct the 'initial' type DA structure,
    to be completed by DAStrFromBasis.
    - idem_pairs is list of pairs of Idempotent. Left idempotent is D side,
    right idempotent is A side.
    - chord_pairs is list of pairs (coeff_d, coeffs_a), where coeff_d is of
    class Strands (for D side output), and coeffs_a is a list of Strands (for A
    side input).

    """
    dastr = SimpleDAStructure(F2, alg1, alg2)
    for i in range(len(idem_pairs)):
        dastr.addGenerator(SimpleDAGenerator(
            dastr, idem_pairs[i][0], idem_pairs[i][1], i))
    for coeff_d, coeffs_a in chord_pairs:
        AddChordToDA(dastr, coeff_d, coeffs_a)
    return dastr
