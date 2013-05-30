"""Defines type DA structures."""

from ddstructure import *

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

    def rmultiply(self, MGen, AGen):
        """Multiply a generator of type DAStructure with an algebra generator
        means forming the tensor.

        """
        return 1*TensorGenerator((Agen, MGen), self.AtensorM)

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

    def testDelta(self):
        """Verify the type DA structure equations."""
        # Convert this to a type DD structure over algebra1 and cobar of
        # algebra2
        cobar2 = CobarAlgebra(self.algebra2)
        ddstr = SimpleDDStructure(self.ring, self.algebra1, cobar2)
        dagen_to_ddgen_map = dict()
        for gen in self.generators:
            ddgen = SimpleDDGenerator(ddstr, gen.idem1, gen.idem2, gen.name)
            dagen_to_ddgen_map[gen] = ddgen
            ddstr.addGenerator(ddgen)
        for (gen_from, coeffs_a), target in self.da_action.items():
            for (coeffs_d, gen_to), ring_coeff in target.items():
                idem = None
                if len(coeffs_a) == 0:
                    idem = gen_from.idem2
                    assert idem == gen_to.idem2
                cobar_gen = TensorStarGenerator(coeffs_a, cobar2, idem)
                ddstr.addDelta(dagen_to_ddgen_map[gen_from],
                               dagen_to_ddgen_map[gen_to],
                               coeffs_d, cobar_gen, ring_coeff)
        return ddstr.testDelta()

    def __str__(self):
        result = "Type DA Structure.\n"
        for (gen_from, coeffs_a), target in self.da_action.items():
            result += "m(%s; %s) = %s\n" % (gen_from, coeffs_a, target)
        return result

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
    return dastr
