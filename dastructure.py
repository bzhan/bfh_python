"""Defines type DA structures."""

from hdiagram import *

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
    an element in the tensor module Tensor((A,M)).

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
    finithe number of type DA operations.

    """
    pass
