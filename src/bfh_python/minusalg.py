"""Strand algebra for the minus theory."""

from .algebra import Generator, SimpleChainComplex
from .algebra import E0
from .pmc import PMC, Strands, StrandAlgebra, StrandDiagram
from .pmc import splitPMC
from .utility import memorize
from .utility import F2

class MinusStrands(Strands):
    """The corresponding Strands class for the strand algebra in the minus
    theory. The main difference with Strands in pmc.py is that the moving
    strands can go through the basepoint, hence the starting point may be
    greater than the ending point (or equal, in the case where the moving strand
    covers the entire PMC).

    """
    def __init__(self, pmc, data):
        # Currently only supporting torus case
        assert pmc == splitPMC(1)

        self.pmc = pmc
        # Compute multiplicity.
        self.multiplicity = [0] * self.pmc.n
        for st in self:
            assert len(st) == 2
            start, end = st[0], st[1]
            if end <= start:
                end += self.pmc.n
            for pos in range(start, end):
                self.multiplicity[pos % self.pmc.n] += 1

class MinusStrandDiagram(StrandDiagram):
    """The corresponding Strands class for the strand algebra in the minus
    theory.

    """
    def __init__(self, parent, left_idem, strands, right_idem = None):
        """Be sure to use MinusStrands to convert strands."""
        if not isinstance(strands, MinusStrands):
            strands = MinusStrands(parent.pmc, strands)

        StrandDiagram.__init__(self, parent, left_idem, strands, right_idem)

class MinusStrandAlgebra(StrandAlgebra):
    """The corresponding Strands class for the strand algebra in the minus
    theory.

    """
    def __init__(self, ring, pmc):
        """Specifies the PMC. Assume multiplicity one and middle idempotent
        size.

        """
        StrandAlgebra.__init__(self, ring, pmc, pmc.genus, mult_one = True)

    @memorize
    def getGenerators(self):
        # Only implemented this case
        assert self.pmc == splitPMC(1)
        n = 4

        algebra = MinusStrandAlgebra(F2, self.pmc)
        result = []
        idems = self.pmc.getIdempotents(algebra.idem_size)
        for idem in idems:
            result.append(MinusStrandDiagram(algebra, idem, []))
        # Only one strand
        for start in range(n):
            for end in range(n):
                strands = MinusStrands(self.pmc, [(start, end)])
                for l_idem in idems:
                    if strands.leftCompatible(l_idem):
                        result.append(
                            MinusStrandDiagram(algebra, l_idem, strands))
        return result        

    @memorize
    def diff(self, gen):
        # Differential is zero in the torus algebra.
        return E0

    @memorize
    def multiply(self, gen1, gen2):
        if not isinstance(gen1, MinusStrandDiagram):
            return NotImplemented
        if not isinstance(gen2, MinusStrandDiagram):
            return NotImplemented
        assert gen1.parent == self and gen2.parent == self, \
            "Algebra not compatible."
        if gen1.right_idem != gen2.left_idem:
            return E0

        # Enforce the multiplicity one condition
        total_mult = [m1+m2 for m1, m2 in zip(gen1.multiplicity,
                                              gen2.multiplicity)]
        if not all([x <= 1 for x in total_mult]):
            return E0

        pmc = gen1.pmc
        new_strands = []

        # Keep track of which strands at right are not yet used.
        strands_right = list(gen2.strands)
        for sd in gen1.strands:
            mid_idem = pmc.pairid[sd[1]]
            possible_match = [sd2 for sd2 in strands_right
                              if pmc.pairid[sd2[0]] == mid_idem]
            if len(possible_match) == 0:
                new_strands.append(sd)
            else: # len(possible_match) == 1
                sd2 = possible_match[0]
                if sd2[0] != sd[1]:
                    return E0
                else:
                    new_strands.append((sd[0], sd2[1]))
                    strands_right.remove(sd2)

        new_strands.extend(strands_right)
        mult_term = MinusStrandDiagram(self, gen1.left_idem, new_strands,
                                       gen2.right_idem)
        # No problem with double crossing in the multiplicity one case.
        return mult_term.elt()

def minusSD(pmc, data):
    """Simple way to obtain a minus strand diagram. Each element of data is
    either an integer or a pair. An integer specifies a double horizontal at
    this position (and its paired position). A pair (p, q) specifies a strand
    from p to q.

    """
    parent = MinusStrandAlgebra(F2, pmc)
    assert parent.idem_size == len(data)
    left_idem = []
    strands = []
    for d in data:
        if isinstance(d, int):
            left_idem.append(pmc.pairid[d])
        else:
            left_idem.append(pmc.pairid[d[0]])
            strands.append(d)
    return MinusStrandDiagram(parent, left_idem, strands)


def getHalfIdComplex():
    """Returns the chain complex underlying the proposed quasi-inverse of the
    type DD bimodule for half-identity.

    """
    class LargeComplexGenerator(Generator, tuple):
        """A generator of the large chain complex. Specified by two strand
        diagrams in the given PMC.

        """
        def __new__(cls, parent, sd_left, sd_right):
            return tuple.__new__(cls, (sd_left, sd_right))

        def __init__(self, parent, sd_left, sd_right):
            "Specifies the two strand diagrams."""
            # Note tuple initialization is automatic
            Generator.__init__(self, parent)

    # Dictionary mapping total multiplicity profile to chain complexes.
    partial_cxs = dict()
    pmc = splitPMC(1)
    alg = MinusStrandAlgebra(F2, pmc)

    # Find the set of generators.
    alg_gens = alg.getGenerators()
    for gen_left in alg_gens:
        for gen_right in alg_gens:
            if gen_left.getLeftIdem() == gen_right.getLeftIdem().comp():
                total_mult = [a + b for a, b in zip(
                    gen_left.multiplicity, gen_right.multiplicity)]
                total_mult = tuple(total_mult)
                if total_mult not in partial_cxs:
                    partial_cxs[total_mult] = SimpleChainComplex(F2)
                cur_gen = LargeComplexGenerator(
                    partial_cxs[total_mult], gen_left, gen_right)
                partial_cxs[total_mult].addGenerator(cur_gen)

    def hasDifferential(gen_from, gen_to):
        """Determine whether there is a differential between two generators
        of the large chain complex.

        """
        left_from, right_from = gen_from
        left_to, right_to = gen_to
        mult_left_from, mult_left_to = \
            left_from.multiplicity, left_to.multiplicity
        diff = [a-b for a, b in zip(mult_left_from, mult_left_to)]
        if all([n == 0 for n in diff]):
            if left_from == left_to and right_to in right_from.diff():
                return True
            if right_from == right_to and left_from in left_to.diff():
                return True
        elif all([n == 0 or n == 1 for n in diff]):
            pos_one = [i for i in range(len(diff)) if diff[i] == 1]
            if len(pos_one) != 1:
                return False
            start, end = pos_one[0], (pos_one[0]+1)%4
            st_move = MinusStrands(pmc, [(start, end)])
            if not st_move.rightCompatible(left_to.getLeftIdem()):
                return False
            left_move = MinusStrandDiagram(
                alg, None, st_move, left_to.getLeftIdem())
            if not st_move.rightCompatible(right_from.getLeftIdem()):
                return False
            right_move = MinusStrandDiagram(
                alg, None, st_move, right_from.getLeftIdem())
            return left_move * left_to == 1*left_from and \
                right_move * right_from == 1*right_to
        else:
            return False

    # Compute differentials
    for total_mult, cx in list(partial_cxs.items()):
        gens = cx.getGenerators()
        for gen_from in gens:
            for gen_to in gens:
                if hasDifferential(gen_from, gen_to):
                    cx.addDifferential(gen_from, gen_to, 1)
        print(total_mult, cx)
        cx.checkDifferential()
        cx.simplify()
