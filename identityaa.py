"""Description of type AA structure for identity cobordism."""

from pmc import *

# Convenient values for specifying sides. Use only in the context of large
# complex generators.
_LEFT, _RIGHT = 0, 1

class HomotopyAA:
    """Contains the description of the large chain complex for type AA
    identity, and the homotopy needed to simplify it.

    """
    def __init__(self, pmc):
        """Specifies the PMC."""
        self.pmc = pmc
        self.pmc_alg = pmc.getAlgebra()

    class LargeComplexGenerator(Generator, tuple):
        """A generator of the large chain complex for type AA identity.
        Specified by two strand diagrams in the given PMC.

        """
        def __new__(cls, parent, sd_left, sd_right):
            return tuple.__new__(cls, (sd_left, sd_right))

        def __init__(self, parent, sd_left, sd_right):
            "Specifies the two strand diagrams."""
            # Note tuple initialization is automatic
            Generator.__init__(self, parent)

        def replaceLeft(self, new_left):
            """Get LargeComplexGenerator with same parent and right strand
            diagram, but new left strand diagram.

            """
            return self.__class__(self.parent, new_left, self[1])

        def replaceRight(self, new_right):
            """Get LargeComplexGenerator with same parent and left strand
            diagram, but new right strand diagram.

            """
            return self.__class__(self.parent, self[0], new_right)

    def _computeLargeChainComplex(self):
        """Computes the large chain complex for type AA identity."""
        # Dictionary mapping total multiplicity profile to chain complexes.
        self.partial_cxs = dict()

        # Find the set of generators.
        alg_gens = self.pmc_alg.getGenerators()
        for gen_left in alg_gens:
            for gen_right in alg_gens:
                if gen_left.getLeftIdem() == gen_right.getLeftIdem().comp():
                    total_mult = [a + b for a, b in zip(
                            gen_left.multiplicity, gen_right.multiplicity)]
                    total_mult = tuple(total_mult)
                    if total_mult not in self.partial_cxs:
                        self.partial_cxs[total_mult] = SimpleChainComplex(F2)
                    cur_gen = self.LargeComplexGenerator(
                        self.partial_cxs[total_mult], gen_left, gen_right)
                    self.partial_cxs[total_mult].addGenerator(cur_gen)

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
                start, end = pos_one[0], pos_one[-1]+1
                if pos_one == range(start, end):
                    st_move = Strands(self.pmc, [(start, end)])
                    if not st_move.rightCompatible(left_to.getLeftIdem()):
                        return False
                    left_move = StrandDiagram(
                        self.pmc_alg, None, st_move, left_to.getLeftIdem())
                    if not st_move.rightCompatible(right_from.getLeftIdem()):
                        return False
                    right_move = StrandDiagram(
                        self.pmc_alg, None, st_move, right_from.getLeftIdem())
                    return left_move * left_to == 1*left_from and \
                        right_move * right_from == 1*right_to
            else:
                return False

        # Compute differentials
        for total_mult, cx in self.partial_cxs.items():
            gens = cx.getGenerators()
            for gen_from in gens:
                for gen_to in gens:
                    if hasDifferential(gen_from, gen_to):
                        cx.addDifferential(gen_from, gen_to, 1)
            cx.checkDifferential()

    def getHomotopyFrom(self, gen_from):
        """Returns the list of generators gen_from has homotopy to."""
        sd_left, sd_right = gen_from
        return [self.LargeComplexGenerator(gen_from.parent, left, right)
                for left, right in homotopyMap(sd_left, sd_right)]

    def _computeHomotopy(self):
        """Computes the homotopy map on the large complex."""
        self.partial_hts = dict()

        for total_mult, cx in self.partial_cxs.items():
            if all([n == 0 for n in total_mult]):
                continue # No homotopy needed
            cur_homotopy = SimpleChainMorphism(cx, cx)
            gens = cx.getGenerators()

            # Now find the homotopy map
            for gen_from in gens:
                homotopy_to = self.getHomotopyFrom(gen_from)
                for gen_to in homotopy_to:
                    cur_homotopy.addMorphism(gen_from, gen_to, 1)
            self.partial_hts[total_mult] = cur_homotopy

    def testHomotopy(self):
        """Test the identity for homotopy."""
        self._computeLargeChainComplex()
        self._computeHomotopy()
        for total_mult, cx in self.partial_cxs.items():
            if all([n == 0 for n in total_mult]):
                continue # No homotopy needed
            cur_homotopy = self.partial_hts[total_mult]
            for gen in cx.getGenerators():
                assert cur_homotopy.apply(gen.diff()) + \
                    cur_homotopy.apply(gen).diff() == 1 * gen

def _getIntervalOrdering(pmc):
    """Get the special ordering of intervals needed for computing the
    homotopy (the order on traversing the glued circle).

    """
    cur_pt = pmc.n - 1
    result = []
    while True:
        cur_pt = pmc.otherp[cur_pt]
        if cur_pt == 0:
            assert len(result) == pmc.n-1
            break
        result.append(cur_pt - 1)
        cur_pt -= 1
    return result

def _moveLeft(sd_pair, a, b):
    """Move a strand (a, b) from right side to the left side. Must be
    successful.

    """
    sd_left, sd_right = sd_pair
    new_right = _factor(sd_right, a, b)
    assert new_right is not None
    to_mult = StrandDiagram(sd_left.parent, None, [(a,b)],
                            sd_left.getLeftIdem())
    new_left = to_mult * sd_left
    assert new_left != E0
    new_left = new_left.getElt()
    return (new_left, new_right)

def _moveRight(sd_pair, a, b):
    """Move a strand (a, b) from left side to the right side. Must be
    successful.

    """
    sd_left, sd_right = sd_pair
    new_left = _factor(sd_left, a, b)
    assert new_left is not None
    to_mult = StrandDiagram(sd_left.parent, None, [(a,b)],
                            sd_right.getLeftIdem())
    new_right = to_mult * sd_right
    assert new_right != E0
    new_right = new_right.getElt()
    return (new_left, new_right)

def _uncross(sd, start1, end1, start2, end2):
    """Uncross the two strands (start1, end1) and (start2, end2) of the given
    strand diagram.

    """
    strand_lst = list(sd.strands)
    if start1 != end1:
        strand_lst.remove((start1, end1))
    if start2 != end2:
        strand_lst.remove((start2, end2))
    strand_lst.append((start1, end2))
    strand_lst.append((start2, end1))
    result_sd = StrandDiagram(sd.parent, sd.getLeftIdem(), strand_lst)
    assert result_sd in sd.diff()
    return result_sd

def _cross(sd, start1, end1, start2, end2):
    """Cross the two strands (start1, end1) and (start2, end2) of the
    strand diagram on the given side of start_gen.

    """
    strand_lst = list(sd.strands)
    strand_lst.remove((start1, end1))
    strand_lst.remove((start2, end2))
    if start1 != end2:
        strand_lst.append((start1, end2))
    if start2 != end1:
        strand_lst.append((start2, end1))
    result_sd = StrandDiagram(sd.parent, sd.getLeftIdem(), strand_lst)
    assert sd in result_sd.diff()
    return result_sd

def _factor(sd, a, b):
    """Try to factor off a strand (a, b) (from the left) from sd. Returns
    the new strand diagram if successful. Otherwise returns None.

    """
    strand_lst = list(sd.strands)
    for p, q in [(p,q) for p,q in strand_lst if p == a and q >= b]:
        strand_lst.remove((p, q))
        if q > b:
            strand_lst.append((b, q))
        result = StrandDiagram(sd.parent, None, strand_lst,
                               sd.getRightIdem())
        to_mult = StrandDiagram(sd.parent, None, [(a,b)],
                                result.getLeftIdem())
        if to_mult * result == 1 * sd:
            return result
        else:
            return None

def homotopyMap(sd_left, sd_right):
    """Implements the homotopy map. Input is a pair of strand diagrams from the
    same PMC. Outputs a list of (from zero to two) pairs that the homotopy maps
    the input pair to.

    """
    # Two helper functions
    def startType(sd, pair):
        """Returns -1 if sd has horizontal strands at pair, ``n`` if sd has
        strand starting at point ``n`` in the pair, -2 otherwise (left
        idempotent of sd is not occupied at pair).

        """
        for p, q in sd.strands:
            if sd.pmc.pairid[p] == pair:
                return p
        if pair in sd.left_idem:
            return -1
        return -2

    def getStrandsAtPoint(sd, point):
        """List of strands (recorded as pair (p,q)) covering the point
        key_pos+1, including any double horizontal at point key_pos+1
        (recorded as (key_pos+1, key_pos+1)).

        """
        result = [(p, q) for p, q in sd.strands if p <= point <= q]
        if sd.pmc.pairid[point] in sd.double_hor:
            result.append((point, point))
        return result

    # To start, we find key_pos and other important locations
    result = []
    pmc = sd_left.pmc
    assert pmc == sd_right.pmc
    total_mult = [a + b for a, b in zip(sd_left.multiplicity,
                                        sd_right.multiplicity)]
    ordering = _getIntervalOrdering(pmc)

    # Index of first interval with multiplicity >= 2 (note multiplicity
    # cannot jump by more than 1.
    lowest_two = find(total_mult, 2)
    # Compute key_pos for two cases
    if lowest_two == -1:
        # Total multiplicity one case
        unoccupied_id = [i for i in range(len(ordering))
                         if total_mult[ordering[i]] != 0]
        assert len(unoccupied_id) > 0
        key_pos = ordering[unoccupied_id[0]]
        # Always use the lowest / highest possible interval, does not appear
        # helpful.
        # key_pos = -1
        # for i in range(len(ordering)):
        #     if total_mult[ordering[i]] != 0 and \
        #             (i == 0 or total_mult[ordering[i-1]] == 0):
        #         if key_pos == -1 or ordering[i] > key_pos:
        #             key_pos = ordering[i]
        total_mult_one = True
    else:
        # Total multiplicity >1 case
        key_pos = lowest_two - 1
        total_mult_one = False
    # Compute key_pair and forbid_pos
    key_pair = pmc.pairid[key_pos+1]
    forbid_pos = pmc.otherp[key_pos+1]
    assert total_mult[key_pos] != 0
    if total_mult_one:
        assert forbid_pos == pmc.n-1 or total_mult[forbid_pos] == 0

    # Now, some more specific information about this generator
    type_left, type_right = \
        startType(sd_left, key_pair), startType(sd_right, key_pair)
    assert type_left != forbid_pos and type_right != forbid_pos

    left_avail = getStrandsAtPoint(sd_left, key_pos+1)
    right_avail = getStrandsAtPoint(sd_right, key_pos+1)
    total_avail = sorted([((p,q), _LEFT) for p,q in left_avail] +
                         [((p,q), _RIGHT) for p,q in right_avail])
    assert len(total_avail) >= 2
    (start1, end1), side1 = total_avail[0]
    (start2, end2), side2 = total_avail[1]

    # Homotopies for the three usual cases
    if side1 == _LEFT and side2 == _LEFT:
        # First case: both lowest available strands are on the left.
        # Un-cross two strands on the left.
        if end1 <= end2:
            return result
        result.append((_uncross(sd_left, start1, end1, start2, end2),sd_right))
    elif side1 == _RIGHT and side2 == _RIGHT:
        # Second case: both lowest available strands are on the right.
        # Cross two strands on the right.
        if end1 >= end2:
            return result
        result.append((sd_left, _cross(sd_right, start1, end1, start2, end2)))
    elif side1 == _RIGHT and side2 == _LEFT:
        # Third (and fourth case in paper). Lowest strand on the right
        # and the second lowest on the left. Move strand to the left.
        result.append(_moveLeft((sd_left, sd_right), start1, start2))
    else:
        return result

    if not total_mult_one:
        return result
    # Three special cases only for the multiplicity one case
    forbid_start = [a for a,b in sd_left.strands if b == forbid_pos]
    if not forbid_start:
        return result
    forbid_start = forbid_start[0]
    gen_start = (sd_left, sd_right)
    if side1 == _RIGHT and side2 == _RIGHT:
        gen_step1 = (sd_left, _cross(sd_right, start1, end1, start2, end2))
        gen_step2 = _moveRight(gen_step1, forbid_start, forbid_pos)
        gen_step3 = _moveLeft(gen_step2, start1, end1)
    elif side1 == _RIGHT and side2 == _LEFT and forbid_start != key_pos+1:
        gen_step1 = _moveLeft(gen_start, start1, key_pos+1)
        gen_step2 = _moveRight(gen_step1, forbid_start, forbid_pos)
        gen_step3 = (_uncross(gen_step2[0], start1, end2, key_pos+1, key_pos+1)
                     ,gen_step2[1])
    elif side1 == _RIGHT and side2 == _LEFT and forbid_start == key_pos+1:
        gen_step1 = _moveLeft(gen_start, start1, key_pos+1)
        gen_step2 = _moveRight(gen_step1, start1, forbid_pos)
        gen_step3 = _moveLeft(gen_step2, start1, key_pos+1)
    else:
        return result
    result.append(gen_step3)
    return result
