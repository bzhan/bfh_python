"""Handles grading groups and grading sets."""

from fractions import Fraction, gcd
from numbers import Number
from linalg import RowSystem
from utility import flatten, grTypeStr, memorize, oppSide, sideStr, tolist
from utility import ACTION_LEFT, ACTION_RIGHT, BIG_GRADING, SMALL_GRADING

class Group:
    """Represents a general group."""
    def multiply(self, elt1, elt2):
        """Returns the product of gen1 and gen2."""
        raise NotImplementedError("Multiply not implemented.")

class GroupElement:
    """Represents an element of a group."""
    def __init__(self, parent):
        """Specifies which group this element is in."""
        self.parent = parent

    def __mul__(self, other):
        """Multiplies this group element with other."""
        return self.parent.multiply(self, other)

class BigGradingGroup(Group):
    """Big grading group associated to a PMC."""
    def __init__(self, pmc):
        self.pmc = pmc
        self.type = BIG_GRADING
        self.spinc_len = self.pmc.n - 1

    def __eq__(self, other):
        return self.pmc == other.pmc

    def __ne__(self, other):
        return not (self == other)

    def __hash__(self):
        return hash((self.pmc, "BigGradingGroup"))

    def multiply(self, elt1, elt2):
        if not isinstance(elt1, BigGradingElement):
            return NotImplemented
        if not isinstance(elt2, BigGradingElement):
            return NotImplemented
        assert elt1.parent == self and elt2.parent == self
        m1 = elt1.spinc
        m2 = elt2.spinc
        new_maslov = 0
        for i in range(len(m1)-1):
            new_maslov += (m1[i] * m2[i+1] - m1[i+1] * m2[i])
        new_maslov /= Fraction(2)
        new_maslov += (elt1.maslov + elt2.maslov)
        new_spinc = [a+b for a, b in zip(m1, m2)]
        return BigGradingElement(self, new_maslov, new_spinc)

    def opp(self):
        """Returns the big grading group associated to the opposite PMC."""
        return BigGradingGroup(self.pmc.opp())

    def zero(self):
        """Returns the zero element of this grading group."""
        return BigGradingElement(self, 0, [0]*self.spinc_len)

    def central(self):
        """Returns the central element (lambda) of this grading group. This has
        maslov component 1 and spinc component zero.

        """
        return BigGradingElement(self, 1, [0]*(self.pmc.n-1))

    def basis(self, i):
        """Returns i'th basis element of the spinc component, with maslov
        component zero.

        """
        spinc_vec = [0] * self.spinc_len
        spinc_vec[i] = 1
        return BigGradingElement(self, 0, spinc_vec)

class BigGradingElement(GroupElement):
    """An element of the big grading group."""
    def __init__(self, parent, maslov, spinc):
        """Specifies the maslov and spinc component of the grading. The spinc
        component is a list of pmc.n-1 multiplicities.

        """
        GroupElement.__init__(self, parent)
        self.maslov = maslov
        self.spinc = list(spinc)
        assert len(self.spinc) == self.parent.pmc.n - 1

    def __eq__(self, other):
        if isinstance(other, int) and other == 0:
            return self.maslov == 0 and all([n == 0 for n in self.spinc])
        return self.parent == other.parent and self.maslov == other.maslov \
            and self.spinc == other.spinc

    def __ne__(self, other):
        return not (self == other)

    def __hash__(self):
        return hash((self.parent, self.maslov, self.spinc))

    def __str__(self):
        return "[%s; %s]" % (str(self.maslov),
                             ", ".join(str(n) for n in self.spinc))

    def __repr__(self):
        return str(self)

    def opp(self):
        """Returns the corresponding grading element in the opposite PMC. Keeps
        the sign of both Maslov and Spin-c components.

        """
        return BigGradingElement(self.parent.opp(), self.maslov,
                                 reversed(self.spinc))

    def Ropp(self):
        """Returns the corresponding grading element in the opposite PMC. Keeps
        the Maslov component and reverses the sign of Spin-c components. (This
        is the R(.) operator in the papers).

        """
        return BigGradingElement(self.parent.opp(), self.maslov,
                                 reversed([-n for n in self.spinc]))

    def inverse(self):
        """Returns the inverse of this grading element. Reverses Maslov and
        Spin-c components.

        """
        return BigGradingElement(self.parent, -self.maslov,
                                 [-n for n in self.spinc])

    def power(self, exp):
        """Returns this grading element raised to the given power (basically
        multiplies the maslov and each spinc component by exp.

        """
        return BigGradingElement(self.parent, self.maslov * exp,
                                 [n * exp for n in self.spinc])

    def toSmallGrading(self):
        """Returns the corresponding small grading element. Should be called
        only for elements already in the small grading group (that is,
        M_*\delta(s) = 0).

        """
        pmc = self.parent.pmc
        gr_group = SmallGradingGroup(pmc)
        mult_tmp = list(self.spinc)
        small_mult = []
        for p, q in pmc.pairs:
            # The first entry of pairs should be in ascending order
            lead_mult = mult_tmp[p]
            small_mult.append(lead_mult)
            for pos in range(p, q):
                mult_tmp[pos] -= lead_mult
        assert all([n == 0 for n in mult_tmp])
        return SmallGradingElement(gr_group, self.maslov, small_mult)

class SmallGradingGroup(Group):
    """Small grading group associated to a PMC."""
    def __init__(self, pmc):
        self.pmc = pmc
        self.type = SMALL_GRADING
        self.spinc_len = self.pmc.num_pair

    def __str__(self):
        return "Small grading group over PMC %s" % str(self.pmc)

    def __repr__(self):
        return str(self)

    def __eq__(self, other):
        return self.pmc == other.pmc

    def __ne__(self, other):
        return not (self == other)

    def __hash__(self):
        return hash((self.pmc, "SmallGradingGroup"))

    def opp(self):
        """Returns the small grading group associated to the opposite PMC."""
        return SmallGradingGroup(self.pmc.opp())

    def multiply(self, elt1, elt2):
        if not isinstance(elt1, SmallGradingElement):
            return NotImplemented
        if not isinstance(elt2, SmallGradingElement):
            return NotImplemented
        assert elt1.parent == self and elt2.parent == self
        m1 = elt1.spinc
        m2 = elt2.spinc
        new_spinc = [a+b for a, b in zip(m1, m2)]
        new_maslov = elt1.maslov + elt2.maslov
        # Note this relies on the ordering of the pairs and the points inside
        # each pair.
        for i in range(self.pmc.num_pair):
            for j in range(i+1, self.pmc.num_pair):
                if self.pmc.pairs[i][1] > self.pmc.pairs[j][0] and \
                   self.pmc.pairs[i][1] < self.pmc.pairs[j][1]:
                    new_maslov += (m1[i] * m2[j] - m1[j] * m2[i])
        return SmallGradingElement(self, new_maslov, new_spinc)

    def zero(self):
        """Returns the zero element of this grading group."""
        return SmallGradingElement(self, 0, [0]*self.pmc.num_pair)

    def central(self):
        """Returns the central element (lambda) of this grading group. This has
        maslov component 1 and spinc component zero.

        """
        return SmallGradingElement(self, 1, [0]*self.pmc.num_pair)

    def basis(self, i):
        """Returns i'th basis element of the spinc component, with maslov
        component zero.

        """
        spinc_vec = [0] * self.spinc_len
        spinc_vec[i] = 1
        return SmallGradingElement(self, 0, spinc_vec)

class SmallGradingElement(GroupElement):
    """An element of the small grading group."""
    def __init__(self, parent, maslov, spinc):
        """Specifies the maslov and spinc component of the grading. The spinc
        component is a list of pmc.num_pair multiplicities.

        """
        GroupElement.__init__(self, parent)
        self.maslov = maslov
        self.spinc = list(spinc)
        assert len(self.spinc) == self.parent.pmc.num_pair

    def __eq__(self, other):
        if isinstance(other, int) and other == 0:
            return self.maslov == 0 and all([n == 0 for n in self.spinc])
        return self.parent == other.parent and self.maslov == other.maslov \
            and self.spinc == other.spinc

    def __ne__(self, other):
        return not (self == other)

    def __hash__(self):
        return hash((self.parent, self.maslov, self.spinc))

    def __str__(self):
        return "[%s; %s]" % (str(self.maslov),
                             ", ".join(str(n) for n in self.spinc))

    def __repr__(self):
        return str(self)

    def opp(self):
        """Returns the corresponding grading element in the opposite PMC. Keeps
        the sign of both Maslov and Spin-c components.

        """
        return self.toBigGrading().opp().toSmallGrading()

    def Ropp(self):
        """Returns the corresponding grading element in the opposite PMC. Keeps
        the Maslov component and reverses the sign of Spin-c components. (This
        is the R(.) operator in the papers).

        """
        return self.toBigGrading().Ropp().toSmallGrading()

    def inverse(self):
        """Returns the inverse of this grading element. Reverses Maslov and
        Spin-c components.

        """
        return SmallGradingElement(self.parent, -self.maslov,
                                   [-n for n in self.spinc])

    def power(self, exp):
        """Returns this grading element raised to the given power (basically
        multiplies the maslov and each spinc component by exp.

        """
        return SmallGradingElement(self.parent, self.maslov * exp,
                                   [n * exp for n in self.spinc])

    def toBigGrading(self):
        """Returns the corresponding big grading element."""
        pmc = self.parent.pmc
        gr_group = BigGradingGroup(pmc)
        big_mult = [0] * (pmc.n - 1)
        for i in range(pmc.num_pair):
            p, q = pmc.pairs[i]
            for pos in range(p, q):
                big_mult[pos] += self.spinc[i]
        return BigGradingElement(gr_group, self.maslov, big_mult)

class GradingRefinement(dict):
    """Represents a fixed grading refinement for a PMC, as a dictionary from
    idempotents to big grading elements.

    """
    def __init__(self, pmc, idem_size, data):
        """data is a dictionary mapping idempotents (in the given PMC) to big
        grading elements.

        """
        dict.__init__(self, data)
        self.pmc = pmc
        self.idem_size = idem_size
        self._checkRefinement()

    def _checkRefinement(self):
        """Check the validity of this refinement. It suffices to check that,
        for each pair of idempotents (s, t), we have
           M_*\delta(\psi(s) - \psi(t)) = s - t,
        where \psi is the refinement map, \delta takes the differential of the
        corresponding one chain on the border, and M_* is the pairing map. It
        is also enough to check this only between one idempotent and all
        others.

        """
        idems = self.pmc.getIdempotents(self.idem_size)
        ref_idem = idems[0]
        for idem in idems:
            gr_diff = self[idem] * self[ref_idem].inverse()
            # Take \delta map
            pt_mult = [-gr_diff.spinc[0]]
            for i in range(1, self.pmc.n-1):
                pt_mult.append(gr_diff.spinc[i-1] - gr_diff.spinc[i])
            pt_mult.append(gr_diff.spinc[-1])
            # Take M_* map
            idem_mult = [0] * self.pmc.num_pair
            for i in range(self.pmc.n):
                idem_mult[self.pmc.pairid[i]] += pt_mult[i]
            # Now compare with idem - ref_idem
            expected_idem_mult = [0] * self.pmc.num_pair
            for pair in idem:
                expected_idem_mult[pair] += 1
            for pair in ref_idem:
                expected_idem_mult[pair] -= 1
            assert idem_mult == expected_idem_mult

@memorize
def standardRefinement(pmc, idem_size = None):
    """Returns the "standard" refinement of the given PMC. For this refinement,
    add together the difference between the lower point for each pair and the
    lower point for each of the idem_size pairs in the PMC.

    """
    # This should not be used as the default refinement, as it does not
    # guarantee sd.small_gr.opp() == sd.opp().small_gr for a strand diagram sd.
    data = dict()
    gr_group = BigGradingGroup(pmc)
    if idem_size is None:
        idem_size = pmc.genus
    for idem in pmc.getIdempotents(idem_size):
        mult = [0] * (pmc.n - 1)
        maslov = 0
        for i in range(idem_size):
            p, q = pmc.pairs[idem[i]]
            p_s, q_s = pmc.pairs[i]
            assert p < q and p_s < q_s
            if p > p_s:
                for pos in range(p_s, p):
                    mult[pos] += 1
            elif p < p_s:
                for pos in range(p, p_s):
                    mult[pos] -= 1
            if idem[i] >= idem_size:
                maslov += Fraction(1, 2)
        data[idem] = BigGradingElement(gr_group, maslov, mult)
    return GradingRefinement(pmc, idem_size, data)

@memorize
def lowerRefinement(pmc, idem_size = None):
    """Returns the "lower" refinement of the given PMC. For this refinement,
    add together the interval from bottom to the lower point for each pair in
    the idempotent.

    """
    # This should not be used as the default refinement, as it does not
    # guarantee sd.small_gr.opp() == sd.opp().small_gr for a strand diagram sd.
    data = dict()
    gr_group = BigGradingGroup(pmc)
    for idem in pmc.getIdempotents(idem_size):
        mult = [0] * (pmc.n - 1)
        for pair in idem:
            p, q = pmc.pairs[pair]
            assert p < q
            for pos in range(0, p):
                mult[pos] += 1
        data[idem] = BigGradingElement(gr_group, 0, mult)
    return GradingRefinement(pmc, idem_size, data)

@memorize
def averageRefinement(pmc, idem_size = None):
    """Returns the "average" refinement of the given PMC. In this refinement,
    for each point in the idempotent, add 1/8 the interval from bottom to that
    point and subtract 1/8 the interval from that point to the top. For points
    not in the idempotent, reverse the signs.

    """
    data = dict()
    gr_group = BigGradingGroup(pmc)
    for idem in pmc.getIdempotents(idem_size):
        mult = [0] * (pmc.n - 1)
        for p in range(0, pmc.n):
            if pmc.pairid[p] in idem:
                for pos in range(0, p):
                    mult[pos] += Fraction(1, 8)
                for pos in range(p, pmc.n-1):
                    mult[pos] -= Fraction(1, 8)
            else:
                for pos in range(0, p):
                    mult[pos] -= Fraction(1, 8)
                for pos in range(p, pmc.n-1):
                    mult[pos] += Fraction(1, 8)
        data[idem] = BigGradingElement(gr_group, 0, mult)
    return GradingRefinement(pmc, idem_size, data)

# Default function for getting refinement data
DEFAULT_REFINEMENT = averageRefinement

class GradingSet:
    """Represents a general grading set. Can specify an arbitrary number of
    group actions (although only up to two is used).

    """
    def __init__(self, actions):
        """actions is a list of tuples (group, side). side must be either
        ACTION_LEFT or ACTION_RIGHT.

        """
        self.actions = actions
        # Either BIG_GRADING or SMALL_GRADING
        # Take type from the first grading group that acts on it. Other grading
        # groups should have the same type.
        self.type = self.actions[0][0].type

    def multiply(self, set_elt, grp_elt):
        """Returns the result of grp_elt acting on set_elt. grp_elt is a list
        of grading group elements, whose length must equal the number of group
        actions.

        """
        raise NotImplementedError("Group action not implemented.")

    def inverse(self):
        """Returns the inverse of this grading set. Change the side of all
        actions and take the inverse of periodic domains.

        """
        raise NotImplementedError("Set inverse not implemented.")

    def opp(self):
        """Returns the opp of this grading set. Take the opp of all algebra
        actions (changing sides at the same time) and all periodic domains.

        """
        raise NotImplementedError("Set opp not implemented.")

    def Ropp(self):
        """Returns the Ropp of this grading set. Take the opp of all algebra
        actions (changing sides at the same time) and take Ropp of all periodic
        domains.

        """
        raise NotImplementedError("Set Ropp not implemented.")

class GradingSetElement:
    """Represents an element of a grading set."""
    def __init__(self, parent):
        """Specify which grading set this element is in."""
        self.parent = parent

    def __mul__(self, other):
        """Returns the result of group action on self."""
        other = tolist(other)
        if len(other) != len(self.parent.actions):
            return NotImplemented
        for comp in other:
            if not isinstance(comp, GroupElement):
                return NotImplemented
        return self.parent.multiply(self, other)

    def __add__(self, other):
        """Adding an integer / rational number ``n`` means adding the maslov
        component by ``n`` (multiplying by lambda^n).

        """
        if not isinstance(other, Number):
            return NotImplemented
        if len(self.parent.actions) == 0:
            return NotImplemented # Override to take care of this
        to_mult = [gr_group.zero() for gr_group, side in self.parent.actions]
        to_mult[0] = self.parent.actions[0][0].central().power(other)
        return self.parent.multiply(self, to_mult)

    def __sub__(self, other):
        """See description of __add__."""
        return self + (-other)

    def inverse(self):
        """Returns the inverse of this element, in the inverse grading set."""
        raise NotImplementedError("Element inverse not implemented.")

    def opp(self):
        """Returns the opp of this element, in the opp grading set."""
        raise NotImplementedError("Element opp not implemented.")

    def Ropp(self):
        """Returns the Ropp of this element, in the Ropp grading set."""
        raise NotImplementedError("Element Ropp not implemented.")

class SimpleGradingSet(GradingSet):
    """Represents a grading set with one action, that can be simply written as
    a grading group modulo a list of gradings of periodic domains.

    """
    def __init__(self, gr_group, side, periodic_domains):
        """periodic_domains is a list of elements of gr_group. Multiply by one
        of these on the side opposite to the acting side results in the same
        element.

        """
        GradingSet.__init__(self, [(gr_group, side)])
        self.gr_group = gr_group
        self.side = side
        self.periodic_domains = periodic_domains
        # Set up row system for periodic domains - use their spinc parts
        row_vecs = [domain.spinc for domain in self.periodic_domains]
        self.row_sys = RowSystem(row_vecs)

    def multiply(self, set_elt, grp_elt):
        if self.side == ACTION_LEFT:
            new_data = grp_elt[0] * set_elt.data
        else: # self.side == ACTION_RIGHT
            new_data = set_elt.data * grp_elt[0]
        return SimpleGradingSetElement(self, new_data)

    @memorize
    def inverse(self):
        new_domains = [domain.inverse() for domain in self.periodic_domains]
        return SimpleGradingSet(self.gr_group, oppSide(self.side), new_domains)

    @memorize
    def opp(self):
        new_domains = [domain.opp() for domain in self.periodic_domains]
        return SimpleGradingSet(self.gr_group.opp(), oppSide(self.side),
                                new_domains)

    @memorize
    def Ropp(self):
        new_domains = [domain.Ropp() for domain in self.periodic_domains]
        return SimpleGradingSet(self.gr_group.opp(), oppSide(self.side),
                                new_domains)

    def eltEquals(self, elt1, elt2):
        """Equality test for grading elements. Uses the row system for periodic
        domains.

        """
        gr1, gr2 = elt1.data, elt2.data
        m1, m2 = gr1.spinc, gr2.spinc
        spinc_diff = [b-a for a, b in zip(m1, m2)]
        # Get the combination of periodic domains that may give elt2 when
        # multiplied with elt1
        comb = self.row_sys.getComb(spinc_diff)
        if comb is None:
            return False
        # Multiply these periodic domains onto elt1 on the side opposite to the
        # action side.
        for i in range(len(self.periodic_domains)):
            to_mult = self.periodic_domains[i].power(comb[i])
            if self.side == ACTION_LEFT:
                gr1 = gr1 * to_mult
            else:
                gr1 = to_mult * gr1
        # At least the spinc part should now be equal
        assert gr1.spinc == gr2.spinc
        return gr1.maslov == gr2.maslov

    def zero(self):
        """Returns the zero element of this grading set."""
        return SimpleGradingSetElement(self, self.gr_group.zero())

    def simplifiedSet(self):
        """Need a consistent interface with GeneralGradingSet."""
        return self

    def simplifiedElt(self, elt):
        """Need a consistent interface with GeneralGradingSet."""
        return elt

    def __str__(self):
        gr_group, side = self.actions[0]
        result = "One-sided %s grading set with PMC %s on %s.\n" \
            % (grTypeStr(self.type), str(gr_group.pmc), sideStr(side))
        result += "Periodic domains:\n"
        for domain in self.periodic_domains:
            result += "%s\n" % str(domain)
        return result

    def __repr__(self):
        return str(self)

class SimpleGradingSetElement(GradingSetElement):
    """Represents an element of a SimpleGradingSet."""
    def __init__(self, parent, data):
        """In addition to parent set, specify data as an element of the grading
        group.

        """
        GradingSetElement.__init__(self, parent)
        self.data = data
        assert isinstance(data, GroupElement)
        assert data.parent == self.parent.actions[0][0]

    def __eq__(self, other):
        # Equality test need to take into account of periodic domains
        return self.parent.eltEquals(self, other)

    def __ne__(self, other):
        return not (self == other)

    def __hash__(self):
        # Dangerous to use hash since equality test is tricky.
        return NotImplemented

    def __str__(self):
        return str(self.data)

    def __repr__(self):
        return repr(self.data)

    def inverse(self):
        return SimpleGradingSetElement(self.parent.inverse(),
                                       self.data.inverse())

    def opp(self):
        return SimpleGradingSetElement(self.parent.opp(), self.data.opp())

    def Ropp(self):
        return SimpleGradingSetElement(self.parent.Ropp(), self.data.Ropp())

    def simplifiedElt(self):
        return self

class SimpleDbGradingSet(GradingSet):
    """Represents a grading set with two actions, that can be written as a pair
    of grading groups modulo a list of gradings of periodic domains.

    """
    def __init__(self, gr_group1, side1, gr_group2, side2, periodic_domains):
        """periodic_domains is a list of pairs of elements of gr_group.
        Multiply by one of these on the side opposite to the acting side (for
        both actions) results in the same element.

        """
        GradingSet.__init__(self, [(gr_group1, side1), (gr_group2, side2)])
        self.gr_group1, self.side1 = gr_group1, side1
        self.gr_group2, self.side2 = gr_group2, side2
        self.periodic_domains = periodic_domains
        # Set up row system for periodic domains - combine the spinc parts
        row_vecs = [d1.spinc + d2.spinc for d1, d2 in self.periodic_domains]
        self.row_sys = RowSystem(row_vecs)

    def multiply(self, set_elt, grp_elt):
        new_data = [None, None]
        for i in range(len(self.actions)):
            gr_group, side = self.actions[i]
            if side == ACTION_LEFT:
                new_data[i] = grp_elt[i] * set_elt.data[i]
            else: # side == ACTION_RIGHT
                new_data[i] = set_elt.data[i] * grp_elt[i]
        return SimpleDbGradingSetElement(self, new_data)

    def oppMultiply(self, set_elt, grp_elt):
        """Multiply on the side opposite to that of action. Usually used for
        multiplying on periodic domains. grp_elt is a list or tuple of grading
        group elements.

        """
        new_data = [None, None]
        for i in range(len(self.actions)):
            gr_group, side = self.actions[i]
            if side == ACTION_LEFT:
                new_data[i] = set_elt.data[i] * grp_elt[i]
            else: # side == ACTION_RIGHT:
                new_data[i] = grp_elt[i] * set_elt.data[i]
        return SimpleDbGradingSetElement(self, new_data)

    @memorize
    def inverse(self):
        new_domains = [(d2.inverse(), d1.inverse())
                       for d1, d2 in self.periodic_domains]
        return SimpleDbGradingSet(self.gr_group2, oppSide(self.side2),
                                  self.gr_group1, oppSide(self.side1),
                                  new_domains)

    @memorize
    def opp(self):
        new_domains = [(d2.opp(), d1.opp()) for d1, d2 in self.periodic_domains]
        return SimpleDbGradingSet(self.gr_group2.opp(), oppSide(self.side2),
                                  self.gr_group1.opp(), oppSide(self.side1),
                                  new_domains)

    @memorize
    def Ropp(self):
        new_domains = [(d2.Ropp(), d1.Ropp())
                       for d1, d2 in self.periodic_domains]
        return SimpleDbGradingSet(self.gr_group2.opp(), oppSide(self.side2),
                                  self.gr_group1.opp(), oppSide(self.side1),
                                  new_domains)

    @memorize
    def partialRopp(self, action_id):
        """Take Ropp on one of the actions."""
        if action_id == 0:
            new_domains = [(d2, d1.Ropp()) for d1, d2 in self.periodic_domains]
            return SimpleDbGradingSet(self.gr_group2, self.side2,
                self.gr_group1.opp(), oppSide(self.side1), new_domains)
        else:
            assert action_id == 1
            new_domains = [(d1, d2.Ropp()) for d1, d2 in self.periodic_domains]
            return SimpleDbGradingSet(
                self.gr_group1, self.side1,
                self.gr_group2.opp(), oppSide(self.side2), new_domains)

    def eltEquals(self, elt1, elt2):
        """Equality test for grading elements. Uses the row system for periodic
        domains.

        """
        (gr11, gr12), (gr21, gr22) = elt1.data, elt2.data
        m1 = gr11.spinc + gr12.spinc
        m2 = gr21.spinc + gr22.spinc
        spinc_diff = [b-a for a, b in zip(m1, m2)]
        # Get the combination of periodic domains that may give elt2 when
        # multiplied with elt1
        # Use the more lenient form of comparison, allowing rational multiple
        # of periodic domains
        comb = self.row_sys.getComb(spinc_diff, use_rational = True)
        if comb is None:
            return False
        # Multiply these periodic domains onto elt1 on the side opposite to the
        # action side.
        for i in range(len(self.periodic_domains)):
            to_mult = [comp.power(comb[i])
                       for comp in self.periodic_domains[i]]
            elt1 = self.oppMultiply(elt1, to_mult)
        # At least the spinc part should now be equal
        assert elt1.data[0].spinc == elt2.data[0].spinc
        assert elt1.data[1].spinc == elt2.data[1].spinc
        # Check the (sum of) maslov part is equal
        return elt1.data[0].maslov + elt1.data[1].maslov == \
            elt2.data[0].maslov + elt2.data[1].maslov

    def shiftRight(self, elt):
        """Use the periodic domains to change elt into a form where the first
        component is zero. Returns the new element if this is possible, and
        None if otherwise.

        """
        spinc = elt.data[0].spinc + elt.data[1].spinc
        comb, reduced_vec = self.row_sys.vecReduce(spinc, use_rational = False)
        if comb is None:
            return None
        for i in range(len(self.periodic_domains)):
            to_mult = [comp.power(-comb[i])
                       for comp in self.periodic_domains[i]]
            elt = self.oppMultiply(elt, to_mult)
        if not all([n == 0 for n in elt.data[0].spinc]):
            return None
        elt.data[1].maslov += elt.data[0].maslov
        elt.data[0].maslov = 0
        return elt

    def zero(self):
        """Returns the zero element of this grading set."""
        return SimpleDbGradingSetElement(self, (self.gr_group1.zero(),
                                                self.gr_group2.zero()))

    def __str__(self):
        result = "Two-sided %s grading set with PMC %s on %s and %s on %s.\n" \
            % (grTypeStr(self.type),
               str(self.gr_group1.pmc), sideStr(self.side1),
               str(self.gr_group2.pmc), sideStr(self.side2))
        result += "Periodic domains:\n"
        for domain in self.periodic_domains:
            result += "%s\n" % strDbGrading(domain)
        return result

    def isAutomorphism(self):
        """Test whether the periodic domains of this grading set determine an
        automorphism of grading groups. This is true if the number of periodic
        domains agrees with the dimension of spinc component (of one side), and
        that the square matrix formed by the periodic domains on each side is
        invertible as an integer matrix.

        Being an automorphism means this grading set can be tensored to another
        without increasing its size.

        """
        for i in range(self.gr_group1.spinc_len):
            elt = SimpleDbGradingSetElement(
                self, [self.gr_group1.basis(i), self.gr_group2.zero()])
            if self.shiftRight(elt) is None:
                return False
        return True

    def simplifiedSet(self):
        return self

    def simplifiedElt(self, elt):
        return elt

    def __repr__(self):
        return str(self)

class SimpleDbGradingSetElement(GradingSetElement):
    """Represents an element of a SimpleDbGradingSet."""
    def __init__(self, parent, data):
        """In addition to parent set, specify data as a pair of elements of
        grading groups.

        """
        GradingSetElement.__init__(self, parent)
        self.data = data
        assert len(data) == 2
        assert isinstance(data[0], GroupElement)
        assert isinstance(data[1], GroupElement)
        assert data[0].parent == self.parent.actions[0][0]
        assert data[1].parent == self.parent.actions[1][0]

    def __eq__(self, other):
        # Equality test need to take into account of periodic domains
        return self.parent.eltEquals(self, other)

    def __ne__(self, other):
        return not (self == other)

    def __hash__(self):
        # Dangerous to use hash since equality test is tricky.
        return NotImplemented

    def __str__(self):
        return strDbGrading(self.data)

    def __repr__(self):
        return str(self)

    def inverse(self):
        d1, d2 = self.data
        return SimpleDbGradingSetElement(self.parent.inverse(),
                                         (d2.inverse(), d1.inverse()))

    def opp(self):
        d1, d2 = self.data
        return SimpleDbGradingSetElement(self.parent.opp(),
                                         (d2.opp(), d1.opp()))

    def Ropp(self):
        d1, d2 = self.data
        return SimpleDbGradingSetElement(self.parent.Ropp(),
                                         (d2.Ropp(), d1.Ropp()))

    def partialRopp(self, action_id):
        """Take Ropp on one of the actions."""
        d1, d2 = self.data
        if action_id == 0:
            return SimpleDbGradingSetElement(self.parent.partialRopp(0),
                                             (d2, d1.Ropp()))
        else:
            assert action_id == 1
            return SimpleDbGradingSetElement(self.parent.partialRopp(1),
                                             (d1, d2.Ropp()))

    def shiftRight(self):
        """See the function of same name in SimpleDbGradingSet for details."""
        return self.parent.shiftRight(self)

    def simplifiedElt(self):
        return self

class GeneralGradingSet(GradingSet):
    """Represents a general grading set formed by tensoring together several
    simple grading sets. The right action of each component should match the
    left action of the next component. The first action of the first component
    and the second action of the last component, if they have two actions, form
    the group action of the compound grading set.

    """
    def __init__(self, comps):
        """Initialize using the list of component grading sets. Group actions
        on this compound grading set can be deduced from these.
        Entries in comps can be GeneralGradingSet themselves. In this case it
        is expanded into its components.

        """
        # Process comps by expanding possible GeneralGradingSet into components
        processed_comps = []
        for comp in comps:
            if isinstance(comp, GeneralGradingSet):
                processed_comps += comp.comps
            else:
                processed_comps.append(comp)
        comps = processed_comps

        self.comps = comps
        self.num_comp = len(self.comps)
        self.actions = []
        self.action_comp_id = []
        # Middle components must be LEFT-RIGHT
        def assertLeftRight(gr_set):
            "Assert gr_set is a LEFT-RIGHT double grading set."""
            assert isinstance(gr_set, SimpleDbGradingSet)
            assert gr_set.side1 == ACTION_LEFT and gr_set.side2 == ACTION_RIGHT
        for i in range(1, self.num_comp-1):
            assertLeftRight(self.comps[i])

        # First component is either ?-RIGHT or RIGHT
        if isinstance(self.comps[0], SimpleGradingSet):
            assert self.comps[0].side == ACTION_RIGHT
        else:
            assert self.comps[0].side2 == ACTION_RIGHT
            self.actions.append((self.comps[0].gr_group1, self.comps[0].side1))
            self.action_comp_id.append((0, 0))
        # Last component is either LEFT-? or LEFT
        if isinstance(self.comps[-1], SimpleGradingSet):
            assert self.comps[-1].side == ACTION_LEFT
        else:
            assert self.comps[-1].side1 == ACTION_LEFT
            self.actions.append((self.comps[-1].gr_group2,
                                 self.comps[-1].side2))
            self.action_comp_id.append((-1, 1))

        # Initialize GradingSet
        self.type = self.comps[0].type

        # Set up row system
        # Each entry of row_action shows the multiplication actions to take
        # when the corresponding row should be added. It is a list of
        # quadruples (ID1, ID2, gr, side), where ID1 is the ID within gr_sets,
        # ID2 is either 0 or 1 giving index within a gr_set (for single grading
        # sets it is always 0). gr is the big/small grading element to be
        # multiplied. side is the side it is multiplied to.

        self.row_action = []
        # Add rows corresponding to periodic domains
        for i in range(self.num_comp):
            gr_set = self.comps[i]
            for domain in gr_set.periodic_domains:
                if i == 0 and isinstance(domain, GroupElement):
                    cur_action = [(0, 0, domain, ACTION_LEFT)]
                elif i == self.num_comp-1 and isinstance(domain, GroupElement):
                    cur_action = [(self.num_comp-1, 0, domain, ACTION_RIGHT)]
                else:
                    cur_action = [(i, 0, domain[0], oppSide(gr_set.side1)),
                                  (i, 1, domain[1], oppSide(gr_set.side2))]
                self.row_action.append(cur_action)
        # Add rows corresponding to tensoring
        for comp_id in range(self.num_comp-1):
            if isinstance(self.comps[comp_id], SimpleGradingSet):
                spinc_len = self.comps[comp_id].gr_group.spinc_len
            else:
                spinc_len = self.comps[comp_id].gr_group2.spinc_len
            for basis_id in range(spinc_len):
                cur_action = self._getTensorAction(comp_id, basis_id)
                self.row_action.append(cur_action)

        row_vecs = [self._vecFromAction(action) for action in self.row_action]
        self.row_sys = RowSystem(row_vecs)

        # Row systems that exclude periodic domains on single grading sets.
        # Useful for simplifying grading elements (see simplifiedElt())
        self.row_action_db_only = [action for action in self.row_action
                                   if len(action) == 2]
        row_vecs_db_only = [self._vecFromAction(action)
                            for action in self.row_action_db_only]
        self.row_sys_db_only = RowSystem(row_vecs_db_only)

        # Find list of zero combinations of rows (these correspond to
        # provincial periodic domains, that may introduce ambiguity in the
        # maslov component.
        self.zero_combs = self.row_sys.getZeroComb()

    def _vecFromAction(self, vec_action):
        """Get row vector for a given action. Basically copy the spinc
        components to the right locations.

        """
        start_pos = dict()
        cur_pos = 0
        for set_id in range(self.num_comp):
            gr_set = self.comps[set_id]
            for action_id in range(len(gr_set.actions)):
                gr_group, side = gr_set.actions[action_id]
                cur_len = len(gr_group.zero().spinc)
                start_pos[(set_id, action_id)] = cur_pos
                cur_pos += cur_len
        total_len = cur_pos
        vec = [0] * total_len
        for set_id, action_id, elt, side in vec_action:
            spinc = elt.spinc
            start = start_pos[(set_id, action_id)]
            for i in range(len(spinc)):
                vec[start+i] += spinc[i]
        return vec

    def _getTensorAction(self, comp_id, basis_id):
        if comp_id == 0 and len(self.comps[0].actions) == 1:
            gr_group = self.comps[0].gr_group
            action_id1 = 0
        else:
            gr_group = self.comps[comp_id].gr_group2
            action_id1 = 1
        basis = gr_group.basis(basis_id)
        return [(comp_id, action_id1, basis, ACTION_RIGHT),
                (comp_id+1, 0, basis.inverse(), ACTION_LEFT)]

    def multiply(self, set_elt, grp_elt):
        elt_data = list(set_elt.data)
        for i in range(len(self.actions)):
            comp, action_id = self.action_comp_id[i]
            # Note group action for set elements are always written on the
            # right in python.
            if action_id == 0:
                if len(self.comps[comp].actions) == 1:
                    elt_data[comp] = elt_data[comp] * grp_elt[i]
                else:
                    to_mult = [grp_elt[i],
                               self.comps[comp].actions[1][0].zero()]
                    elt_data[comp] = elt_data[comp] * to_mult
            else: # action_id == 1
                to_mult = [self.comps[comp].actions[0][0].zero(), grp_elt[i]]
                elt_data[comp] = elt_data[comp] * to_mult
        return GeneralGradingSetElement(self, elt_data)

    @memorize
    def inverse(self):
        return GeneralGradingSet(
            reversed([comp.inverse() for comp in self.comps]))

    @memorize
    def opp(self):
        return GeneralGradingSet(reversed([comp.opp() for comp in self.comps]))

    @memorize
    def Ropp(self):
        return GeneralGradingSet(
            reversed([comp.Ropp() for comp in self.comps]))

    @memorize
    def partialRopp(self, action_id):
        """Take Ropp on one of the actions."""
        assert len(self.actions) == 2
        if action_id == 0:
            raise NotImplementedError("partialRopp on first not implemented.")
        else:
            assert action_id == 1
            return GeneralGradingSet(self.comps[0:-1] + \
                                         [self.comps[-1].partialRopp(1)])

    def _performActionListForm(self, elt_list, comb, db_only = False):
        """Perform the set of row actions given by comb on elt_list (list form
        given by elt.listForm(). Return the resulting list (original list may
        also be changed).

        If db_only is set to True, use the db_only row system (exclude periodic
        domains on single grading sets).

        """
        if db_only:
            action_list = self.row_action_db_only
        else:
            action_list = self.row_action
        for i in range(len(action_list)):
            vec_action = action_list[i]
            for set_id, action_id, elt, side in vec_action:
                to_mult = elt.power(comb[i])
                if side == ACTION_LEFT:
                    elt_list[set_id][action_id] = \
                        to_mult * elt_list[set_id][action_id]
                else: # side == ACTION_RIGHT
                    elt_list[set_id][action_id] = \
                        elt_list[set_id][action_id] * to_mult
        return elt_list

    def eltEquals(self, elt1, elt2):
        """Equality test for grading elements. Use row system to identify
        actions to take and use self.row_action to carry them out.

        """
        # Use more lenient form of comparison for double grading sets
        if len(self.actions) == 2:
            use_rational = True
        else:
            use_rational = False
        maslov_diff, mod, spinc_diff = \
            self.eltDiffShortForm(elt1, elt2, use_rational)
        if spinc_diff is None:
            return False
        if not all([n == 0 for n in spinc_diff]):
            return False
        if mod == 0:
            return maslov_diff == 0
        if isinstance(maslov_diff, Fraction):
            if maslov_diff.denominator != 1: return False
            return maslov_diff.numerator % mod == 0
        else: # maslov_diff is int
            return maslov_diff % mod == 0

    def eltDiffShortForm(self, elt1, elt2, use_rational = False):
        """Returns the difference (maslov, mod, spinc) between two elements in
        short form. Measures the grading of elt1 with reference to elt2
        (not the other way around!).

        """
        elt1.clearFirst()
        elt2.clearFirst()
        # Form vectors corresponding to the two elements
        m1 = sum([comp.spinc for comp in flatten(elt1.listForm())], [])
        m2 = sum([comp.spinc for comp in flatten(elt2.listForm())], [])
        spinc_diff = [b-a for a, b in zip(m1, m2)]

        # Get the combination of row actions that give elt2, up to a difference
        # in standard form, when performed on elt1
        if not use_rational:
            for n in spinc_diff:
                if isinstance(n, Fraction) and n.denominator != 1:
                    return (0, 0, None)
        comb, reduced_vec = self.row_sys.vecReduce(spinc_diff, use_rational)

        # Perform these actions
        elt1_list = self._performActionListForm(elt1.listForm(), comb)

        # Find and compare the sum of maslov parts.
        elt2_list = elt2.listForm()
        sum1 = sum(comp.maslov for comp in flatten(elt1_list))
        sum2 = sum(comp.maslov for comp in flatten(elt2_list))
        maslov_diff = sum1 - sum2

        # Are there zero combinations that could introduce an ambiguity?
        cur_gcd = 0
        for zero_comb in self.zero_combs:
            elt1_list_after = self._performActionListForm(elt1_list, zero_comb)
            sum_after = sum(comp.maslov for comp in flatten(elt1_list_after))
            cur_gcd = gcd(cur_gcd, abs(sum_after - sum1))
        mod = cur_gcd

        spinc_short = self.row_sys.shortForm(reduced_vec, use_rational)
        return (maslov_diff, mod, [-n for n in spinc_short])

    def eltAbsoluteGrading(self, elt):
        """Returns the absolute grading of element. Element must be in the zero
        spinc class.

        """
        maslov, mod, spinc_short = \
            self.eltDiffShortForm(elt, self.zero(), use_rational = True)
        assert mod == 0 and spinc_short is not None and \
            all([n == 0 for n in spinc_short])
        return maslov

    def zero(self):
        zero_data = [gr_set.zero() for gr_set in self.comps]
        return GeneralGradingSetElement(self, zero_data)

    @memorize
    def canSimplify(self):
        """Returns whether the general grading set can be simplified down to one
        SimpleGradingSet or SimpleDbGradingSet. If the grading set has one
        action, all of its double components must be automorphisms. If the
        grading set has two actions, all but the first component must be
        automorphisms.

        """
        assert len(self.actions) > 0
        if isinstance(self.comps[0], SimpleGradingSet):
            return all([gr_set.isAutomorphism() for gr_set in self.comps[1:]])
        elif isinstance(self.comps[-1], SimpleGradingSet):
            return all([gr_set.isAutomorphism() for gr_set in self.comps[0:-1]])
        else:
            return all([gr_set.isAutomorphism() for gr_set in self.comps[1:]])
        
    @memorize
    def simplifiedSet(self):
        """Returns an equivalent SimpleGradingSet or SimpleDbGradingSet using
        the automorphism property of its double components.

        """
        assert len(self.actions) > 0
        new_pdomains = []
        if not self.canSimplify():
            return self

        if isinstance(self.comps[0], SimpleGradingSet):
            for p_domain in self.comps[0].periodic_domains:
                elt = p_domain
                for i in range(1, len(self.comps)):
                    elt = SimpleDbGradingSetElement(
                        self.comps[i], [elt, self.comps[i].gr_group2.zero()])
                    elt = elt.shiftRight()
                    elt = elt.data[1]
                new_pdomains.append(elt)
            gr_group, side = self.comps[-1].actions[1]
            return SimpleGradingSet(gr_group, side, new_pdomains)
        elif isinstance(self.comps[-1], SimpleGradingSet):
            return self.inverse().simplifiedSet().inverse()
        else: # Double grading set
            for p_domain in self.comps[0].periodic_domains:
                elt1, elt2 = p_domain
                for i in range(1, len(self.comps)):
                    elt2 = SimpleDbGradingSetElement(
                        self.comps[i], [elt2, self.comps[i].gr_group2.zero()])
                    elt2 = elt2.shiftRight()
                    elt2 = elt2.data[1]
                new_pdomains.append([elt1, elt2])
            gr_group1, side1 = self.comps[0].actions[0]
            gr_group2, side2 = self.comps[-1].actions[1]
            return SimpleDbGradingSet(
                gr_group1, side1, gr_group2, side2, new_pdomains)

    def simplifiedElt(self, elt):
        """Returns the version of elt in simplified form of this set."""
        assert len(self.actions) > 0
        if not self.canSimplify():
            return elt

        if isinstance(self.comps[-1], SimpleGradingSet):
            return self.inverse().simplifiedElt(elt.inverse()).inverse()

        spinc = sum([comp.spinc for comp in flatten(elt.listForm())], [])
        comb, reduced_vec = self.row_sys_db_only.vecReduce(
            spinc, use_rational = True)
        comb = [-n for n in comb]
        elt_list = flatten(
            self._performActionListForm(elt.listForm(), comb, db_only = True))
        prev_maslov_sum = 0
        for comp in elt_list[:-1]:
            assert all([n == 0 for n in comp.spinc])
            prev_maslov_sum += comp.maslov
        elt_list[-1].maslov += prev_maslov_sum
        simpl_set = self.simplifiedSet()
        if isinstance(self.comps[0], SimpleGradingSet):
            return SimpleGradingSetElement(simpl_set, elt_list[-1])
        else:
            return SimpleDbGradingSetElement(
                simpl_set, [simpl_set.gr_group1.zero(), elt_list[-1]])

    def __str__(self):
        result = "Composite grading set.\n"
        result += "\n".join([str(gr_set) for gr_set in self.comps])
        return result

    def __repr__(self):
        return str(self)

class GeneralGradingSetElement(GradingSetElement):
    """Represents an element of a GeneralGradingSet."""
    def __init__(self, parent, data):
        """In addition to parent set, specify data as a list of
        SimpleGradingSetElement or SimpleDbGradingSetElement, one for each
        element of parent.gr_sets.
        Entries in data can be GeneralGradingSetElement themselves. In this
        case it is expanded into its components.

        """
        GradingSetElement.__init__(self, parent)
        processed_data = []
        for comp in data:
            if isinstance(comp, GeneralGradingSetElement):
                processed_data += comp.data
            else:
                processed_data.append(comp)
        data = processed_data
        for i in range(len(data)):
            assert data[i].parent == self.parent.comps[i]
        self.data = data

    def clearFirst(self):
        """Clear out the first part of grading in all components starting at
        the second, using the tensor relations. This is useful in avoiding
        fractions when comparing elements and putting their differences in
        short form.

        """
        for i in range(len(self.data)-1):
            if isinstance(self.data[i+1], SimpleDbGradingSetElement):
                data = self.data[i+1].data[0]
            else:
                data = self.data[i+1].data

            if isinstance(self.data[i], SimpleDbGradingSetElement):
                to_mult = [self.data[i].parent.gr_group1.zero(), data]
            else:
                to_mult = [data]
            self.data[i] = self.data[i] * to_mult

            if isinstance(self.data[i+1], SimpleDbGradingSetElement):
                to_mult = [data.inverse(),
                           self.data[i+1].parent.gr_group2.zero()]
            else:
                to_mult = [data.inverse()]
            self.data[i+1] = self.data[i+1] * to_mult

    def listForm(self):
        """Present grading data in this element as a list of lists of grading
        group elements, one list for each Simple(Db)GradingSetElement.

        """
        result = []
        for elt in self.data:
            if isinstance(elt.data, GroupElement):
                result.append([elt.data])
            else:
                result.append(list(elt.data))
        return result

    def __add__(self, other):
        # Overridden to take care of the case of no group actions
        if not isinstance(other, Number):
            return NotImplemented
        return GeneralGradingSetElement(self.parent,
                                        [self.data[0] + other] + self.data[1:])

    def __eq__(self, other):
        # Equality test need to take into account of tensor product as well as
        # periodic domains.
        return self.parent.eltEquals(self, other)

    def __ne__(self, other):
        return not (self == other)

    def __hash__(self):
        # Dangerous to use hash since equality test is tricky.
        return NotImplemented

    def __str__(self):
        return str(self.data)

    def __repr__(self):
        return repr(self.data)

    def inverse(self):
        return GeneralGradingSetElement(
            self.parent.inverse(), reversed([d.inverse() for d in self.data]))

    def opp(self):
        return GeneralGradingSetElement(
            self.parent.opp(), reversed([d.opp() for d in self.data]))

    def Ropp(self):
        return GeneralGradingSetElement(
            self.parent.Ropp(), reversed([d.Ropp() for d in self.data]))

    def partialRopp(self, action_id):
        """Take Ropp on one of the actions."""
        assert len(self.parent.actions) == 2
        if action_id == 0:
            raise NotImplementedError("partialRopp on first not implemented.")
        else:
            assert action_id == 1
            return GeneralGradingSetElement(
                self.parent.partialRopp(1),
                self.data[0:-1] + [self.data[-1].partialRopp(1)])

    def simplifiedElt(self):
        """Returns the equivalent element in the simplified version of parent.

        """
        return self.parent.simplifiedElt(self)

def strDbGrading(data):
    """data is a pair of gradings. Print in the [maslov; spinc1; spinc2]
    format.

    """
    total_maslov = data[0].maslov + data[1].maslov
    spinc1, spinc2 = data[0].spinc, data[1].spinc
    return "[%s; %s; %s]" % (str(total_maslov),
                             ", ".join(str(n) for n in spinc1),
                             ", ".join(str(n) for n in spinc2))
