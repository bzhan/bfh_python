"""This module offers minimal support for PMC with boundaries and unmatched
points. A normal PMC (defined in pmc.py) can be split into one or several
local PMC's like this.

"""

from algebra import DGAlgebra, Element, Generator
from algebra import E0
from pmc import Strands, StrandDiagram
from utility import memorize, memorizeHash, subset
from utility import F2
import itertools

class LocalPMC:
    """Represents a pointed matched circle with boundaries and unmatched
    points.

    """
    def __init__(self, n, matching, endpoints):
        """Creates a pointed matched circle with n points (including
        endpoints).

        - matching: a list of tuples, with each tuple containing either one
        (for unpaired point) or two (for matched pair) points.
        - endpoints: list of endpoints. Must be disjoint from numbers involved
        in matching.

        """
        self.n = n
        self.endpoints = endpoints

        # Total number of pairs and single points.
        self.num_pair = len(matching)
        # otherp[i] is the point paired to i. Equals -1 for boundary point, and
        # i for unpaired points.
        self.otherp = [-1] * self.n
        for pair in matching:
            if len(pair) == 2:
                self.otherp[pair[0]] = pair[1]
                self.otherp[pair[1]] = pair[0]
            else:  # len(pair) == 1
                self.otherp[pair[0]] = pair[0]
        for endpoint in endpoints:
            assert self.otherp[endpoint] == -1
        # pairid[i] is the ID of the pair containing i.
        # 0 <= pairid[i] < num_pair if i is not an endpoint.
        # pairid[i] = -1 if i is an endpoint.
        self.pairid = [-1] * self.n
        # pairs[i] is the pair of points with ID i (0 <= i < num_pair)
        self.pairs = []
        pair_count = 0
        for pos in range(self.n):
            if pos in endpoints: continue
            if self.pairid[pos] == -1:
                self.pairid[pos] = self.pairid[self.otherp[pos]] = pair_count
                if pos == self.otherp[pos]:
                    self.pairs.append((pos,))
                else:
                    self.pairs.append((pos, self.otherp[pos]))
                pair_count += 1
        assert pair_count == self.num_pair

    def __eq__(self, other):
        return self.otherp == other.otherp

    def __ne__(self, other):
        return not (self == other)

    def __hash__(self):
        return hash(tuple(self.otherp))

    def __str__(self):
        return str(self.pairs)

    def __repr__(self):
        return "LocalPMC(%s)" % str(self)

    def sd(self, data):
        """Simple way to obtain a local strand diagram for this local PMC. Each
        element of data is either an integer or a pair. An integer specifies a
        single or double horizontal at this position (and its paired position,
        if any). A pair (p, q) specifies a strand from p to q.

        """
        parent = self.getAlgebra()
        left_idem = []
        strands = []
        for d in data:
            if isinstance(d, int):
                assert self.pairid[d] != -1
                left_idem.append(self.pairid[d])
            else:
                if self.pairid[d[0]] != -1:
                    left_idem.append(self.pairid[d[0]])
                strands.append(d)
        return LocalStrandDiagram(parent, left_idem, strands)

    def getAlgebra(self):
        """Returns the local strand algebra for this local PMC."""
        return LocalStrandAlgebra(F2, self)

    def getSingleIdems(self):
        """Return the list of indices of pairs / single points that are single
        points.

        """
        return [i for i in range(self.num_pair) if len(self.pairs[i]) == 1]

    def getIdempotents(self):
        """Get the list of all idempotents (no restriction on size)."""
        return [LocalIdempotent(self, data) for sz in range(self.num_pair+1)
                for data in itertools.combinations(range(self.num_pair), sz)]

    def getStrandDiagrams(self):
        """Returns the list of generators of the local strand algebra. Note we
        automatically impose the multiplicity-one condition, and there are no
        constraints on the size of idempotents.

        """
        algebra = self.getAlgebra()
        result = []
        def search(cur_strands):
            """Search starting with the given list of strands. May only add
            strands after the end position of the last strand.

            """
            # First, check if the current list of strands is valid.
            left_occupied = [0] * self.num_pair
            right_occupied = [0] * self.num_pair
            for start, end in cur_strands:
                start_id, end_id = self.pairid[start], self.pairid[end]
                if start_id != -1:
                    left_occupied[start_id] += 1
                if end_id != -1:
                    right_occupied[end_id] += 1
            if any([n >= 2 for n in left_occupied + right_occupied]):
                # There should not be two strands starting or ending at points
                # in the same pair.
                return

            # Enumerate all possible ways of adding idempotents.
            empty_idems = [i for i in range(self.num_pair)
                           if left_occupied[i] == 0 and right_occupied[i] == 0]
            left_idem = [i for i in range(self.num_pair)
                         if left_occupied[i] > 0]
            for idems_to_add in subset(empty_idems):
                result.append(LocalStrandDiagram(
                    algebra, left_idem + list(idems_to_add), cur_strands))

            # Now enumerate all ways of adding more strands.
            last_end = 0
            if len(cur_strands) > 0:
                last_end = cur_strands[-1][1]
            for start in range(last_end, self.n):
                for end in range(start + 1, self.n):
                    if self.pairid[start] == -1 and self.pairid[end] == -1 and \
                       end == start + 1:
                        # Exclude cases where a strand goes from an
                        # end-boundary-point to a start-boundary-point
                        break
                    search(cur_strands + [(start, end)])
                    if self.pairid[end] == -1:
                        # No strand should go beyond an end-boundary-point.
                        break
        search([])
        return result
                
class LocalIdempotent(tuple):
    """Represents a local idempotent in a certain local PMC. Stored as a tuple
    of occupied pairs.

    """
    def __new__(cls, local_pmc, data):
        return tuple.__new__(cls, tuple(sorted(data)))

    def __init__(self, local_pmc, data):
        self.local_pmc = local_pmc

    def __eq__(self, other):
        if isinstance(other, LocalIdempotent):
            return self.local_pmc == other.local_pmc and \
                tuple.__eq__(self, other)
        else:
            return False

    def __ne__(self, other):
        return not (self == other)

    def __hash__(self):
        return hash((self.local_pmc, tuple(self), "LocalIdempotent"))

    def __str__(self):
        return repr(self)

    def __repr__(self):
        return "(%s)" % ",".join(str(self.local_pmc.pairs[i]) for i in self)

    def removeSingleHor(self, idems = None):
        """Returns a new idempotent with single points removed.

        If idems is None (default case): all single idempotents will be removed.
        Otherwise, idems must be a list containing single idems to be removed.

        """
        return LocalIdempotent(
            self.local_pmc, [i for i in self
                             if len(self.local_pmc.pairs[i]) == 2 or
                             (idems is not None and i not in idems)])

    def toAlgElt(self):
        """Get the local strand diagram corresponding to this idempotent (the
        strand algebra is uniquely specified from the local PMC.

        """
        return LocalStrandDiagram(self.local_pmc.getAlgebra(), self, [])

class LocalStrands(tuple):
    """Represents a fixed list of strands in a local PMC. Stored as a tuple of
    pairs.

    """
    def __new__(cls, local_pmc, data):
        return tuple.__new__(cls, tuple(sorted(data)))

    def __init__(self, local_pmc, data):
        self.local_pmc = local_pmc
        # Compute multiplicity. multiplicity[i] represents the multiplicity
        # on the interval (i, i+1). Intervals that are gaps between two
        # boundary points should never be occupied, but we do not check for it
        # here.
        self.multiplicity = [0] * (self.local_pmc.n - 1)
        for st in self:
            assert len(st) == 2 and st[0] < st[1]
            for pos in range(st[0], st[1]):
                self.multiplicity[pos] += 1

    def propagateRight(self, left_idem):
        """Find the right_idem given left_idem and strand info. This is similar
        to propagateRight in pmc.py, except we should note that idempotents can
        appear or disappear due to strands going into and out of the boundary.

        """
        pmc = self.local_pmc
        idem_count = [0] * pmc.num_pair
        for pair in left_idem:
            idem_count[pair] += 1
        for st in self:
            if st[0] not in pmc.endpoints:
                if idem_count[pmc.pairid[st[0]]] == 0: return None
                idem_count[pmc.pairid[st[0]]] -= 1
        for st in self:
            if st[1] not in pmc.endpoints:
                if idem_count[pmc.pairid[st[1]]] == 1: return None
                idem_count[pmc.pairid[st[1]]] += 1
        right_idem = [i for i in range(pmc.num_pair) if idem_count[i] == 1]
        return LocalIdempotent(pmc, right_idem)

class LocalStrandDiagram(Generator):
    """A strand diagram in an local PMC.

    Multiplicity-one condition is hard-coded in. Dealing with multiplicity
    greater than one when there are strands going to the boundary can be more
    subtle.

    """
    def __init__(self, parent, left_idem, strands):
        """Specifies the parent algebra (which contains the local PMC), left
        idempotent, and strands.

        Input parameters:
        - parent: must be an object of LocalStrandAlgebra.
        - left_idem: tuple containing IDs of occupied pairs.
        - strands: tuple of pairs specifying strands.

        """
        Generator.__init__(self, parent)
        self.local_pmc = parent.local_pmc
        self.left_idem = left_idem
        if not isinstance(self.left_idem, LocalIdempotent):
            self.left_idem = LocalIdempotent(self.local_pmc, self.left_idem)
        if not isinstance(strands, LocalStrands):
            strands = LocalStrands(self.local_pmc, strands)
        self.strands = strands

        # Get right_idem and multiplicity from strands
        self.right_idem = self.strands.propagateRight(self.left_idem)
        self.multiplicity = self.strands.multiplicity

        # Enumerate single and double horizontals
        self.all_hor = list(self.left_idem)
        for st in self.strands:
            start_idem = self.local_pmc.pairid[st[0]]
            if start_idem != -1:
                if start_idem not in self.all_hor:
                    print left_idem, strands
                self.all_hor.remove(start_idem)
        self.all_hor = tuple(self.all_hor)
        self.single_hor = tuple([i for i in self.all_hor
                                 if len(self.local_pmc.pairs[i]) == 1])
        self.double_hor = tuple([i for i in self.all_hor
                                 if len(self.local_pmc.pairs[i]) == 2])

    def isIdempotent(self):
        """Tests whether this generator is an idempotent."""
        return len(self.strands) == 0

    def getLeftIdem(self):
        """Return the left idempotent."""
        return self.left_idem

    def getRightIdem(self):
        """Return the right idempotent."""
        return self.right_idem

    def __str__(self):
        return "[%s]" % \
            ",".join([str(self.local_pmc.pairs[i]) for i in self.all_hor] +
                     ["%s->%s" % (p, q) for (p, q) in self.strands])

    def __repr__(self):
        return str(self)

    def inputForm(self):
        """Output in the form used for sd."""
        return "[%s]" % \
            ", ".join([str(self.local_pmc.pairs[i][0]) for i in self.all_hor] +
                      ["(%s, %s)" % (p, q) for (p, q) in self.strands])

    def __eq__(self, other):
        return self.parent == other.parent and \
            self.left_idem == other.left_idem and \
            self.strands == other.strands

    def __ne__(self, other):
        return not (self == other)

    @memorizeHash
    def __hash__(self):
        return hash((self.parent, tuple(self.left_idem),
                     tuple(self.strands)))

    @memorize
    def removeSingleHor(self, idems = None):
        """Return a local strand diagram that is just like this, except with
        some single horizontal lines removed.

        If idems is None (default case): all single idempotents will be removed.
        Otherwise, idems must be a list containing single idems to be removed.

        """
        new_left_idem = list(self.left_idem)
        for i in self.single_hor:
            if idems is None or i in idems:
                new_left_idem.remove(i)
        return LocalStrandDiagram(self.parent, new_left_idem, self.strands)

    @memorize
    def addSingleHor(self, idems):
        """Opposite of removeSingleHor. Add some single idempotents. The
        conditions on inputs are more strict than in removeSingleHor. idems must
        be provided and must be idempotents not already in self.

        """
        new_left_idem = list(self.left_idem)
        for i in idems:
            assert i not in self.left_idem
            assert len(self.parent.local_pmc.pairs[i]) == 1
            new_left_idem.append(i)
        return LocalStrandDiagram(self.parent, new_left_idem, self.strands)

# Don't need anything beyond Element at this time.
LocalStrandDiagram.ELT_CLASS = Element

class LocalStrandAlgebra(DGAlgebra):
    """Represents the strand algebra of a local PMC."""

    def __init__(self, ring, local_pmc):
        """Specifies the local PMC. Note that unlike StrandAlgebra, the size of
        idempotent is variable, and we only implement the multiplicity one case.

        """
        DGAlgebra.__init__(self, ring)
        self.local_pmc = local_pmc

    def __str__(self):
        return "Local strand algebra over %s" % str(self.local_pmc)

    def __eq__(self, other):
        return self.local_pmc == other.local_pmc

    def __ne__(self, other):
        return not (self == other)

    @memorizeHash
    def __hash__(self):
        return hash(("LocalStrandAlgebra", self.local_pmc))

    @memorize
    def diff(self, gen):
        cur_strands = gen.strands
        result = E0
        # In multiplicity one case, only need to worry about uncrossing a moving
        # strand with a horizontal. Also, no need to worry about double
        # crossing.
        for st in cur_strands:
            for i in gen.all_hor:
                for p in gen.local_pmc.pairs[i]:
                    if st[0] <= p and p <= st[1]:
                        new_strands = list(cur_strands)
                        new_strands.remove(st)
                        new_strands.extend([(st[0], p), (p, st[1])])
                        result += LocalStrandDiagram(
                            self, gen.left_idem, new_strands).elt()
        return result

    @memorize
    def getGenerators(self):
        return self.local_pmc.getStrandDiagrams()

    def multiply(self, gen1, gen2):
        if not isinstance(gen1, LocalStrandDiagram):
            return NotImplemented
        if not isinstance(gen2, LocalStrandDiagram):
            return NotImplemented
        assert gen1.parent == self and gen2.parent == self, \
            "Algebra not compatible."

        if gen1.right_idem != gen2.left_idem:
            return E0

        pmc = self.local_pmc
        # Multiplicity-one condition
        total_mult = [m1+m2 for m1, m2 in zip(gen1.multiplicity,
                                              gen2.multiplicity)]
        if not all([x <= 1 for x in total_mult]):
            return E0

        new_strands = []

        # Keep track of which strands at right are not yet used.
        strands_right = list(gen2.strands)
        for sd in gen1.strands:
            mid_idem = pmc.pairid[sd[1]]
            if mid_idem == -1:
                # Strands going to the boundary go to the product
                new_strands.append((sd[0], sd[1]))
                continue
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
        # Since we are in the multiplicity-one case, no need to worry about
        # double-crossing. Can return now.
        return LocalStrandDiagram(self, gen1.left_idem, new_strands).elt()

class PMCSplitting:
    """Contains information about a splitting of a full PMC into two local PMCs.

    """

    def __init__(self, pmc, intervals):
        """Given a full PMC and a list of intervals (specified by pairs),
        construct the splitting where local_pmc is the restriction of pmc to the
        intervals, and outer_pmc is the restriction of pmc to the complement of
        the intervals.

        - intervals: must be ordered, disjoint intervals. Each interval is
        specified in the format (start, end), which represents the interval
        between these two points. Intervals with start > end are ignored.

        """
        self.pmc = pmc
        self.intervals = tuple(intervals)
        outer_intervals = PMCSplitting.complementIntervals(self.pmc, intervals)
        self.local_pmc, self.local_mapping = \
            PMCSplitting.restrictPMC(self.pmc, intervals)
        self.outer_pmc, self.outer_mapping = \
            PMCSplitting.restrictPMC(self.pmc, outer_intervals)

    def __str__(self):
        return "Splitting of %s with local intervals %s.\n" % (self.pmc,
                                                               self.intervals)

    def __eq__(self, other):
        return self.pmc == other.pmc and self.intervals == other.intervals

    def __ne__(self, other):
        return not (self == other)

    @memorizeHash
    def __hash__(self):
        return hash(("PMCSplitting", self.pmc, self.intervals))

    def restrictStrandDiagramLocal(self, sd):
        """Returns the local strand diagram that is the restriction of sd to
        local_pmc.

        """
        # Memorize within sd.
        if not hasattr(sd, "restrictLocal"):
            sd.restrictLocal = dict()

        if self not in sd.restrictLocal:
            rv = PMCSplitting.restrictStrandDiagram(
                self.pmc, sd, self.local_pmc, self.local_mapping)
            sd.restrictLocal[self] = rv
            return rv
        else:
            return sd.restrictLocal[self]

    def restrictStrandDiagramOuter(self, sd):
        """Returns the local strand diagram that is the restriction of sd to
        outer_pmc.

        """
        # Memorize within sd.
        if not hasattr(sd, "restrictOuter"):
            sd.restrictOuter = dict()

        if self not in sd.restrictOuter:
            rv = PMCSplitting.restrictStrandDiagram(
                self.pmc, sd, self.outer_pmc, self.outer_mapping)
            sd.restrictOuter[self] = rv
            return rv
        else:
            return sd.restrictOuter[self]

    def restrictIdempotentLocal(self, idem):
        """Returns the local idempotent that is the restriction of idem to
        local_pmc.

        """
        return PMCSplitting.restrictIdempotent(
            self.pmc, idem, self.local_pmc, self.local_mapping)

    def restrictIdempotentOuter(self, idem):
        """Returns the local idempotent that is the restriction of idem to
        outer_pmc.

        """
        return PMCSplitting.restrictIdempotent(
            self.pmc, idem, self.outer_pmc, self.outer_mapping)

    def joinStrandDiagram(self, sd1, sd2):
        """Joins two strand diagrams. sd1 is in local_pmc and sd2 is in
        outer_pmc.
        Returns None if these two local strand diagrams cannot be joined.
        Otherwise returns the strand diagram in a full PMC.

        """
        assert sd1.parent.local_pmc == self.local_pmc
        assert sd2.parent.local_pmc == self.outer_pmc

        # Check idempotent is OK. Create joined left_idem.
        left_idem = []
        mappings = {1 : self.local_mapping, 2 : self.outer_mapping}
        local_pmcs = {1 : self.local_pmc, 2 : self.outer_pmc}
        local_left_idem = {1 : sd1.left_idem, 2 : sd2.left_idem}
        for pairid in range(self.pmc.num_pair):
            p, q = self.pmc.pairs[pairid]
            for i in (1, 2):
                if p in mappings[i]:
                    has_idem_p = (local_pmcs[i].pairid[mappings[i][p]]
                                  in local_left_idem[i])
                if q in mappings[i]:
                    has_idem_q = (local_pmcs[i].pairid[mappings[i][q]]
                                  in local_left_idem[i])
            if has_idem_p or has_idem_q:
                left_idem.append(pairid)

        # Construct inverse mapping from points in local_pmc to points in pmc.
        inv_mappings = {1 : {}, 2 : {}}
        for pt, local_pt in mappings[1].items():
            inv_mappings[1][local_pt] = pt
        for pt, local_pt in mappings[2].items():
            inv_mappings[2][local_pt] = pt

        # Create joined strands.
        # First create list of local strands in sorted order in original pmc.
        local_strands = {1 : sd1.strands, 2 : sd2.strands}
        all_local_strands = []
        for i in (1, 2):
            for start, end in local_strands[i]:
                if start in local_pmcs[i].endpoints:
                    assert start+1 not in local_pmcs[i].endpoints
                    start_pos = inv_mappings[i][start + 1]
                else:
                    start_pos = inv_mappings[i][start]
                all_local_strands.append((start_pos, (start, end), i))
        all_local_strands = sorted(all_local_strands)

        # First check that every loose end is matched. Otherwise return.
        for strands_id in range(len(all_local_strands)):
            start_pos, (start, end), pmc_id = all_local_strands[strands_id]
            if start in local_pmcs[pmc_id].endpoints:
                if strands_id == 0:
                    return None
                prev_pos, (prev_start, prev_end), prev_pmc_id = \
                    all_local_strands[strands_id - 1]
                if prev_end not in local_pmcs[prev_pmc_id].endpoints:
                    return None
                # Check boundaries match in one of the two directions
                if prev_pmc_id == pmc_id:
                    return None
                if inv_mappings[pmc_id][start+1] - 1 != \
                   inv_mappings[prev_pmc_id][prev_end-1]:
                    return None
            if end in local_pmcs[pmc_id].endpoints:
                if strands_id == len(all_local_strands) - 1:
                    return None
                next_pos, (next_start, next_end), next_pmc_id = \
                    all_local_strands[strands_id + 1]
                if next_start not in local_pmcs[next_pmc_id].endpoints:
                    return None

        # Having made sure that every loose end is closed, we can simply take
        # the sequence of non-endpoints.
        all_strand_boundaries = []
        for start_pos, (start, end), pmc_id in all_local_strands:
            for pt in (start, end):
                if pt not in local_pmcs[pmc_id].endpoints:
                    all_strand_boundaries.append(inv_mappings[pmc_id][pt])
        # Now, simply take pairs:
        start_pts = all_strand_boundaries[::2]
        end_pts = all_strand_boundaries[1::2]
        for pt in start_pts:
            if self.pmc.otherp[pt] in start_pts:
                return None
        for pt in end_pts:
            if self.pmc.otherp[pt] in end_pts:
                return None
        strands = zip(start_pts, end_pts)
        if not Strands(self.pmc, strands).leftCompatible(left_idem):
            return None
        return self.pmc.getAlgebra(idem_size = len(left_idem),
                                   mult_one = True).getStrandDiagram(
                                       tuple(left_idem), tuple(strands))

    def joinIdempotent(self, idem1, idem2):
        """Join two local idempotents. """
        return self.joinStrandDiagram(
            idem1.toAlgElt(), idem2.toAlgElt()).left_idem

    @staticmethod
    def complementIntervals(pmc, intervals):
        """Given a full PMC and a list of intervals, return the list of
        intervals that forms the complement.

        """
        result = []
        if len(intervals) == 0:
            result.append((0, pmc.n - 1))
            return result

        if intervals[0][0] != 0:
            result.append((0, intervals[0][0] - 1))
        if intervals[-1][1] != pmc.n - 1:
            result.append((intervals[-1][1] + 1, pmc.n - 1))
        for i in range(len(intervals) - 1):
            # Consider this case later
            assert intervals[i][1] != intervals[i+1][0], \
                "Intervals with no points in the middle is not implemented."
            result.append((intervals[i][1] + 1, intervals[i+1][0] - 1))
        return sorted(result)

    @staticmethod
    def restrictPMC(pmc, intervals):
        """Given a full PMC and a list of intervals, return a pair
        (local_pmc, mapping) where
        - local_pmc is the local PMC that is the restriction of pmc to the
        intervals.
        - mapping is a dictionary mapping from points in pmc to points in
        local_pmc.

        Example:
        PMCSplitting.restrictPMC(PMC([(0, 2),(1, 3)]), [(0, 2)])
          => (LocalPMC(4, [(0, 2), (1,)], [3]), {0:0, 1:1, 2:2})

        """
        # For each interval, add the appropriate endpoints. Fill in endpoints
        # and mapping first.
        num_local_points = 0
        endpoints = []
        mapping = {}
        for start, end in intervals:
            if start > end:
                continue
            length = end - start + 1
            if start != 0:
                endpoints.append(num_local_points)
                num_local_points += 1
            for pt in range(start, end+1):
                mapping[pt] = num_local_points
                num_local_points += 1
            if end != pmc.n - 1:
                endpoints.append(num_local_points)
                num_local_points += 1
        # Now compute the matching.
        matching = []
        for p, q in pmc.pairs:
            if p in mapping and q in mapping:
                matching.append((mapping[p], mapping[q]))
            elif p in mapping and q not in mapping:
                matching.append((mapping[p],))
            elif p not in mapping and q in mapping:
                matching.append((mapping[q],))
        return (LocalPMC(num_local_points, matching, endpoints), mapping)

    @staticmethod
    def restrictStrandDiagram(pmc, sd, local_pmc, mapping):
        """Restrict the given strand diagram to the local_pmc, using mapping as
        the dictionary from points in pmc to points in local_pmc.

        """
        assert sd.parent.pmc == pmc
        # First construct the left idempotent.
        local_left_idem = []
        for (start, end) in sd.strands:
            if start in mapping:
                local_left_idem.append(local_pmc.pairid[mapping[start]])
        for pairid in sd.double_hor:
            p, q = pmc.pairs[pairid]
            if p in mapping:
                local_left_idem.append(local_pmc.pairid[mapping[p]])
            elif q in mapping:
                local_left_idem.append(local_pmc.pairid[mapping[q]])

        # Next, construct strands. For each strand in sd, construct zero or more
        # child strands.
        local_strands = []
        for start, end in sd.strands:
            # Whether to extend the previous item on local_strands. Otherwise,
            # will add new local_strand when needed.
            extend_prev = False
            for pt in range(start, end):
                if pt in mapping:
                    local_pt = mapping[pt]
                    # The interval (pt, pt+1) corresponds to
                    # (local_pt, local_pt+1) in the local PMC. Note local_pt+1
                    # may be an endpoint.
                    if extend_prev:
                        prev_start, prev_end = local_strands[-1]
                        assert prev_end == local_pt
                        local_strands[-1] = (prev_start, local_pt+1)
                    else:
                        if pt != start:
                            assert local_pt - 1 in local_pmc.endpoints
                            local_strands.append((local_pt-1, local_pt+1))
                        else:
                            local_strands.append((local_pt, local_pt+1))
                        extend_prev = True
                    # Turn off extend_prev when endpoint in local_pmc is
                    # reached.
                    if local_pt + 1 in local_pmc.endpoints:
                        extend_prev = False
            # Special case at the end.
            if end - 1 not in mapping and end in mapping:
                local_end = mapping[end]
                assert local_end - 1 in local_pmc.endpoints
                local_strands.append((local_end - 1, local_end))

        return LocalStrandDiagram(local_pmc.getAlgebra(),
                                  local_left_idem, local_strands)

    @staticmethod
    def restrictIdempotent(pmc, idem, local_pmc, mapping):
        """Restrict the given idempotent to the local_pmc, using mapping as the
        dictionary from points in pmc to points in local_pmc.

        """
        local_idem = []
        for pairid in idem:
            p, q = pmc.pairs[pairid]
            if p in mapping:
                local_idem.append(local_pmc.pairid[mapping[p]])
            elif q in mapping:
                local_idem.append(local_pmc.pairid[mapping[q]])
        return LocalIdempotent(local_pmc, local_idem)
