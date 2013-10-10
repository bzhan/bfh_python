"""This module offers minimal support for PMC with boundaries and unmatched
points. A normal PMC (defined in pmc.py) can be split into one or several
local PMC's like this.

"""

from pmc import StrandDiagram

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

def restrictPMC(pmc, intervals):
    """Given a normal PMC, and a list of intervals (specified by pairs), return
    a pair (local_pmc, mapping) where
    - local_pmc is the local PMC that is the restriction of pmc to the
    intervals.
    - mapping is a dictionary mapping from points in pmc to points in
    local_pmc.

    Input parameters:
    - pmc: must be an object of PMC.
    - intervals: must be ordered, disjoint intervals. Each interval is
    specified in the format (start, end), which represents the interval between
    these two points.

    Example:
    restrictPMC(PMC([(0, 2),(1, 3)]), [(0, 2)])
      => (LocalPMC(4, [(0, 2), (1,)], [3]), {0:0, 1:1, 2:2})

    """
    # For each interval, add the appropriate endpoints. Fill in endpoints and
    # mapping first.
    num_local_points = 0
    endpoints = []
    mapping = {}
    for start, end in intervals:
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

class LocalStrands:
    """Represents a fixed list of strands in an local PMC."""
    def __init__(self, local_pmc, data):
        self.local_pmc = local_pmc
        self.data = tuple(sorted(data))
        # Compute multiplicity. multiplicity[i] represents the multiplicity
        # on the interval (i, i+1). Intervals that are gaps between two
        # boundary points should never be occupied, but we do not check for it
        # here.
        self.multiplicity = [0] * (self.local_pmc.n - 1)
        for st in data:
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
        for st in self.data:
            if st[0] not in pmc.endpoints:
                if idem_count[pmc.pairid[st[0]]] == 0: return None
                idem_count[pmc.pairid[st[0]]] -= 1
        for st in self.data:
            if st[1] not in pmc.endpoints:
                if idem_count[pmc.pairid[st[1]]] == 1: return None
                idem_count[pmc.pairid[st[1]]] += 1
        right_idem = [i for i in range(pmc.num_pair) if idem_count[i] == 1]
        return right_idem

class LocalStrandDiagram:
    """A strand diagram in an local PMC. Mainly has to handle multiplication
    two such strand diagrams.

    Multiplicity-one condition is hard-coded in. Dealing with multiplicity
    greater than one when there are strands going to the boundary can be more
    subtle.

    """
    def __init__(self, local_pmc, left_idem, strands):
        """Specifies the pmc, left idempotent, and strands.

        Input parameters:
        - local_pmc: must be an object of LocalPMC.
        - left_idem: tuple containing IDs of occupied pairs.
        - strands: tuple of pairs specifying strands.

        """
        self.local_pmc = local_pmc
        self.left_idem = left_idem
        if not isinstance(strands, LocalStrands):
            strands = LocalStrands(self.local_pmc, strands)
        self.strands = strands

        # Get right_idem and multiplicity from strands
        self.right_idem = self.strands.propagateRight(self.left_idem)
        self.multiplicity = self.strands.multiplicity

    def __eq__(self, other):
        return self.left_idem == other.left_idem and \
            self.strands.data == other.strands.data

    def __ne__(self, other):
        return not (self == other)

    def __hash__(self):
        return hash((self.local_pmc, tuple(self.left_idem),
                     tuple(self.strands.data)))

    def multiply(self, sd2):
        """Multiply two local strand diagrams. Returns None if the product is
        zero. Otherwise return the product as an local strand diagram.

        """
        if not isinstance(sd2, LocalStrandDiagram):
            return NotImplemented
        assert self.local_pmc == sd2.local_pmc, \
            "Local PMC must be the same for multiplication."
        # Idempotent matches
        if self.right_idem != sd2.left_idem:
            return None
        # Multiplicity-one condition
        total_mult = [m1+m2 for m1, m2 in zip(self.multiplicity,
                                              sd2.multiplicity)]
        if not all([x <= 1 for x in total_mult]):
            return None

        pmc = self.local_pmc
        new_strands = []

        # Keep track of which strands at right are not yet used.
        strands_right = list(sd2.strands.data)
        for sd in self.strands.data:
            mid_idem = pmc.pairid[sd[1]]
            possible_match = [sd2 for sd2 in strands_right
                              if pmc.pairid[sd2[0]] == mid_idem]
            if len(possible_match) == 0:
                new_strands.append(sd)
            else: # len(possible_match) == 1
                sd2 = possible_match[0]
                if sd2[0] != sd[1]:
                    return None
                else:
                    new_strands.append((sd[0], sd2[1]))
                    strands_right.remove(sd2)

        new_strands.extend(strands_right)
        # Since we are in the multiplicity-one case, no need to worry about
        # double-crossing. Can return now.
        return LocalStrandDiagram(self.local_pmc,
                                     self.left_idem, new_strands)

    def join(self, sd2, pmc, mapping1, mapping2):
        """Join with another local strand diagram. Returns None if these two
        local strand diagrams cannot be joined. Otherwise returns the (normal)
        strand diagram in a normal PMC.

        Input parameters:
        - pmc: PMC for the joined strand diagram.
        - mapping1 and mapping2 are dictionaries mapping points in pmc to
        points in self.local_pmc and sd2.local_pmc.

        """
        # First make sure that the two local PMC's can be joined into pmc.
        # Can add more checks...
        for pt in range(pmc.n):
            assert (pt in mapping1 and pt not in mapping2) \
                or (pt not in mapping1 and pt in mapping2)

        # Check idempotent is OK. Create joined left_idem.
        left_idem = []
        mapping = {1 : mapping1, 2 : mapping2}
        local_pmc = {1 : self.local_pmc, 2 : sd2.local_pmc}
        local_left_idem = {1 : self.left_idem, 2 : sd2.left_idem}
        for pairid in range(pmc.num_pair):
            p, q = pmc.pairs[pairid]
            for i in (1, 2):
                if p in mapping[i]:
                    has_idem_p = (local_pmc[i].pairid[mapping[i][p]]
                                  in local_left_idem[i])
                if q in mapping[i]:
                    has_idem_q = (local_pmc[i].pairid[mapping[i][q]]
                                  in local_left_idem[i])
            if has_idem_p != has_idem_q:
                return None
            if has_idem_p:
                left_idem.append(pairid)

        # Construct inverse mapping from points in local_pmc to points in pmc.
        inv_mapping = {1 : {}, 2 : {}}
        for pt, local_pt in mapping[1].items():
            inv_mapping[1][local_pt] = pt
        for pt, local_pt in mapping[2].items():
            inv_mapping[2][local_pt] = pt

        # Create joined strands.
        # First create list of local strands in sorted order in original pmc.
        local_strands = {1 : self.strands.data, 2 : sd2.strands.data}
        all_local_strands = []
        for i in (1, 2):
            for start, end in local_strands[i]:
                if start in local_pmc[i].endpoints:
                    assert start+1 not in local_pmc[i].endpoints
                    start_pos = inv_mapping[i][start + 1]
                else:
                    start_pos = inv_mapping[i][start]
                all_local_strands.append((start_pos, (start, end), i))
        all_local_strands = sorted(all_local_strands)

        # First check that every loose end is matched. Otherwise return.
        for strands_id in range(len(all_local_strands)):
            start_pos, (start, end), pmc_id = all_local_strands[strands_id]
            if start in local_pmc[pmc_id].endpoints:
                if strands_id == 0:
                    return None
                prev_pos, (prev_start, prev_end), prev_pmc_id = \
                    all_local_strands[strands_id - 1]
                if prev_end not in local_pmc[prev_pmc_id].endpoints:
                    return None
                # Check boundaries match in one of the two directions
                if prev_pmc_id == pmc_id:
                    return None
                if inv_mapping[pmc_id][start+1] - 1 != \
                   inv_mapping[prev_pmc_id][prev_end-1]:
                    return None
            if end in local_pmc[pmc_id].endpoints:
                if strands_id == len(all_local_strands) - 1:
                    return None
                next_pos, (next_start, next_end), next_pmc_id = \
                    all_local_strands[strands_id + 1]
                if next_start not in local_pmc[next_pmc_id].endpoints:
                    return None

        # Having made sure that every loose end is closed, we can simply take
        # the sequence of non-endpoints.
        all_strand_boundaries = []
        for start_pos, (start, end), pmc_id in all_local_strands:
            for pt in (start, end):
                if pt not in local_pmc[pmc_id].endpoints:
                    all_strand_boundaries.append(inv_mapping[pmc_id][pt])
        # Now, simply take pairs:
        start_pts = all_strand_boundaries[::2]
        end_pts = all_strand_boundaries[1::2]
        for pt in start_pts:
            if pmc.otherp[pt] in start_pts:
                return None
        for pt in end_pts:
            if pmc.otherp[pt] in end_pts:
                return None
        strands = zip(start_pts, end_pts)
        return StrandDiagram(pmc.getAlgebra(idem_size = len(left_idem),
                                            mult_one = True),
                             left_idem, strands)

def restrictStrandDiagram(pmc, sd, local_pmc, mapping):
    """Restrict the given strand diagram to the local_pmc, using mapping as the
    dictionary from points in pmc to points in local_pmc.

    """
    # First construct the left idempotent.
    local_left_idem = []
    for pairid in sd.left_idem:
        p, q = pmc.pairs[pairid]
        if p in mapping:
            local_left_idem.append(local_pmc.pairid[mapping[p]])
        elif q in mapping:
            local_left_idem.append(local_pmc.pairid[mapping[q]])

    # Next, construct strands. For each strand in sd, construct zero or more
    # child strands.
    local_strands = []
    for start, end in sd.strands:
        # Whether to extend the previous item on local_strands. Otherwise, will
        # add new local_strand when needed.
        extend_prev = False
        for pt in range(start, end):
            if pt in mapping:
                local_pt = mapping[pt]
                # The interval (pt, pt+1) corresponds to (local_pt, local_pt+1)
                # in the local PMC. Note local_pt+1 may be an endpoint.
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
                # Turn off extend_prev when endpoint in local_pmc is reached.
                if local_pt + 1 in local_pmc.endpoints:
                    extend_prev = False
        # Special case at the end.
        if end - 1 not in mapping and end in mapping:
            local_end = mapping[end]
            assert local_end - 1 in local_pmc.endpoints
            local_strands.append((local_end - 1, local_end))

    return LocalStrandDiagram(local_pmc, local_left_idem, local_strands)
