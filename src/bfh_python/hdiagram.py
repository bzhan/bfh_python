"""Code for handling Heegaard diagrams."""

from fractions import Fraction
from .grading import BigGradingGroup, SimpleDbGradingSet, \
    SimpleDbGradingSetElement, SimpleGradingSet, SimpleGradingSetElement, \
    SmallGradingGroup
from .grading import DEFAULT_REFINEMENT
from .linalg import RowSystem
from .pmc import Idempotent, PMC
from .utility import SummableDict
from .utility import tolist
from .utility import ACTION_LEFT, ACTION_RIGHT, BIG_GRADING, DEFAULT_GRADING, \
    NEG, POS

"""Constants for each type of segment."""
ALPHA, BETA, BORDER = list(range(3))

def opp(orientation):
    """Return NEG if argument is POS, and POS otherwise."""
    if orientation == POS: return NEG
    else: return POS

class _Point(object):
    """Represents a point."""
    def __init__(self, name):
        self.name = name

    def __str__(self):
        return str(self.name)

    def __repr__(self):
        return "pt(%s)" % str(self.name)

class _Segment(object):
    """Represents a segment."""
    def __init__(self, name, start, end):
        """Specify the name, start point, and end point of this segment."""
        self.name = name
        self.start = start
        self.end = end

    def __str__(self):
        return "%s->%s" % (str(self.start), str(self.end))

    def __repr__(self):
        return "seg(%s: %s->%s)" % \
            (str(self.name), str(self.start), str(self.end))

    def oseg(self):
        """Get the oriented segment with POS direction."""
        return _OrientedSegment(self, POS)

    def orseg(self):
        """Get the oriented segment with NEG direction."""
        return _OrientedSegment(self, NEG)

    def bdZeroChain(self):
        """Returns the zero-chain (end - start)."""
        return _ZeroChain({self.end : 1}) + _ZeroChain({self.start : -1})

class _OrientedSegment(object):
    """Represents an oriented segment as a segment together with an
    orientation.

    """
    def __init__(self, seg, orientation):
        """Specify the underlying segment as well as orientation (either POS or
        NEG

        """
        self.seg = seg
        assert orientation == POS or orientation == NEG, \
            "Orientation cannot be %d" % orientation
        self.orientation = orientation
        if orientation == POS:
            self.start = seg.start
            self.end = seg.end
            self.longname = str(seg.name)
        else:
            self.start = seg.end
            self.end = seg.start
            self.longname = str(seg.name) + "_r"

    def __str__(self):
        return "%s->%s" % (str(self.start), str(self.end))

    def __repr__(self):
        return "oseg(%s: %s->%s)" % \
            (str(self.longname), str(self.start), str(self.end))

    def __eq__(self, other):
        return self.seg == other.seg and self.orientation == other.orientation

    def __ne__(self, other):
        return not (self == other)

    def __hash__(self):
        return hash((self.seg, self.orientation))

    def toOneChain(self):
        """Return the one-chain of this oriented segment."""
        return _OneChain({self.seg : self.orientation})

    def opp(self):
        """Return the opposite oriented segment."""
        return _OrientedSegment(self.seg, opp(self.orientation))

class _Path(list):
    """A path represented as a list of OrientedSegment."""
    def __init__(self, data = None, name = "", iscycle = False):
        """Specify the data as a list of OrientedSegment. Optionally, provide
        a name for the path, and/or whether this path should be considered as a
        cycle.

        """
        if data is None:
            data = []
        list.__init__(self, data)
        self.name = name
        self.iscycle = iscycle
        if not self._checkValidity():
            raise ValueError("Invalid path")

    def _checkValidity(self):
        """Performs the following checks:
        #. Whether the end point of each segment agrees with the start point of
        the next one.
        #. If ``iscycle`` is set to ``True``, then check whether the end point
        of the last segment agrees with the start point of the first.

        """
        for i in range(len(self)-1):
            if self[i].end != self[i+1].start:
                return False
        if self.iscycle and self[-1].end != self[0].start:
            return False
        return True

    def __str__(self):
        result = str(self[0].start)
        result += "".join(["->%s" % oseg.end for oseg in self])
        return result

    def __repr__(self):
        if self.name == "":
            return "path(%s)" % str(self)
        else:
            return "path(%s: %s)" % (str(self.name), str(self))

    def toOneChain(self):
        """Return the one-chain of this path."""
        return _OneChain().accumulate([oseg.toOneChain() for oseg in self])

    def opp(self, name = ""):
        """Return the opposite path."""
        return _Path(reversed([oseg.opp() for oseg in self]),
                     name, self.iscycle)

    def subPath(self, start, end):
        """Returns a part of this path between start and end. If start <= end,
        then the returned path is self[start:end]. If start > end, then the
        returned path is the reverse of self[end:start]. One can imagine
        labeling the start of segment i in the path as point i, then this
        function returns the part of the path between points start and end,
        with the appropriate orientation.

        This function will correctly handle -1 (same as length of path).

        """
        if start == -1: start = len(self)
        if end == -1: end = len(self)
        if start <= end:
            return _Path(self[start:end])
        else:
            return _Path(self[end:start]).opp()

class _Cell(object):
    """Represents an oriented cell."""
    def __init__(self, name, boundary):
        """Specify the name of the cell, and the boundary as a Path (or a list
        of Path) with ``iscycle`` set to ``True``.

        """
        self.name = name
        if isinstance(boundary, _Path):
            self.boundary = [boundary]
        else:
            self.boundary = boundary
        for b in self.boundary:
            assert b.iscycle, "Boundary of cell must be a cycle"

    def __str__(self):
        return str(self.boundary[0])

    def __repr__(self):
        if self.name == "":
            return "cell(%s)" % str(self.boundary[0])
        else:
            return "cell(%s: %s)" % (str(self.name), str(self.boundary[0]))

    def toDomain(self):
        """Return the domain of this cell."""
        return _Domain({self : 1})

    def bdOneChain(self):
        """Return the one chain represented by its boundary."""
        return _OneChain().accumulate(
            [bd.toOneChain() for bd in self.boundary])

class _Domain(SummableDict):
    """A signed sum of oriented cells, represented as a dictionary from Cell
    to integers.

    """
    def diff(self):
        """Return the boundary of this domain as a one-chain."""
        return _OneChain().accumulate([coeff * cell.bdOneChain()
                                       for cell, coeff in list(self.items())])

class _OneChain(SummableDict):
    """A signed sum of oriented segments, represented as a dictionary from
    Segment to integers.

    """
    def diff(self):
        """Return the boundary of this one-chain as a zero-chain."""
        return _ZeroChain().accumulate([coeff * seg.bdZeroChain()
                                        for seg, coeff in list(self.items())])

class _ZeroChain(SummableDict):
    """A signed sum of points, represented as a dictionary from Point to
    integers.

    """
    pass

class HFGenerator(object):
    """Represents a generator of the Heegaard Floer chain complex for a certain
    Heegaard diagram.

    """
    def __init__(self, parent, points):
        """Specifies the parent Heegaard diagram, and a list of points in that
        diagram.

        """
        self.parent = parent
        self.points = points

    def getIdem(self):
        """Returns the idempotent of this generator. For corresponding function
        in HeegaardDiagram class for details.

        """
        return self.parent.getGeneratorIdem(self)

    def getIdemSize(self):
        """Returns the size of the type A idempotent on the first PMC. In the
        one border case this is always the genus. In the two border cases
        generators are grouped by this size.

        """
        idems = self.getIdem()
        if len(idems) == 0:
            return 0  # closed diagram.
        return len(idems[0])

    def getDIdem(self):
        """Returns the complement of the idempotent of this generator, in the
        opposite PMC. (used in type D and type DD structures)

        """
        idem = self.getIdem()
        if isinstance(idem, Idempotent):
            return idem.opp().comp()
        else:
            return [i.opp().comp() for i in idem]

    def toZeroChain(self):
        """Returns the zero-chain of this generator."""
        return _ZeroChain().accumulate([_ZeroChain({pt : 1})
                                        for pt in self.points])

class HeegaardDiagram(object):
    """Represents a possibly bordered Heegaard Diagram."""
    def __init__(self, name, points, segments, cells, alpha, beta, border,
                 basept):
        """Specify the information needed to build a bordered Heegaard diagram.
        Elements of border, and elements of alpha and beta representing cycles
        must be Path with ``iscycle`` set to ``True``. Elements of alpha and
        beta representing arcs must be Path with ``iscycle`` set to ``False``.

        """
        self.name = name
        self.points = points
        self.segments = segments
        self.cells = cells
        self.alpha = tolist(alpha)
        self.beta = tolist(beta)
        self.border = tolist(border)
        self.basept = tolist(basept)

        self.num_point = len(self.points)
        self.num_segment = len(self.segments)
        self.num_cell = len(self.cells)

        # Separate both alpha and beta lists into cycles and arcs
        self.alpha_cycles, self.alpha_arcs = [], []
        self.beta_cycles, self.beta_arcs = [], []
        for path in self.alpha:
            if path.iscycle: self.alpha_cycles.append(path)
            else: self.alpha_arcs.append(path)
        for path in self.beta:
            if path.iscycle: self.beta_cycles.append(path)
            else: self.beta_arcs.append(path)
        assert len(self.alpha_arcs) % 2 == 0 and len(self.beta_arcs) % 2 == 0

        # Compute genus = number of cycles + (number of arcs / 2)
        self.genus = len(self.alpha_cycles) + (len(self.alpha_arcs) // 2)
        assert self.genus == len(self.beta_cycles) + (len(self.beta_arcs) // 2)
        # Compute pmc_genus = length of border / 4
        assert all([len(path) % 4 == 0 for path in self.border])
        self.pmc_genus = [len(path) // 4 for path in self.border]

        # For each point, compute its position in the alpha, beta, or border
        # paths. Represent them as tuples in the alpha_pos, beta_pos, and
        # border_pos fields of each point. The field is None if not applicable
        # for that point. Each tuple is (path index, point index within path).
        # Both indices start at zero.
        # For each segment, record its type (one of ALPHA, BETA, and BORDER),
        # and info: the tuple (path index, index within path, and orientation).
        for pt in points:
            pt.alpha_info = pt.beta_info = pt.border_info = None
        for alpha_id in range(len(self.alpha)):
            alpha_path = self.alpha[alpha_id]
            for oseg_id in range(len(alpha_path)):
                oseg = alpha_path[oseg_id]
                oseg.seg.type = ALPHA
                oseg.seg.info = (alpha_id, oseg_id, oseg.orientation)
                oseg.start.alpha_info = (alpha_id, oseg_id)
            if not alpha_path.iscycle:
                alpha_path[-1].end.alpha_info = (alpha_id, len(alpha_path))
        for beta_id in range(len(self.beta)):
            beta_path = self.beta[beta_id]
            for oseg_id in range(len(beta_path)):
                oseg = beta_path[oseg_id]
                oseg.seg.type = BETA
                oseg.seg.info = (beta_id, oseg_id, oseg.orientation)
                oseg.start.beta_info = (beta_id, oseg_id)
            if not beta_path.iscycle:
                beta_path[-1].end.beta_info = (beta_id, len(beta_path))
        for border_id in range(len(self.border)):
            border_path = self.border[border_id]
            for oseg_id in range(len(border_path)):
                oseg = border_path[oseg_id]
                oseg.seg.type = BORDER
                oseg.seg.info = (border_id, oseg_id, oseg.orientation)
                oseg.start.border_info = (border_id, oseg_id)
            assert border_path.iscycle

        # Compute the PMC on each border. In addition, associate to each arc
        # the tuple (border_id, idem_pair) in attribute idem_info
        self.num_pmc = len(self.border)
        pmc_matchings = [[] for i in range(self.num_pmc)]
        for path in self.alpha_arcs + self.beta_arcs:
            start_pt = path[0].start
            end_pt = path[-1].end
            start_border_id, start_border_pos = start_pt.border_info
            end_border_id, end_border_pos = end_pt.border_info
            assert start_border_id == end_border_id
            idem_pair = (start_border_pos, end_border_pos)
            pmc_matchings[start_border_id].append(idem_pair)
            path.idem_info = (start_border_id, idem_pair)

        assert all([self.pmc_genus[i] == len(pmc_matchings[i]) // 2
                    for i in range(len(pmc_matchings))])
        self.pmc_list = [PMC(matching) for matching in pmc_matchings]

        # Compute the list of generators
        self._computeHFGenerators()

        # Construct row system of domains
        self._computeRowSystem()

    def _computeHFGenerators(self):
        """Compute the list of generators. Produces _generator_list attribute
        in self. This is a list of generators if there are less than two
        borders, and a list of lists of generators, grouped by size of
        idempotent on the first border.

        """
        self.two_border_case = (len(self.border) == 2)

        if not self.two_border_case:
            self._generator_list = []
        else:
            self._generator_list = [[] for i in range(self.pmc_genus[0]*2+1)]

        def helper(pts_chosen, pt_pos, alpha_occupied, beta_occupied):
            """pts_chosen is the list of points in the current partial
            generator. pt_pos is the index of the next point to be considered
            for inclusion. alpha_occupied is a list of booleans on whether each
            alpha cycle / arc is already used. Likewise for beta_occupied.

            """
            if len(pts_chosen) == self.genus:
                # Check each alpha cycle and beta cycle are used
                for i in range(len(self.alpha)):
                    if self.alpha[i].iscycle and not alpha_occupied[i]:
                        return
                for i in range(len(self.beta)):
                    if self.beta[i].iscycle and not beta_occupied[i]:
                        return
                # Can add this generator
                hf_gen = HFGenerator(self, pts_chosen)
                if not self.two_border_case:
                    self._generator_list.append(hf_gen)
                    return
                else:
                    self._generator_list[hf_gen.getIdemSize()].append(hf_gen)
                    return
            # Need to add more points
            if pt_pos == len(self.points):
                return
            helper(pts_chosen, pt_pos+1, alpha_occupied,
                            beta_occupied)
            cur_pt = self.points[pt_pos]
            if cur_pt.border_info != None:
                # Cannot add points on the border
                return
            cur_alpha_id = cur_pt.alpha_info[0]
            cur_beta_id = cur_pt.beta_info[0]
            if alpha_occupied[cur_alpha_id] or beta_occupied[cur_beta_id]:
                return
            # Attempt to add this point
            alpha_occupied[cur_alpha_id] = beta_occupied[cur_beta_id] = True
            helper(pts_chosen+[cur_pt], pt_pos+1, alpha_occupied,
                            beta_occupied)
            alpha_occupied[cur_alpha_id] = beta_occupied[cur_beta_id] = False

        # Finally call with initial conditions
        helper([], 0, [False]*len(self.alpha), [False]*len(self.beta))

    def _computeRowSystem(self):
        """Construct the row system of segments associated to the diagram. Each
        row is a vector in R^n, where n is the number of interior segments.
        Segments on the border are ignored. Each row corresponds to either a
        non-basepoint cell or an alpha/beta path. This is used to find both
        periodic domains (linear relations on the rows) and domains connecting
        two generators (represent the boundary as a linear combination of
        rows).

        Besides the row system self.row_sys, this function constructs the
        following:
        *. used_cells: list of non-basepoint cells. The indices of cells in
        this array agree with their indices in the row system.
        *. num_used_cell: number of non-basepoint cells.
        *. interior_segs: non-border segments. The indices of segments in this
        array agree with their (column) indices in the row system.
        *. seg_to_id: translating segments to their index in interior_segs.
        *. num_interior_seg: number of interior segments.

        """
        # Construct list of non-basepoint cells.
        self.used_cells = list(set(self.cells) - set(self.basept))
        self.num_used_cell = len(self.used_cells)

        # Construct list of interior segments, as well as a dictionary
        # translating segments to indices.
        self.interior_segs = [seg for seg in self.segments
                              if seg.type != BORDER]
        self.seg_to_id = dict()
        self.num_interior_seg = len(self.interior_segs)
        for i in range(self.num_interior_seg):
            self.seg_to_id[self.interior_segs[i]] = i

        rows = [cell.bdOneChain() for cell in self.used_cells] + \
            [path.toOneChain() for path in self.alpha + self.beta]
        vec_rows = []
        for row in rows:
            cur_vec = [0] * self.num_interior_seg
            for seg, coeff in list(row.items()):
                if seg.type != BORDER:
                    cur_vec[self.seg_to_id[seg]] += coeff
            vec_rows.append(cur_vec)
        self.row_sys = RowSystem(vec_rows)

    def __str__(self):
        return "HeegaardDiagram(%s)" % str(self.name)

    def __repr__(self):
        output = str(self)
        output += "\n\nPoints: %s" % str(self.points)
        output += "\n\nSegments: %s" % str(self.segments)
        output += "\n\nCells: %s" % str(self.cells)
        output += "\n\nBase point: %s" % repr(self.basept[0])
        output += "\n"
        return output

    def getPMCs(self):
        """Return the PMC's (zero, one, or two) associated to each border."""
        return self.pmc_list

    def restrictOneChain(self, one_chain, seg_type):
        """Restrict the one chain to include segments of only the given segment
        type.

        """
        return _OneChain().accumulate([_OneChain({seg : coeff})
                                       for seg, coeff in list(one_chain.items())
                                       if seg.type == seg_type])

    def restrictZeroChain(self, zero_chain):
        """Restrict the zero chain to include only points in the interior."""
        return _ZeroChain().accumulate([_ZeroChain({pt : coeff})
                                        for pt, coeff in list(zero_chain.items())
                                        if pt.border_info is None])

    def getHFGenerators(self, idem_size = None):
        """Get generators of the Heegaard Floer chain complex associated to
        this Heegaard diagram. Each generator is a tuple of self.genus interior
        points. To be a valid generator, each cycle must contain exactly one
        point, and each arc must contain at most one point. By default, exactly
        half of the arcs on each side contain points. This can be changed by
        specifying idem_size, which is the number of arcs on the first border
        (usually considered the left border) that contain points.

        """
        if not self.two_border_case:
            # idem_size is meaningless then
            assert idem_size is None
            return self._generator_list
        else:
            if idem_size is None:
                # Defaults to half of the arcs occupied
                idem_size = self.pmc_genus[0]
            return self._generator_list[idem_size]

    def getGeneratorByIdem(self, idem, D_idem = False):
        """Returns the first generator with the given idempotent. Returns None
        if none exist. Note this takes time linear in the number of generators.
        (can be improved)

        """
        if D_idem:
            idem_size = self.pmc_list[0].num_pair - len(idem[0])
        else:
            idem_size = len(idem[0])
        for gen in self.getHFGenerators(idem_size):
            if not D_idem and gen.getIdem() == idem:
                return gen
            if D_idem and gen.getDIdem() == idem:
                return gen
        return None

    def getGeneratorIdem(self, generator):
        """Find the idempotent of this generator. This gives the idempotent of
        arcs occupied by this generator (that is, type A idempotents), in the
        same PMC's as that will be returned by getPMCs().

        """
        num_border = len(self.border)
        raw_idems = [[] for i in range(num_border)]
        for pt in generator.points:
            assert pt.alpha_info is not None and pt.beta_info is not None
            alpha_path = self.alpha[pt.alpha_info[0]]
            beta_path = self.beta[pt.beta_info[0]]
            for path in (alpha_path, beta_path):
                if not path.iscycle: # is an arc that goes to the boundary
                    border_id, pair = path.idem_info
                    pairid = self.pmc_list[border_id].pairid[pair[0]]
                    raw_idems[border_id].append(pairid)
        idems = []
        for i in range(num_border):
            idems.append(Idempotent(self.pmc_list[i], raw_idems[i]))
        return idems

    def getPeriodicDomains(self):
        """Get the list of periodic domains."""
        vec_domains = self.row_sys.getZeroComb()
        domains = []
        for vec_domain in vec_domains:
            domains.append(_Domain().accumulate(
                    [self.used_cells[i].toDomain() * vec_domain[i]
                     for i in range(self.num_used_cell)]))
        return domains

    def getConnectingDomain(self, gen1, gen2):
        """Get a domain connecting two generators. Note the returned value is
        not uniquely specified (only unique up to periodic domains).

        The alpha part of the boundary of this domain goes from gen1 to gen2,
        while the beta part goes from gen2 to gen1.

        """
        bdChain = _OneChain()
        alpha_id_to_pos = dict()
        beta_id_to_pos = dict()
        # For each point in gen1: if that point is on an arc, immediately add
        # the appropriate partial arc. Otherwise record the position on the
        # cycle.
        for pt in gen1.points:
            alpha_id, oseg_id = pt.alpha_info
            alpha_path = self.alpha[alpha_id]
            if alpha_path.iscycle:
                alpha_id_to_pos[alpha_id] = oseg_id
            else:
                bdChain += alpha_path.subPath(oseg_id, -1).toOneChain()
            beta_id, oseg_id = pt.beta_info
            beta_path = self.beta[beta_id]
            if beta_path.iscycle:
                beta_id_to_pos[beta_id] = oseg_id
            else:
                bdChain += beta_path.subPath(0, oseg_id).toOneChain()

        # For each point on gen2: if that point is on an arc, add the partial
        # arc like before. Otherwise, use the dictionaries constructed for gen1
        # to find the part of cycle to add
        for pt in gen2.points:
            alpha_id, oseg_id = pt.alpha_info
            alpha_path = self.alpha[alpha_id]
            if alpha_path.iscycle:
                bdChain += alpha_path.subPath(
                    alpha_id_to_pos[alpha_id], oseg_id).toOneChain()
            else:
                bdChain += alpha_path.subPath(0, oseg_id).toOneChain()
            beta_id, oseg_id = pt.beta_info
            beta_path = self.beta[beta_id]
            if beta_path.iscycle:
                bdChain += beta_path.subPath(
                    oseg_id, beta_id_to_pos[beta_id]).toOneChain()
            else:
                bdChain += beta_path.subPath(oseg_id, -1).toOneChain()

        # Now translate bdChain to vectors in the row system.
        bd_vec = [0] * self.num_interior_seg
        for seg, coeff in list(bdChain.items()):
            assert seg.type != BORDER
            bd_vec[self.seg_to_id[seg]] += coeff
        vec_domain = self.row_sys.getComb(bd_vec)
        if vec_domain is None:
            return None
        else:
            return _Domain().accumulate(
                [self.used_cells[i].toDomain() * vec_domain[i]
                 for i in range(self.num_used_cell)])

    def getMaslov(self, domain, x, y):
        """Returns the Maslov grading of a domain connecting generators x and
        y. This is defined as:
           \mu(B) = -e(B) - n_x(B) - n_y(B),
        where e(B) is the Euler measure, and n_x, n_y are average
        multiplicities of points in the generators on the domain.

        """
        # We will use the fact that each corner of an individual cell counts as
        # a right angle (add 1/4 to n_x and n_y).
        maslov = Fraction(0)
        xpts, ypts = x.points, y.points
        for cell, coeff in list(domain.items()):
            cell_maslov = 1 - len(cell.boundary) # correction to e(B)
            for boundary in cell.boundary:
                cell_maslov += (4-len(boundary)) / Fraction(4) # e(B)
                for oseg in boundary:
                    if oseg.end in xpts:
                        cell_maslov += Fraction(1,4) # n_x(B)
                    if oseg.end in ypts:
                        cell_maslov += Fraction(1,4) # n_y(B)
            maslov -= cell_maslov * coeff
        return maslov

    def getBigGrading(self, domain, x, y):
        """Returns the big grading (that is, not refined) of a domain
        connecting generators x and y. This has a Maslov part and a Spin-c
        part that records the part of boundary of domain on the border. This
        should be used only when the diagram has at least one border. The
        returned value is a list of BigGradingElement, one for each border (and
        PMC). The first grading contains the Maslov part, the remaining ones
        have Maslov part zero.

        """
        maslov = self.getMaslov(domain, x, y)
        multiplicity = []
        pmcs = self.getPMCs()
        for pmc in pmcs:
            multiplicity.append([0] * (pmc.n-1))
        bd_border = self.restrictOneChain(domain.diff(), BORDER)
        for seg, coeff in list(bd_border.items()):
            assert seg.type == BORDER
            border_id, oseg_id, orientation = seg.info
            # Cannot be the edge with the basepoint
            assert oseg_id != pmcs[border_id].n - 1
            multiplicity[border_id][oseg_id] += orientation * coeff
        result = []
        for i in range(len(pmcs)):
            result.append(pmcs[i].big_gr(0, multiplicity[i]))
        result[0].maslov = maslov
        return result

    def getBigGradingD(self, domain, x, y):
        """Returns the big grading of domain used in type D structures (with
        Ropp()).

        """
        return [a.Ropp() for a in self.getBigGrading(domain, x, y)]

    def getSmallGradingD(self, domain, x, y):
        """Returns the small (refined) grading of a domain connecting
        generators x and y. Similar to and uses getBigGrading as a first step.

        """
        big_gr = self.getBigGradingD(domain, x, y)
        small_gr = []
        for i in range(len(big_gr)):
            if len(x.getIdem()) == 2:
                # In 2 boundary component case, find size of the current
                # (D-side) idempotent.
                idem_size = len(x.getDIdem()[i])
            else:
                idem_size = None
            cur_big_gr = big_gr[i]
            cur_pmc_opp = self.pmc_list[i].opp()
            refine_data = DEFAULT_REFINEMENT(cur_pmc_opp, idem_size)
            phix, phiy = [refine_data[p.getDIdem()[i]] for p in (x, y)]
            small_gr.append(
                (phiy * cur_big_gr * phix.inverse()).toSmallGrading())
        return small_gr

    def computeDGrading(self, base_gen, base_gr = None):
        """Compute type D grading for each generator, with the grading of
        generator base_gen set to base_gr (defaults to zero grading).

        This function assumes that the diagram has one border, and that all
        generators are linked by domains. The return value is the pair
        consisting of the grading set and the dictionary mapping generators to
        elements of a SimpleGradingSet.

        """
        if DEFAULT_GRADING == BIG_GRADING:
            gr_fun, gr_group_cls = self.getBigGradingD, BigGradingGroup
        else:
            gr_fun, gr_group_cls = self.getSmallGradingD, SmallGradingGroup
        # First construct grading set (mainly construct gradings of periodic
        # domains)
        assert len(self.pmc_list) == 1
        pmc = self.pmc_list[0]
        pmc_opp = pmc.opp()
        if base_gr is None:
            base_gr = gr_group_cls(pmc_opp).zero()
        periodic_domains = self.getPeriodicDomains()
        domains_gr = [gr_fun(domain, base_gen, base_gen)[0]
                      for domain in periodic_domains]
        domains_gr = [base_gr.inverse() * domain_gr * base_gr
                      for domain_gr in domains_gr]
        gr_set = SimpleGradingSet(gr_group_cls(pmc_opp), ACTION_LEFT,
                                  domains_gr)
        # Now compute grading of each generator by finding domains connecting
        # it to the base generator.
        result = dict()
        for gen in self.getHFGenerators():
            conn_domain = self.getConnectingDomain(base_gen, gen)
            domain_gr = gr_fun(conn_domain, base_gen, gen)[0]
            domain_gr = domain_gr * base_gr
            result[gen] = SimpleGradingSetElement(gr_set, domain_gr)
        return (gr_set, result)

    def computeDDGrading(self, base_gen, base_gr = None):
        """Compute type DD grading for each generator, with the grading of
        generator base_gen set to base_gr (defaults to zero grading).

        This function assumes that the diagram has two borders, and that all
        generators are linked by domains. The return value is the pair
        consisting of the grading set and the dictionary mapping generators to
        elements of a SimpleDbGradingSet.

        """
        if DEFAULT_GRADING == BIG_GRADING:
            gr_fun, gr_group_cls = self.getBigGradingD, BigGradingGroup
        else:
            gr_fun, gr_group_cls = self.getSmallGradingD, SmallGradingGroup
        # First construct grading set (mainly construct gradings of periodic
        # domains)
        assert len(self.pmc_list) == 2
        pmc1, pmc2 = self.pmc_list
        pmc1_opp, pmc2_opp = pmc1.opp(), pmc2.opp()
        if base_gr is None:
            base_gr = [gr_group_cls(pmc1_opp).zero(),
                       gr_group_cls(pmc2_opp).zero()]
        periodic_domains = self.getPeriodicDomains()
        domains_gr = [gr_fun(domain, base_gen, base_gen)
                      for domain in periodic_domains]
        domains_gr = [[base_gr[0].inverse() * domain_gr0 * base_gr[0],
                       base_gr[1].inverse() * domain_gr1 * base_gr[1]]
                      for domain_gr0, domain_gr1 in domains_gr]
        gr_set = SimpleDbGradingSet(gr_group_cls(pmc1_opp), ACTION_LEFT,
                                    gr_group_cls(pmc2_opp), ACTION_LEFT,
                                    domains_gr)
        # Now compute grading of each generator by finding domains connecting
        # it to the base generator.
        result = dict()
        for gen in self.getHFGenerators(base_gen.getIdemSize()):
            conn_domain = self.getConnectingDomain(base_gen, gen)
            domain_gr0, domain_gr1 = gr_fun(conn_domain, base_gen, gen)
            domain_gr = [domain_gr0 * base_gr[0], domain_gr1 * base_gr[1]]
            result[gen] = SimpleDbGradingSetElement(gr_set, domain_gr)
        return (gr_set, result)

    def computeDAGrading(self, base_gen, base_gr = None):
        """Compute type DA grading for each generator. The type DA grading is
        related in a straightforward way to the type DD grading. So see the
        function computeDDGrading for details.

        """
        if base_gr is not None:
            return NotImplemented  # code the other case later

        dd_gr_set, dd_result = self.computeDDGrading(base_gen, base_gr)
        lr_domains = [(d1, d2.Ropp()) for d1, d2 in dd_gr_set.periodic_domains]
        gr_set = SimpleDbGradingSet(
            dd_gr_set.gr_group1, ACTION_LEFT,
            dd_gr_set.gr_group2.opp(), ACTION_RIGHT, lr_domains)
        result = dict()
        for x, gr_x in list(dd_result.items()):
            gr_x1, gr_x2 = gr_x.data
            gr_x_lr = SimpleDbGradingSetElement(gr_set, (gr_x1, gr_x2.Ropp()))
            result[x] = gr_x_lr
        return (gr_set, result)

def diagramFromCycleInfo(name, num_interior_point = 0, length_border = [],
                         alpha_cycles = [], alpha_arcs = [],
                         beta_cycles = [], beta_arcs = [],
                         crossing_orientation = [], cell_info = [],
                         basept = []):
    """Building a Heegaard diagram by specifying points on each cycle. The
    meanings of each argument is as follows:
    *. ``num_interior_point`` - number of interior points.
    *. ``length_border`` - list consisting of number of points on each border.
    The points will be numbered in the counter-clockwise direction.
    *. ``alpha_cycles`` - list of list of points on each alpha cycle.
    *. ``alpha_arcs`` - list of list of points on each alpha arc, except for
    the first and last element of each list, which is a pair ``(b, n)`` where
    ``b`` is the index of the border and ``n`` is the index of the point on the
    border. Indices start at 0 and increases in the counterclockwise direction.
    *. ``beta_cycles`` - similar to ``alpha_cycles``.
    *. ``beta_arcs`` - similar to ``alpha_arcs``.
    *. ``crossing_orientation`` - orientation at each interior point. +1 means
    if the alpha cycle/arc goes from left to right, then the beta cycle/arc
    goes from bottom to top. Otherwise -1.
    *. ``cell_info`` - give names for a selection of cells. Each element is
    one of the following:
        **. ((3, 5), "a") - points 3 and 5 are consecutive points on the
        boundary of the cell "a", going counter-clockwise.
        **. ([(3, 5), (4, 6)], "a") - cell "a" has two boundaries.
        **. (("b", 4), "a") - part of cell "a" is the first border, from point
        4 to point 5.
        **. (("b2", 4), "a") - part of cell "a" is the second border. "b1" is
        also recognized.
    *. ``basept`` - Specify the basepoint cell. See ``cell_info`` for giving
    the cell. Can also give name of the cell in place of boundary information.
    For bordered diagrams the obvious basepoint cell is automatically set.

    """

    # Construct the set of points
    interior_points = [_Point("p%d" % i) for i in range(num_interior_point)]
    border_points = [[_Point("b%d,%d" % (i, j))
                      for j in range(length_border[i])]
                     for i in range(len(length_border))]
    points = interior_points + sum(border_points, [])

    def getPt(param):
        # Single integer for interior points, and pair for boundary points
        if isinstance(param, int):
            return interior_points[param]
        else:
            return border_points[param[0]][param[1]]

    # Construct the set of segments and paths
    interior_segs = []
    border_segs = []
    alpha = []
    beta = []
    border = []

    def processPath(ptIDlist, prefix, iscycle):
        # Given a list of points, generate segments connecting consecutive
        # points, as well as a path connecting all.
        ptlist = [getPt(param) for param in ptIDlist]
        newsegs = []
        for i in range(len(ptlist)-1):
            newsegs.append(_Segment("%s,%d" % (prefix, i),
                                    ptlist[i], ptlist[i+1]))
        if iscycle:
            newsegs.append(_Segment("%s,%d" % (prefix, len(ptlist)-1),
                                    ptlist[-1], ptlist[0]))
        newpath = _Path([seg.oseg() for seg in newsegs], prefix, iscycle)
        return newsegs, newpath

    for i in range(len(alpha_cycles)):
        newsegs, newpath = processPath(alpha_cycles[i], "ac%d"%i, True)
        interior_segs.extend(newsegs)
        alpha.append(newpath)

    for i in range(len(alpha_arcs)):
        newsegs, newpath = processPath(alpha_arcs[i], "aa%d"%i, False)
        interior_segs.extend(newsegs)
        alpha.append(newpath)

    for i in range(len(beta_cycles)):
        newsegs, newpath = processPath(beta_cycles[i], "bc%d"%i, True)
        interior_segs.extend(newsegs)
        beta.append(newpath)

    for i in range(len(beta_arcs)):
        newsegs, newpath = processPath(beta_arcs[i], "ba%d"%i, False)
        interior_segs.extend(newsegs)
        beta.append(newpath)

    for i in range(len(length_border)):
        curpath = [(i, j) for j in range(length_border[i])]
        newsegs, newpath = processPath(curpath, "bd%d"%i, True)
        border_segs.append(newsegs)
        border.append(newpath)

    segments = interior_segs + sum(border_segs, [])

    # For each point, compute its three or four neighbors.
    alphain, alphaout, betain, betaout, borderin, borderout = \
        ({},{},{},{},{},{})
    for seg in segments:
        if seg.name[1] == 'd': # border
            borderout[seg.start] = seg
            borderin[seg.end] = seg
        elif seg.name[0] == 'a': # alpha
            alphaout[seg.start] = seg
            alphain[seg.end] = seg
        else: # beta
            betaout[seg.start] = seg
            betain[seg.end] = seg

    # Construct the map from oriented segments to the next one in the same cell
    # in the counterclockwise direction.
    nextseg = {}
    for i in range(num_interior_point):
        pt = interior_points[i]
        ori = crossing_orientation[i]
        if ori == 1:
            nextseg[alphain[pt].oseg()] = betaout[pt].oseg()
            nextseg[betaout[pt].orseg()] = alphaout[pt].oseg()
            nextseg[alphaout[pt].orseg()] = betain[pt].orseg()
            nextseg[betain[pt].oseg()] = alphain[pt].orseg()
        else: # ori == -1
            nextseg[alphain[pt].oseg()] = betain[pt].orseg()
            nextseg[betain[pt].oseg()] = alphaout[pt].oseg()
            nextseg[alphaout[pt].orseg()] = betaout[pt].oseg()
            nextseg[betaout[pt].orseg()] = alphain[pt].orseg()
    for cur_border in border_points:
        for pt in cur_border:
            segoutborder = None # Oriented segment pointing out of border
            if pt in alphaout:
                segoutborder = alphaout[pt].oseg()
            elif pt in alphain:
                segoutborder = alphain[pt].orseg()
            elif pt in betaout:
                segoutborder = betaout[pt].oseg()
            else: # betain.has_key(pt)
                segoutborder = betain[pt].orseg()
            nextseg[borderin[pt].oseg()] = segoutborder
            nextseg[segoutborder.opp()] = borderout[pt].oseg()

    # Construct cells
    cells = []
    # Start with a list of oriented segments that can be used as boundary of
    # cells.
    bdcells = set([seg.oseg() for seg in segments] + \
                      [seg.orseg() for seg in interior_segs])

    while bdcells:
        startbd = bdcells.pop()
        curboundary = [startbd]
        nextbd = nextseg[startbd]
        while nextbd != startbd:
            curboundary.append(nextbd)
            bdcells.remove(nextbd)
            nextbd = nextseg[nextbd]
        cells.append(_Cell("c%d" % (len(cells)+1),
                           _Path(curboundary, iscycle = True)))

    # Construct dictionary from oriented segment to the unique cell for which
    # it is in the counterclockwise boundary.
    segtocell = {}
    for cell in cells:
        for boundary in cell.boundary:
            for oseg in boundary:
                segtocell[oseg] = cell

    # Find the basepoint cell
    allbasept = []
    if len(border) > 0:
        allbasept.append(segtocell[border[0][-1]])

    return HeegaardDiagram(name, points, segments, cells, alpha, beta, border,
                           allbasept)

def getZeroFrameDiagram(genus):
    """Get Heegaard diagram for the 0-framed handlebody of the given genus. The
    left (and only) boundary will have splitPMC(genus) (really the opposite of
    that, but that is the same thing).

    """
    a_arcs = []
    b_cycles = []
    for i in range(genus):
        a_arcs.append([(0,4*i), i, (0,4*i+2)])
        a_arcs.append([(0,4*i+1), (0,4*i+3)])
        b_cycles.append([i])
    return diagramFromCycleInfo("0-framed handlebody of genus %d" % genus,
                                num_interior_point = genus,
                                length_border = [4*genus],
                                alpha_arcs = a_arcs, beta_cycles = b_cycles,
                                crossing_orientation = [1]*genus)

def getZeroFrameAdmDiagram(genus):
    """Get Heegaard diagram for the 0-framed handlebody of a given genus, with
    each beta circle isotopied to cross another alpha arc two times.

    """
    # TODO: figure out grading of the associated type D structure from this
    # diagram? (Can also figure out grading using algebra maps there).
    a_arcs = []
    b_cycles = []
    for i in range(genus):
        a_arcs.append([(0,4*i), 3*i, (0,4*i+2)])
        a_arcs.append([(0,4*i+1), 3*i+2, 3*i+1, (0,4*i+3)])
        b_cycles.append([3*i,3*i+1,3*i+2])
    return diagramFromCycleInfo("0-framed handlebody of genus %d" % genus,
                                num_interior_point = 3*genus,
                                length_border = [4*genus],
                                alpha_arcs = a_arcs, beta_cycles = b_cycles,
                                crossing_orientation = [1,1,-1]*genus)

def getInfFrameDiagram(genus):
    """Get Heegaard diagram for the inf-framed handlebody of the given genus.
    The boundary is as in the 0-framed handlebody.

    """
    a_arcs = []
    b_cycles = []
    for i in range(genus):
        a_arcs.append([(0,4*i), (0,4*i+2)])
        a_arcs.append([(0,4*i+1), i, (0,4*i+3)])
        b_cycles.append([i])
    return diagramFromCycleInfo("Inf-framed handlebody of genus %d" % genus,
                                num_interior_point = genus,
                                length_border = [4*genus],
                                alpha_arcs = a_arcs, beta_cycles = b_cycles,
                                crossing_orientation = [1]*genus)

def getPlatDiagram(genus):
    """Get Heegaard diagram for the plat handlebody of the given genus. The
    left boundary will have linearPMC(genus).

    """
    a_arcs = []
    b_cycles = []
    a_arcs.append([(0,0), (0,2)])
    a_arcs.append([(0,4*genus-3), 2*genus-2, (0,4*genus-1)])
    b_cycles.append([0])
    for i in range(genus-1):
        a_arcs.append([(0,4*i+1), 2*i, 2*i+1, (0,4*i+4)])
        a_arcs.append([(0,4*i+3), (0,4*i+6)])
        b_cycles.append([2*i+1, 2*i+2])
    return diagramFromCycleInfo("Plat handlebody of genus %d" % genus,
                                num_interior_point = 2*genus-1,
                                length_border = [4*genus],
                                alpha_arcs = a_arcs, beta_cycles = b_cycles,
                                crossing_orientation = [1,-1]*(genus-1)+[1])

def getIdentityDiagram(pmc):
    """Get Heegaard diagram for the identity cobordism of the given PMC. The
    left (first) boundary will match the opposite of pmc, and the right
    (second) boundary will match pmc itself.

    """
    n = pmc.n
    a_arcs = []
    b_cycles = []
    for i in range(pmc.num_pair):
        p, q = pmc.pairs[i]
        op, oq = n-1-p, n-1-q
        pt1, pt2 = 2*i, 2*i+1
        a_arcs.append([(1,p), pt1, (1,q)])
        a_arcs.append([(0,oq), pt2, (0,op)])
        b_cycles.append([pt1, pt2])
    return diagramFromCycleInfo("Identity cobordism for %s" % repr(pmc),
                                num_interior_point = n, length_border = [n,n],
                                alpha_arcs = a_arcs, beta_cycles = b_cycles,
                                crossing_orientation = [1,-1]*pmc.num_pair)

def getArcslideDiagram(slide):
    """Get Heegaard diagram for the given arcslide. The left (first) boundary
    has the opposite of starting pmc. The right (second) boundary has ending
    pmc.

    """
    pmc1, pmc2 = slide.start_pmc, slide.end_pmc
    b_pair, c_pair = slide.b_pair, slide.c_pair
    n, num_pair = pmc1.n, pmc1.num_pair
    a_arcs = []
    b_cycles = []
    for i in range(num_pair):
        if i not in (b_pair, c_pair):
            # Draw arcs for ordinary pairs
            p, q = pmc1.pairs[i]   # index of points at left, from bottom
            rp, rq = slide.to_r[p], slide.to_r[q] # index at right, from bottom
            assert rp < rq
            pt1, pt2 = 2*i, 2*i+1
            a_arcs.append([(1,rp), pt1, (1,rq)])
            op, oq = n-1-p, n-1-q   # index at left, from top
            a_arcs.append([(0,oq), pt2, (0,op)])
            b_cycles.append([pt1, pt2])
    # Now draw the B and C pairs.
    bpt1, bpt2 = 2*b_pair, 2*b_pair+1   # two ordinary points on B pair
    cpt1, cpt2 = 2*c_pair, 2*c_pair+1   # two ordinary points on C pair
    spt = 2*num_pair   # the one point special for this diagram
    b1, b2, c1, c2 = slide.b1, slide.b2, slide.c1, slide.c2 # left, from bottom
    rb1, rb2, rc1, rc2 = (slide.to_r[b1], slide.to_r[b2],
                          slide.to_r[c1], slide.to_r[c2])  # right, from bottom
    ob1, ob2, oc1, oc2 = n-1-b1, n-1-b2, n-1-c1, n-1-c2
    a_arcs.append([(1,rb1), bpt1, (1,rb2)])
    a_arcs.append([(0,ob2), bpt2, spt, (0,ob1)])
    a_arcs.append([(1,rc1), cpt1, (1,rc2)])
    a_arcs.append([(0,oc2), cpt2, (0,oc1)])
    b_cycles.append([bpt1, bpt2])
    if b1 > c1:   # special cycle on c2 -> special pair at bottom
        b_cycles.append([cpt1, cpt2, spt])
    else:
        b_cycles.append([cpt1, spt, cpt2])
    # Note orientation at spt is always -1
    return diagramFromCycleInfo("Diagram for %s" % repr(slide),
                                num_interior_point = n+1,
                                length_border = [n,n],
                                alpha_arcs = a_arcs, beta_cycles = b_cycles,
                                crossing_orientation = [1,-1]*num_pair+[-1])

def getSimpleCobordismDiagram(start_pmc, insert_pos):
    """Get Heegaard diagram for a simple cobordism (see SimpleCobordismDA in
    cobordismda.py).

    """
    n, num_pair = start_pmc.n, start_pmc.num_pair
    a_arcs, b_cycles = [], []
    for i in range(num_pair):
        pt1, pt2 = 2*i, 2*i+1
        p, q = start_pmc.pairs[i]
        rp, rq = p, q
        if p >= insert_pos:
            rp += 4
        if q >= insert_pos:
            rq += 4
        op, oq = n-1-p, n-1-q  # index at left, from top
        a_arcs.append([(1, rp), pt1, (1, rq)])
        a_arcs.append([(0, oq), pt2, (0, op)])
        b_cycles.append([pt1, pt2])
    # Now draw the cobordism part, using points 2*num_pair+(0, 1, 2)
    pt0, pt1, pt2 = [2*num_pair+i for i in (0, 1, 2)]
    a_arcs.append([(1, insert_pos), pt2, pt1, (1, insert_pos+2)])
    a_arcs.append([(1, insert_pos+1), pt0, (1, insert_pos+3)])
    b_cycles.append([pt0, pt1, pt2])
    return diagramFromCycleInfo(
        "Diagram for cobordism starting at %s inserting at position %s" % (
            start_pmc, insert_pos),
        num_interior_point = 2 * num_pair + 3,
        length_border = [n, n + 4],
        alpha_arcs = a_arcs, beta_cycles = b_cycles,
        crossing_orientation = [1, -1]*num_pair + [-1, 1, -1])

def getCobordismDiagramLeft(cob):
    """Get Heegaard diagram for a cobordism on the linear PMC, with larger PMC
    on the left. cob is of type Cobordism (cobordism.py)

    """
    genus_l = cob.genus
    genus_r = genus_l - 1
    num_pair_l = genus_l * 2
    start_pmc, end_pmc = cob.large_pmc, cob.small_pmc
    start_n = start_pmc.n
    pt_count = [0]
    a_arcs, b_cycles = [], []

    def process_pair(left_pair):
        # Given the index of a pair on the left, draw the same a_arcs and
        # b_cycles as in the identity diagram for that pair.
        p, q = start_pmc.pairs[left_pair]
        rp, rq = cob.to_s[p], cob.to_s[q]
        op, oq = start_n-1-p, start_n-1-q
        pt1, pt2 = pt_count[0], pt_count[0] + 1
        pt_count[0] += 2
        a_arcs.append([(1, rp), pt1, (1, rq)])
        a_arcs.append([(0, oq), pt2, (0, op)])
        b_cycles.append([pt1, pt2])

    if cob.is_degenerate:
        c_pair, p_pair = cob.c_pair, cob.p_pair
        for i in range(num_pair_l):
            if i not in (c_pair, p_pair):
                process_pair(i)
        # Now consider the c-pair and p-pair
        c1, c2 = start_pmc.pairs[c_pair]
        p1, p2 = start_pmc.pairs[p_pair]
        oc1, oc2, op1, op2 = [start_n-1-p for p in (c1, c2, p1, p2)]
        pt = pt_count[0]
        pt_count[0] += 1
        a_arcs.append([(0, oc2), (0, oc1)])
        a_arcs.append([(0, op2), pt, (0, op1)])
        b_cycles.append([pt])
        crossing_orientation = [1, -1] * (num_pair_l - 2) + [1]

    else:  # cob is not degenerate
        c_pair, d_pair, u_pair = cob.c_pair, cob.d_pair, cob.u_pair
        du_pair = cob.du_pair  # for the right side
        for i in range(num_pair_l):
            if i not in (c_pair, d_pair, u_pair):
                process_pair(i)
        # Now consider the remaining pairs
        c1, c2 = start_pmc.pairs[c_pair]
        d1, d2 = start_pmc.pairs[d_pair]
        u1, u2 = start_pmc.pairs[u_pair]
        rdu1, rdu2 = end_pmc.pairs[du_pair]
        oc1, oc2, od1, od2, ou1, ou2 = [start_n-1-p
                                        for p in (c1, c2, d1, d2, u1, u2)]
        pt0, pt1, pt2, pt3 = [pt_count[0] + i for i in (0, 1, 2, 3)]
        pt_count[0] += 4
        a_arcs.append([(0, oc2), (0, oc1)])
        a_arcs.append([(0, ou2), pt0, pt2, (0, ou1)])
        a_arcs.append([(0, od2), pt3, (0, od1)])
        a_arcs.append([(1, rdu1), pt1, (1, rdu2)])
        b_cycles.append([pt0, pt1])
        b_cycles.append([pt2, pt3])
        crossing_orientation = [1, -1] * (num_pair_l - 3) + [-1, 1, 1, -1]

    return diagramFromCycleInfo(
        "Diagram for cobordism on linear PMC starting at genus %s and "
        "reducing pair %s " % (genus_l, c_pair),
        num_interior_point = pt_count[0],
        length_border = [start_n, start_n - 4],
        alpha_arcs = a_arcs, beta_cycles = b_cycles,
        crossing_orientation = crossing_orientation)
