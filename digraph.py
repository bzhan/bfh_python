"""Handles things related to directed graphs."""

from identityaa import *
from ddstructure import *

class DiGraph:
    """Interface for a general directed graph."""
    def getOutEdges(self, node):
        """Get the list of outward edges from a node."""
        raise NotImplementedError("Get list of outward edges not implemented.")

class DiGraphNode():
    """A general node in a digraph."""
    def __init__(self, parent):
        """Every node needs a parent graph."""
        self.parent = parent

class DiGraphEdge():
    """A general edge in a digraph."""
    def __init__(self, source, target):
        """Every edge needs a source and target (nodes from the same digraph).

        """
        assert source.parent == target.parent
        self.source = source
        self.target = target

class ConcreteDiGraph(DiGraph):
    """Represents a directed graph as a list of nodes and lists of edges."""
    def __init__(self):
        """Initializes an empty concrete digraph."""
        self.nodes = []
        self.edges = dict() # dictionary from nodes to list of edges

    def __str__(self):
        result = "Di-graph with %d nodes." % len(self.nodes)
        for node, outs in self.edges.items():
            node_str = "Node %s with outward edges:\n" % str(node)
            node_str += "\n".join([str(edge) for edge in outs])
            result += ("\n" + node_str)
        return result

    def addNode(self, node):
        """Add a node."""
        self.nodes.append(node)
        assert node not in self.edges
        self.edges[node] = []
        assert node.parent == self

    def addEdge(self, edge):
        """Add an edge between two nodes in this graph."""
        self.edges[edge.source].append(edge)

    def getNodes(self):
        return self.nodes
        
    def getOutEdges(self, node):
        return self.edges[node]

class TypeDGraphNode(DiGraphNode):
    """A node in a type D graph. Stores the generator of type D structure
    corresponding to this node.

    """
    def __init__(self, parent, dgen):
        """Specifies generator of type D structure corresponding to this node.

        """
        DiGraphNode.__init__(self, parent)
        self.dgen = dgen
        self.idem = self.dgen.idem

    def __str__(self):
        return str(self.dgen)

    def __repr__(self):
        return str(self)

    def __eq__(self, other):
        return self.dgen == other.dgen

    def __ne__(self, other):
        return not (self == other)

    def __hash__(self):
        return hash((self.dgen, "TypeDGraphNode"))

class TypeDGraphEdge(DiGraphEdge):
    """An edge in a type D graph. Corresponding to a type D operation. Stores
    the algebra coefficient.

    """
    def __init__(self, source, target, coeff):
        """Specifies the algebra coefficient."""
        DiGraphEdge.__init__(self, source, target)
        self.coeff = coeff

    def __str__(self):
        return "Edge from %s to %s with coeff %s" % \
            (str(self.source), str(self.target), str(self.coeff))

    def __repr__(self):
        return str(self)

class TypeDGraph(ConcreteDiGraph):
    """Corresponding to a type D structure. Nodes corresponds to generators and
    edges corresponds to type D operations.

    """
    def __init__(self, dstr):
        """Creates a type D graph from a type D structure."""
        ConcreteDiGraph.__init__(self)
        # Maintain a dictionary from DGenerator to nodes in graph
        self.graph_node = dict()
        self.algebra = dstr.algebra

        # Add nodes
        for dgen in dstr.getGenerators():
            cur_node = TypeDGraphNode(self, dgen)
            self.addNode(cur_node)
            self.graph_node[dgen] = cur_node
        # Add edges
        for gen_from in dstr.getGenerators():
            for (alg_coeff, gen_to), ring_coeff in gen_from.delta().items():
                self.addEdge(TypeDGraphEdge(
                        self.graph_node[gen_from], self.graph_node[gen_to],
                        alg_coeff))

class UniversalDiGraphNode(DiGraphNode, tuple):
    """A node in the universal digraph of an algebra. A sequence of algebra
    generators.

    """
    def __new__(cls, parent, data):
        return tuple.__new__(cls, data)

    def __init__(self, parent, data):
        """Specifies parent DiGraph. data is the list of generators."""
        self.parent = parent
        # Note tuple initialization is automatic
    
class UniversalDiGraph(DiGraph):
    """The universal digraph of an algebra (usually strand algebra of a PMC.)
    Nodes correspond to ordered sequences of algebra generators. From each
    node the set of outward edges is the set of algebra generators, appending
    that generator onto the sequence.

    """
    def __init__(self, algebra):
        """Create the universal digraph for the given algebra."""
        self.algebra = algebra

    def getOutEdges(self, gen_from):
        result = []
        for alg_gen in self.algebra.getGenerators():
            if not alg_gen.isIdempotent():
                gen_to = UniversalDiGraphNode(self, gen_from + (alg_gen,))
                result.append(TypeDGraphEdge(gen_from, gen_to, alg_gen))
        return result

    def getInitialNode(self):
        """Return the starting node, consisting of the empty sequence of
        generators.

        """
        return UniversalDiGraphNode(self, tuple())

class TypeDDGraphNode(DiGraphNode):
    """A node in a type DD graph. Stores the generator of type DD structure as
    well as a strand diagram whose right idempotent agrees with the left
    idempotent of that generator.

    """
    def __init__(self, parent, ddgen, sd):
        """Specifies generator of type DD structure and the strand diagram."""
        DiGraphNode.__init__(self, parent)
        self.ddgen = ddgen
        self.sd = sd
        self.idem1, self.idem2 = ddgen.idem1, ddgen.idem2

    def __str__(self):
        return "%s,%s" % (str(self.ddgen), str(self.sd))

    def __repr__(self):
        return str(self)

    def __eq__(self, other):
        return self.ddgen == other.ddgen

    def __ne__(self, other):
        return not (self == other)

    def __hash__(self):
        return hash((self.ddgen, self.sd, "TypeDDGraphNode"))

TypeDDGraphEdge = TypeDGraphEdge

class TypeDDGraph(ConcreteDiGraph):
    """Corresponding to a type DD structure. Nodes correspond to generators
    plus a strand diagram, edges correspond to type DD operations. tensor_side
    specifies which algebra action (the one by algebra1 or algebra2) is
    tensored with the type A or AA structure. The other side gives type D
    operations.

    """
    def __init__(self, ddstr, tensor_side):
        """Creates a type DD graph from a type DD structure."""
        ConcreteDiGraph.__init__(self)
        self.tensor_side = tensor_side
        # Maintain a dictionary from (DDGenerator, StrandDiagram) to nodes in
        # graph.
        self.graph_node = dict()
        # Dictionary from DDGenerator to directly corresponding node in graph
        # (the one with idempotent strand diagram).
        self.ddgen_node = dict()
        self.algebra1, self.algebra2 = ddstr.algebra1, ddstr.algebra2

        if tensor_side == 2:
            alg_gens = self.algebra1.getGenerators()
        else:
            alg_gens = self.algebra2.getGenerators()
        # Add nodes
        for ddgen in ddstr.getGenerators():
            for agen in alg_gens:
                if tensor_side == 2:
                    idem_to_match = ddgen.idem1
                else:
                    idem_to_match = ddgen.idem2
                if agen.getRightIdem() == idem_to_match:
                    cur_node = TypeDDGraphNode(self, ddgen, agen)
                    self.addNode(cur_node)
                    self.graph_node[(ddgen, agen)] = cur_node
                    if agen.isIdempotent():
                        self.ddgen_node[ddgen] = cur_node

    @memorize
    def getOutEdges(self, gen_from):
        result = []
        x, sd = gen_from.ddgen, gen_from.sd
        for (a, b, y), ring_coeff in x.delta().items():
            if self.tensor_side == 2:
                output_a, tensor_a = a, b
            else:
                output_a, tensor_a = b, a
            new_sd = sd * output_a
            if new_sd != E0:
                new_sd = new_sd.getElt()
                gen_to = self.graph_node[(y, new_sd)]
                result.append(TypeDDGraphEdge(gen_from, gen_to, tensor_a))
        return result

class TypeAAGraphNode(DiGraphNode, tuple):
    """A node in a type AA graph. Stores two strand diagrams of the same PMC.

    """
    def __new__(cls, parent, sd_left, sd_right):
        return tuple.__new__(cls, (sd_left, sd_right))

    def __init__(self, parent, sd_left, sd_right):
        """Specifies the two strand diagrams. Compute idempotents as the right
        idempotents of the strand diagrams.

        """
        # Note tuple initialization is automatic
        DiGraphNode.__init__(self, parent)
        self.idem1 = sd_left.getRightIdem()
        self.idem2 = sd_right.getRightIdem()

    def isHomology(self):
        """Returns whether this node is a homology node (consisting of
        idempotents).

        """
        return self[0].isIdempotent() and self[1].isIdempotent()

class TypeAAGraphEdge(DiGraphEdge):
    """An edge in a type AA graph. Stores type of edge and the algebra
    coefficient (if applicable).

    """
    # Possible values of edge_type
    ALG_LEFT, ALG_RIGHT, HOMOTOPY = range(3)

    def __init__(self, source, target, edge_type, coeff = None):
        """Specifies type of edge and algebra coefficient. coeff should be None
        if and only if edge_type is HOMOTOPY.

        """
        DiGraphEdge.__init__(self, source, target)
        self.edge_type = edge_type
        self.coeff = coeff

    def __str__(self):
        def strEdgeType(edge_type):
            if edge_type == self.ALG_LEFT: return "Left algebra"
            elif edge_type == self.ALG_RIGHT: return "Right algebra"
            else: return "Homotopy"
        result = "%s edge from %s to %s" % \
            (strEdgeType(self.edge_type), self.source, self.target)
        if self.edge_type != self.HOMOTOPY:
            result += " with coefficient %s" % str(self.coeff)
        return result

    def __repr__(self):
        return str(self)

class ATensorDGenerator(Generator, tuple):
    """Generator of a chain complex formed by tensoring a type A structure and
    a type D structure (actually, the result of tensoring D1 * CFAA(Id) * D2).

    """
    def __new__(cls, parent, gen_left, gen_right):
        return tuple.__new__(cls, (gen_left, gen_right))

    def __init__(self, parent, gen_left, gen_right):
        """Specify generators on the two sides of the tensor. Both are type D
        generators because the AA generator in between is assumed.

        """
        # Note tuple initialization is automatic
        Generator.__init__(self, parent)

class DATensorDGenerator(DGenerator, tuple):
    """Generator of a type D structure formed by tensoring a type DA structure
    and a type D structure (actually, the result of tensoring
    DD * CFAA(Id) * D).

    """
    def __new__(cls, parent, gen_left, gen_right):
        return tuple.__new__(cls, (gen_left, gen_right))

    def __init__(self, parent, gen_left, gen_right):
        """Specify generators on two sides of the tensor (DD and D generators).

        """
        # Note tuple initialization is automatic
        DGenerator.__init__(self, parent, gen_left.idem1)

class ATensorDDGenerator(DGenerator, tuple):
    """Generator of a type D structure formed by tensoring a type A structure
    and a type DD structure (actually, the result of tensoring
    D * CFAA(Id) * DD).

    """
    def __new__(cls, parent, gen_left, gen_right):
        return tuple.__new__(cls, (gen_left, gen_right))

    def __init__(self, parent, gen_left, gen_right):
        """Specify generators on two sides of the tensor (D and DD generators).

        """
        # Note tuple initialization is automatic
        DGenerator.__init__(self, parent, gen_right.idem2)

class DATensorDDGenerator(DDGenerator, tuple):
    """Generator of a type DD structure formed by tensoring a type DA structure
    and a type DD structure (actually, the result of tensoring
    DD * CFAA(Id) * DD).

    """
    def __new__(cls, parent, gen_left, gen_right):
        return tuple.__new__(cls, (gen_left, gen_right))

    def __init__(self, parent, gen_left, gen_right):
        """Specify generators on two sides of the tensor (both DD generators).

        """
        # Note tuple initialization is automatic
        DDGenerator.__init__(self, parent, gen_left.idem1, gen_right.idem2)

class TypeAAGraph(DiGraph):
    """Digraph used for simplifying a type AA structure. Nodes correspond to
    generators of the large chain complex. Edges can be either homotopies or
    multiplications on either side of the generator.

    """
    ALG_LEFT = TypeAAGraphEdge.ALG_LEFT
    ALG_RIGHT = TypeAAGraphEdge.ALG_RIGHT
    HOMOTOPY = TypeAAGraphEdge.HOMOTOPY
    
    def __init__(self, pmc):
        """Creates a type AA graph simplifying the type AA of identity for the
        given PMC. Uses the chain complex and homotopy calculated in class
        HomotopyAA.

        """
        self.pmc_alg = pmc.getAlgebra()

        # Dictionary mapping idempotents to the homology node whose idem1 is
        # that idempotent.
        self.homology_node = dict()
        for idem in pmc.getIdempotents():
            self.homology_node[idem] = TypeAAGraphNode(
                self, idem.toAlgElt(self.pmc_alg),
                idem.comp().toAlgElt(self.pmc_alg))

    @memorize
    def algFactor(self, sd1, sd2):
        """Factor out sd2 from sd1 on the right. Attempt to find sd such that
        sd * sd2 = sd1. Returns sd if it exists. Otherwise return None.

        """
        if sd2.getRightIdem() != sd1.getRightIdem():
            return None
        # Remove strands in sd2 one by one. Keep track of remaining strands.
        st_remain = list(sd1.strands)
        for a, b in sd2.strands:
            found = False
            for i in range(len(st_remain)):
                p, q = st_remain[i]
                if q == b and p <= a:
                    found = True
                    st_remain.remove((p, q))
                    if p < a:
                        st_remain.append((p, a))
                    break
            if not found:
                return None
        if not Strands(sd1.pmc, st_remain).rightCompatible(sd2.getLeftIdem()):
            return None
        result = StrandDiagram(sd1.parent, None, st_remain, sd2.getLeftIdem())
        if result * sd2 == 1 * sd1:
            return result
        else:
            return None

    @memorize
    def getHomotopyEdges(self, source):
        """Get the list of homotopy edges starting at source."""
        return [TypeAAGraphEdge(source,
                                TypeAAGraphNode(self, *target), self.HOMOTOPY)
                for target in homotopyMap(*source)]

    def getAlgLeftTarget(self, source, coeff):
        """If there is a left algebra action edge starting at source and with
        coefficient coeff, return the target node of that edge. Otherwise,
        return None.

        """
        factor = self.algFactor(source[0], coeff)
        if factor is None:
            return None
        return TypeAAGraphNode(self, factor, source[1])

    def getAlgRightTarget(self, source, coeff):
        """If there is a right algebra action edge starting at source and with
        coefficient coeff, return the target node of that edge. Otherwise,
        return None.

        """
        prod = source[1] * coeff
        if prod == E0:
            return None
        prod = prod.getElt()
        return TypeAAGraphNode(self, source[0], prod)

    def _searchDoubleD(self, d_graph1, d_graph2, start_pos):
        """Search for paths in d_graph1*self*d_graph2 with the given list of
        starting positions (each element of start_pos is a tuple (d1_pos,
        d2_pos, aa_pos)). Returns a list of lists of end positions.

        """
        def search(d1_pos, d2_pos, aa_pos, is_homotopy):
            """Helper function performing a one step search, starting at the
            given locations in the three graphs. If is_homotopy is set, the
            next move must be a homotopy. Otherwise, it must be an algebra
            action (on either side). Returns the list of end states with aa_pos
            at a homology node.

            """
            result = []
            if is_homotopy:
                if aa_pos.isHomology():
                    result.append((d1_pos, d2_pos, aa_pos))
                else:
                    for edge in self.getHomotopyEdges(aa_pos):
                        result += search(d1_pos, d2_pos, edge.target, False)
            else:
                for d1_edge in d_graph1.getOutEdges(d1_pos):
                    target = self.getAlgLeftTarget(aa_pos, d1_edge.coeff.opp())
                    if target is not None:
                        result += search(d1_edge.target, d2_pos, target, True)
                for d2_edge in d_graph2.getOutEdges(d2_pos):
                    target = self.getAlgRightTarget(aa_pos, d2_edge.coeff)
                    if target is not None:
                        result += search(d1_pos, d2_edge.target, target, True)
            return result

        full_result = []
        for d1_pos, d2_pos, aa_pos in start_pos:
            full_result.append(search(d1_pos, d2_pos, aa_pos, False))
        return full_result

    def tensorDoubleD(self, d_graph1, d_graph2):
        """Computes the chain complex D1 * CFAA(Id) * D2, where D1, D2 are type
        D structures with graphs d_graph1, d_graph2, and CFAA(Id) is
        represented by this graph. Both D1 and D2 are left type D structures.
        The algebra acting on D1 is opposite of self.pmc_alg, and the algebra
        acting on D2 is the same as self.pmc_alg

        """
        assert d_graph1.algebra.opp() == self.pmc_alg
        assert d_graph2.algebra == self.pmc_alg
        cx = SimpleChainComplex(F2)
        # Generators of the chain complex:
        for node1 in d_graph1.getNodes():
            for node2 in d_graph2.getNodes():
                if node1.idem == node2.idem.opp().comp():
                    cur_gen = ATensorDGenerator(cx, node1.dgen, node2.dgen)
                    cx.addGenerator(cur_gen)

        # Search the graphs for edges in the chain complex
        for gen_start in cx.getGenerators():
            dgen1, dgen2 = gen_start
            d1_pos = d_graph1.graph_node[dgen1]
            d2_pos = d_graph2.graph_node[dgen2]
            aa_pos = self.homology_node[dgen1.idem.opp()]
            pos = [(d1_pos, d2_pos, aa_pos)]
            end_states = self._searchDoubleD(d_graph1, d_graph2, pos)[0]
            for d1_end, d2_end, aa_end in end_states:
                gen_end = ATensorDGenerator(cx, d1_end.dgen, d2_end.dgen)
                cx.addDifferential(gen_start, gen_end, 1)
        return cx

    def tensorDDandD(self, dd_graph, d_graph):
        """Computes the type D structure DD1 * CFAA(Id) * D2, where DD1 is a
        type DD structure with graph dd_graph, CFAA(Id) is represented by this
        graph, and D2 is a type D structure with graph d_graph.

        """
        assert dd_graph.tensor_side == 2
        assert dd_graph.algebra2.opp() == self.pmc_alg
        assert d_graph.algebra == self.pmc_alg
        dstr = SimpleDStructure(F2, dd_graph.algebra1)
        # Generators of the type D structure:
        for ddgen, node1 in dd_graph.ddgen_node.items():
            for node2 in d_graph.getNodes():
                if node1.idem2 == node2.idem.opp().comp():
                    cur_gen = DATensorDGenerator(dstr, ddgen, node2.dgen)
                    dstr.addGenerator(cur_gen)

        # Search the graphs for type D operations
        for gen_start in dstr.getGenerators():
            ddgen, dgen = gen_start
            d1_pos = dd_graph.ddgen_node[ddgen]
            d2_pos = d_graph.graph_node[dgen]
            aa_pos = self.homology_node[ddgen.idem2.opp()]
            pos = [(d1_pos, d2_pos, aa_pos)]
            end_states = self._searchDoubleD(dd_graph, d_graph, pos)[0]
            for d1_end, d2_end, aa_end in end_states:
                gen_end = DATensorDGenerator(dstr, d1_end.ddgen, d2_end.dgen)
                dstr.addDelta(gen_start, gen_end, d1_end.sd, 1)
        return dstr

    def tensorDandDD(self, d_graph, dd_graph):
        """Computes the type D structure D1 * CFAA(Id) * DD2, where D1 is a
        type D structure with graph d_graph, CFAA(Id) is represented by this
        graph, and DD2 is a type DD structure with graph dd_graph.

        """
        assert dd_graph.tensor_side == 1
        assert d_graph.algebra.opp() == self.pmc_alg
        assert dd_graph.algebra1 == self.pmc_alg
        dstr = SimpleDStructure(F2, dd_graph.algebra2)
        # Generators of the type D structure:
        for node1 in d_graph.getNodes():
            for ddgen, node2 in dd_graph.ddgen_node.items():
                if node1.idem == node2.idem1.opp().comp():
                    cur_gen = ATensorDDGenerator(dstr, node1.dgen, ddgen)
                    dstr.addGenerator(cur_gen)

        # Search the graphs for type D operations
        for gen_start in dstr.getGenerators():
            dgen, ddgen = gen_start
            d1_pos = d_graph.graph_node[dgen]
            d2_pos = dd_graph.ddgen_node[ddgen]
            aa_pos = self.homology_node[dgen.idem.opp()]
            pos = [(d1_pos, d2_pos, aa_pos)]
            end_states = self._searchDoubleD(d_graph, dd_graph, pos)[0]
            for d1_end, d2_end, aa_end in end_states:
                gen_end = ATensorDDGenerator(dstr, d1_end.dgen, d2_end.ddgen)
                dstr.addDelta(gen_start, gen_end, d2_end.sd, 1)
        return dstr

    def tensorDoubleDD(self, dd_graph1, dd_graph2):
        """Compute the type DD structure DD1 * CFAA(Id) * DD2."""
        assert dd_graph1.tensor_side == 2 and dd_graph2.tensor_side == 1
        assert dd_graph1.algebra2.opp() == self.pmc_alg
        assert dd_graph2.algebra1 == self.pmc_alg
        ddstr = SimpleDDStructure(F2, dd_graph1.algebra1, dd_graph2.algebra2)
        # Generators of the type DD structure:
        for ddgen1, node1 in dd_graph1.ddgen_node.items():
            for ddgen2, node2 in dd_graph2.ddgen_node.items():
                if node1.idem2 == node2.idem1.opp().comp():
                    cur_gen = DATensorDDGenerator(ddstr, ddgen1, ddgen2)
                    ddstr.addGenerator(cur_gen)

        # Search the graphs for type DD operations
        for gen_start in ddstr.getGenerators():
            ddgen1, ddgen2 = gen_start
            d1_pos = dd_graph1.ddgen_node[ddgen1]
            d2_pos = dd_graph2.ddgen_node[ddgen2]
            aa_pos = self.homology_node[ddgen1.idem2.opp()]
            pos = [(d1_pos, d2_pos, aa_pos)]
            end_states = self._searchDoubleD(dd_graph1, dd_graph2, pos)[0]
            for d1_end, d2_end, aa_end in end_states:
                gen_end = DATensorDDGenerator(ddstr,
                                              d1_end.ddgen, d2_end.ddgen)
                ddstr.addDelta(gen_start, gen_end, d1_end.sd, d2_end.sd, 1)
        return ddstr

    def getTypeDA(self, dd_graph):
        """Returns the type DA structure formed by tensoring dd_graph with
        self. For now self is placed at left, and dd_graph is placed at right.

        """
        assert dd_graph.tensor_side == 1
        assert dd_graph.algebra1 == self.pmc_alg
        univ_digraph = UniversalDiGraph(self.pmc_alg.opp())
        for ddgen, d2_pos in dd_graph.ddgen_node.items():
            d1_pos = univ_digraph.getInitialNode()
            aa_pos = self.homology_node[ddgen.idem1.comp()]
            pos = [(d1_pos, d2_pos, aa_pos)]
            end_states = self._searchDoubleD(univ_digraph, dd_graph, pos)[0]
            for d1_end, d2_end, aa_end in end_states:
                print d1_end, d2_end.sd

@memorize
def getTypeAAGraph(pmc):
    """Memorized version of constructor of TypeAAGraph."""
    return TypeAAGraph(pmc)

def computeATensorD(dstr1, dstr2):
    """Compute the tensor product dstr1 * CFAA(Id) * dstr2. dstr1 and dstr2 are
    left type D structures over opposite PMC's. Result is a chain complex.

    """
    aa_graph = getTypeAAGraph(dstr2.algebra.pmc)
    cx = aa_graph.tensorDoubleD(TypeDGraph(dstr1), TypeDGraph(dstr2))
    # Add gradings if necessary
    if hasattr(dstr1, "gr_set") and hasattr(dstr2, "gr_set"):
        cx.gr_set = GeneralGradingSet([dstr1.gr_set.Ropp(), dstr2.gr_set])
        cx.grading = dict()
        for gen in cx.getGenerators():
            dgen1, dgen2 = gen
            cx.grading[gen] = GeneralGradingSetElement(
                cx.gr_set, [dstr1.grading[dgen1].Ropp(), dstr2.grading[dgen2]])
    return cx

def computeDATensorD(ddstr1, dstr2):
    """Compute the tensor product ddstr1 * CFAA(Id) * dstr2. ddstr1 is a left
    type DD structure and dstr2 is a left type D structure. The algebra of the
    second action of ddstr1 is opposite to the algebra of the action on dstr2.

    """
    aa_graph = getTypeAAGraph(dstr2.algebra.pmc)
    dstr = aa_graph.tensorDDandD(TypeDDGraph(ddstr1, 2), TypeDGraph(dstr2))
    # Add gradings if necessary
    if hasattr(ddstr1, "gr_set") and hasattr(dstr2, "gr_set"):
        dstr.gr_set = GeneralGradingSet(
            [ddstr1.gr_set.partialRopp(1), dstr2.gr_set])
        dstr.grading = dict()
        for dgen in dstr.getGenerators():
            ddgen1, dgen2 = dgen
            dstr.grading[dgen] = GeneralGradingSetElement(
                dstr.gr_set, [ddstr1.grading[ddgen1].partialRopp(1),
                              dstr2.grading[dgen2]])
    return dstr

def computeATensorDD(dstr1, ddstr2):
    """Compute the tensor product dstr1 * CFAA(Id) * ddstr2. dstr1 is a left
    type D structure and ddstr2 is a left type DD structure. The algebra of the
    action on dstr1 is the opposite to the algebra of the first action on
    dstr2.

    """
    aa_graph = getTypeAAGraph(ddstr2.algebra1.pmc)
    dstr = aa_graph.tensorDandDD(TypeDGraph(dstr1), TypeDDGraph(ddstr2, 1))
    # Add gradings if necessary
    if hasattr(dstr1, "gr_set") and hasattr(ddstr2, "gr_set"):
        dstr.gr_set = GeneralGradingSet(
            [ddstr2.gr_set.partialRopp(0), dstr1.gr_set])
        dstr.grading = dict()
        for dgen in dstr.getGenerators():
            dgen1, ddgen2 = dgen
            dstr.grading[dgen] = GeneralGradingSetElement(
                dstr.gr_set, [ddstr2.grading[ddgen2].partialRopp(0),
                              dstr1.grading[dgen1]])
    return dstr

def computeDATensorDD(ddstr1, ddstr2):
    """Compute the tensor product ddstr1 * CFAA(Id) * ddstr2. Both ddstr1 and
    ddstr2 are left type DD structures. The algebra of the second action on
    ddstr1 is the opposite to the algebra of the first action on ddstr2.

    """
    aa_graph = getTypeAAGraph(ddstr2.algebra1.pmc)
    ddstr = aa_graph.tensorDoubleDD(TypeDDGraph(ddstr1, 2),
                                    TypeDDGraph(ddstr2, 1))
    # Add gradings if necessary
    if hasattr(ddstr1, "gr_set") and hasattr(ddstr2, "gr_set"):
        ddstr.gr_set = GeneralGradingSet(
            [ddstr1.gr_set.partialRopp(1), ddstr2.gr_set])
        ddstr.grading = dict()
        for ddgen in ddstr.getGenerators():
            ddgen1, ddgen2 = ddgen
            ddstr.grading[ddgen] = GeneralGradingSetElement(
                ddstr.gr_set, [ddstr1.grading[ddgen1].partialRopp(1),
                               ddstr2.grading[ddgen2]])
    return ddstr
