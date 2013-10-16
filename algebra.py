"""Definitions of core algebraic objects.

Some design decisions:

Each generator can belong to only one chain complex or module. This allows for
the possibility of ``a.diff()`` where ``a`` is a Generator.

"""

import heapq
from numbers import Number
from utility import NamedObject, SummableDict
from utility import fracToInt, memorize, memorizeHash, safeMultiply
from utility import F2

class FreeModule:
    """Represents a free module over some ring."""
    def __init__(self, ring):
        """Specifies the ring. Should be either a subclass of Ring or a python
        Number.

        """
        self.ring = ring

    def getGenerators(self):
        """Returns a list of all generators. Need not be implemented for every
        module.

        """
        raise NotImplementedError("Get set of generators not implemented.")

class Generator:
    """Represents a generator of a free module. By default, generators are
    distinguished by (python) identity. Implement __eq__, __ne__, and __hash__
    for custom behavior.

    """

    def __init__(self, parent):
        """Only information that every generator needs is the parent module."""
        self.parent = parent

    def elt(self, coeff = 1):
        """Returns the element coeff * self."""
        return self.ELT_CLASS({self : coeff})

    def diff(self):
        """Returns the differential of this generator. Make sense only if
        parent module implements diff().

        """
        return self.parent.diff(self)

    def antiDiff(self):
        """Computes the dual of differential. Make sense only if parent module
        implements antiDiff().

        """
        return self.parent.antiDiff(self)

    def factor(self):
        """Find all ways to factor this generator into a product of two
        generators. Make sense only if parent module implements factor().

        """
        return self.parent.factor(self)

    def delta(self):
        """Returns the delta of this generator (for type D and DD structures).
        Make sense only if parent module implements delta().

        """
        return self.parent.delta(self)

    def __mul__(self, other):
        """Multiplies this generator by ``other`` on the left. Can be either
        algebra multiplication or module action. Usually returns an element
        rather than a generator.

        """
        if isinstance(other, Number):
            return self.elt(other)
        elif hasattr(self.parent, "multiply"):
            return self.parent.multiply(self, other)
        else:
            return NotImplemented

    def __rmul__(self, other):
        """Multiplies this generator by ``other`` on the right. Usually
        represents a left module action (with ``other`` an algebra generator).

        """
        if isinstance(other, Number):
            return self.elt(other)
        elif hasattr(self.parent, "rmultiply"):
            return self.parent.rmultiply(self, other)
        else:
            return NotImplemented

    def toSimpleGenerator(self, name):
        """Convert to a SimpleGenerator with the given name. All fields are
        preserved, except ``name`` which is overwritten, and _hash_val which is
        removed, if present.

        """
        new_obj = SimpleGenerator(self.parent, name)
        new_obj.__dict__.update(self.__dict__)
        new_obj.name = name # to make sure original name is overwritten
        if hasattr(new_obj, '_hash_val'):
            del new_obj._hash_val # reset hash value
        return new_obj

class SimpleGenerator(Generator, NamedObject):
    """Generator has a name. Distinguished by name."""
    def __init__(self, parent, name):
        """Specifies name and parent module."""
        Generator.__init__(self, parent)
        NamedObject.__init__(self, name)

class Element(SummableDict):
    """Represents an element of a free module, as a dictionary from generators
    to coefficients. For example: a+2b will be represented as {a:1,b:2}.

    """
    def __init__(self, data = None):
        """Corrects type of coefficients if necessary."""
        if data is None:
            data = {}
        SummableDict.__init__(self, data)
        if self:
            convert = self.getElt().parent.ring.convert
            for key, value in self.items():
                self[key] = convert(value)

    def __str__(self):
        if self == 0:
            return "0"
        terms = []
        for gen, coeff in self.items():
            if coeff == 1:
                terms.append(str(gen))
            else:
                terms.append(str(coeff)+"*"+str(gen))
        return "+".join(terms)

    def __repr__(self):
        return str(self)

    def __mul__(self, other):
        # First try multiplying each coefficient with other, using the function
        # in SummableDict.
        result = SummableDict.__mul__(self, other)
        if result != NotImplemented:
            return result

        # Now try to multiply each key by other on the left.
        result = E0
        for k, v in self.items():
            prod = safeMultiply(k, other)
            if prod is NotImplemented:
                return NotImplemented
            result += [term * (v * coeff) for term, coeff in prod.items()]
        return result

    def __rmul__(self, other):
        # First try multiplying each coefficient with other, using the function
        # in SummableDict.
        result = SummableDict.__rmul__(self, other)
        if result != NotImplemented:
            return result

        # Now try to multiply key by other on the left.
        result = E0
        for k, v in self.items():
            prod = safeMultiply(other, k)
            if prod is NotImplemented:
                return NotImplemented
            result += [term * (v * coeff) for term, coeff in prod.items()]
        return result

    def diff(self):
        """Returns the differential of this element."""
        return sum([coeff * gen.diff() for gen, coeff in self.items()], E0)

# Name of the class for elements containing this generator
Generator.ELT_CLASS = Element
# Short-hand for empty element
E0 = Element()

class ChainComplex(FreeModule):
    """Represents a general chain complex."""
    def diff(self, gen):
        """Returns the differential of a generator. """
        raise NotImplementedError("Differential not implemented.")

    @memorize
    def _getAntiDiffMap(self):
        """Helper function generating tables of dual of differential, for calls
        to antiDiff.

        """
        gen_list = self.getGenerators()
        antiDiffMap = {}
        for gen in gen_list:
            antiDiffMap[gen] = E0
        for gen in gen_list:
            for dgen, coeff in gen.diff().items():
                antiDiffMap[dgen] += coeff * gen
        return antiDiffMap

    def antiDiff(self, gen):
        """Returns the dual of the differential of gen, as an element of this
        algebra. The element is a sum of all terms c*y, for which gen appears
        in the dy with coefficient c.
        By default, need getGenerators() to be implemented. antiDiff of all
        generators is computed at once.

        """
        return self._getAntiDiffMap()[gen]

class SimpleChainComplex(ChainComplex):
    """Represents a chain complex with a finite number of generators, with
    explicitly stored generating set and differential. The generating set is
    stored as a python set, and differential is stored as a dictionary mapping
    from generators to elements. Each generator must be a key in the dictionary
    (even if its differential is zero).

    """
    def __init__(self, ring):
        """Initialize an empty chain complex."""
        ChainComplex.__init__(self, ring)
        self.generators = set()
        self.differential = dict()

    def __str__(self):
        result = "Chain complex.\n"
        for k, v in self.differential.items():
            result += "d(%s) = %s\n" % (k, v)
        return result

    def __repr__(self):
        return str(self)

    def __len__(self):
        return len(self.generators)

    def diff(self, generator):
        return self.differential[generator]

    def getGenerators(self):
        return list(self.generators)

    def reindex(self):
        """Replace the generators by simple generators indexed by integers. The
        names of the new generators are 'g1', 'g2', etc.

        """
        gen_list = list(self.generators)
        new_gen_list = []
        # Dictionary mapping original generators to new ones
        translate_dict = dict()
        for i in range(len(gen_list)):
            new_gen = gen_list[i].toSimpleGenerator("g%d"%(i+1))
            new_gen_list.append(new_gen)
            translate_dict[gen_list[i]] = new_gen
        self.generators = set(new_gen_list)
        new_diff = dict()
        for gen, dgen in self.differential.items():
            new_diff[translate_dict[gen]] = dgen.translateKey(translate_dict)
        self.differential = new_diff
        if hasattr(self, "grading"):
            new_grading = dict()
            for gen, gr in self.grading.items():
                if gen in translate_dict: # gen is still in chain complex
                    new_grading[translate_dict[gen]] = gr
            self.grading = new_grading

    def addGenerator(self, generator):
        """Add a generator. No effect if the generator already exists."""
        assert generator.parent == self
        self.generators.add(generator)
        if not self.differential.has_key(generator):
            self.differential[generator] = E0

    def addDifferential(self, gen_from, gen_to, coeff):
        """Add coeff * gen_to to the differential of gen_from. Both gen_from
        and gen_to should be generators of this complex.

        """
        assert gen_from.parent == self and gen_to.parent == self
        self.differential[gen_from] += coeff * gen_to

    def simplify(self, find_homology_basis = False):
        """Simplify a chain complex using cancellation lemma."""
        # Build dictionary of coefficients
        arrows = dict()
        for x in self.generators:
            arrows[x] = dict()
        for x in self.generators:
            for y, coeff in self.differential[x].items():
                arrows[x][y] = coeff

        arrows = simplifyComplex(arrows, F2.zero, find_homology_basis)

        # Now rebuild the chain complex
        self.generators = set(arrows.keys())
        self.differential = dict()
        for x in arrows:
            self.differential[x] = E0
            for y, coeff in arrows[x].items():
                self.differential[x] += coeff * y

    def checkDifferential(self):
        """Checks the relation d^2 for differentials."""
        for gen in self.generators:
            assert gen.diff().diff() == 0

    def checkGrading(self):
        """Check grading is consistent with differentials."""
        for x in self.generators:
            for y, coeff in x.diff().items():
                assert self.grading[x] - 1 == self.grading[y]

    def getGradingInfo(self):
        """Shows the distribution of gradings in an easy-to-read format."""
        distr_by_spinc = dict()
        for x in self.generators:
            ref_grading = self.grading[x]
            break
        for x in self.generators:
            maslov, mod, spinc = \
                self.gr_set.eltDiffShortForm(self.grading[x], ref_grading)
            spinc = tuple(spinc)
            if spinc not in distr_by_spinc:
                distr_by_spinc[spinc] = []
            distr_by_spinc[spinc].append(maslov)
        distr_count = dict()
        for spinc, distr in distr_by_spinc.items():
            min_gr, max_gr = min(distr), max(distr)
            cur_count = [0] * fracToInt(max_gr - min_gr + 1)
            for gr in distr:
                cur_count[fracToInt(gr - min_gr)] += 1
            cur_count = tuple(cur_count)
            if cur_count not in distr_count:
                distr_count[cur_count] = 0
            distr_count[cur_count] += 1
        return distr_count

    def getAbsGradingInfo(self):
        """Returns the list of absolute gradings of generators."""
        return sorted([self.gr_set.eltAbsoluteGrading(self.grading[gen])
                       for gen in self.generators])

class SimpleChainMorphism:
    """Represents a morphism between two simple chain complexes (which may be
    the same). Need not be a chain map (so can be used to represent homotopy,
    for example). Represented explicitly.

    """
    def __init__(self, cx_from, cx_to):
        """Gives the source and target chain complexes of the morphism."""
        self.cx_from = cx_from
        self.cx_to = cx_to
        self.morphism = dict()

    def addMorphism(self, gen_from, gen_to, coeff):
        """gen_from is a generator of cx_from, and gen_to is a generator of
        cx_to. Adds coeff * gen_to to f(gen_from).

        """
        assert gen_from.parent == self.cx_from and gen_to.parent == self.cx_to
        if gen_from not in self.morphism:
            self.morphism[gen_from] = E0
        self.morphism[gen_from] += coeff * gen_to

    def apply(self, x):
        """Computes f(x), where x is either a generator or an element of
        cx_from. The returned value is always an Element of cx_to.

        """
        if isinstance(x, Generator):
            assert x.parent == self.cx_from
            if x not in self.morphism:
                return E0
            return self.morphism[x]
        else:
            assert isinstance(x, Element)
            return sum([coeff * self.apply(gen)
                        for gen, coeff in x.items()], E0)

    def __str__(self):
        result = "Morphism between two chain complexes.\n"
        for k, v in self.morphism.items():
            result += "f(%s) = %s\n" % (k, v)
        return result

    def __repr__(self):
        return str(self)

class DGAlgebra(ChainComplex):
    """Represents a general differential-graded algebra."""
    def multiply(self, gen1, gen2):
        """Returns the product of gen1 and gen2, as an algebra element."""
        raise NotImplementedError("Multiply not implemented.")

    @memorize
    def _getFactorMap(self):
        """Helper function generating tables for factoring generators, for
        calls to factor.

        """
        factorMap = {}
        gen_list = self.getGenerators()
        for gen in gen_list:
            factorMap[gen] = E0
        # Resulting element lies in the tensor product of A with itself
        parent = Tensor((self, self))
        for gen1 in gen_list:
            for gen2 in gen_list:
                for prod, coeff in (gen1*gen2).items():
                    tensor_gen = TensorGenerator((gen1, gen2), parent)
                    factorMap[prod] += coeff * tensor_gen
        return factorMap

    def factor(self, gen):
        """Returns an element of (A tensor A), where A is the present algebra.
        The element is the sum of all terms c*(p tensor q), for which gen is a
        term in p*q with coefficient c.
        By default, need getGenerators() to be implemented. factor of all
        generators is computed at once.

        """
        return self._getFactorMap()[gen]

class Tensor(FreeModule, tuple):
    """Represents a free module whose generating set is the product of the
    generating set of the two sides.

    """
    def __init__(self, data):
        """Specifies the left and right modules."""
        for d in data[1:]:
            assert d.ring == data[0].ring
        # Note tuple initialization is automatic
        FreeModule.__init__(self, data[0].ring)

class TensorGenerator(Generator, tuple):
    """Represents a generator of a free module that is a tensor product of two
    or more free modules. Works as a tuple of the components.

    """
    def __new__(cls, data, parent = None):
        return tuple.__new__(cls, tuple(data))

    def __init__(self, data, parent = None):
        """If parent is None, a default is used."""
        if parent is None:
            parent = Tensor(tuple([comp.parent for comp in data]))
        # Note tuple initialization is automatic
        Generator.__init__(self, parent)

    def __eq__(self, other):
        return tuple.__eq__(self, other) and self.parent == other.parent

    def __ne__(self, other):
        return not (self == other)

    @memorizeHash
    def __hash__(self):
        return hash((tuple(self), "Tensor"))

    def __str__(self):
        return "**".join(str(comp) for comp in self)

    def getLeftIdem(self):
        """Get the left idempotent. Only works if same function is implemented
        in each component.

        """
        return TensorIdempotent(tuple([comp.getLeftIdem() for comp in self]))

    def getRightIdem(self):
        """Get the right idempotent. Only works if same function is implemented
        in each part.

        """
        return TensorIdempotent(tuple([comp.getRightIdem() for comp in self]))

class TensorElement(Element):
    """Represents an element of the tensor product of two or more modules.

    TODO: Add support for quickly collecting terms by one of the components.

    """
    def __init__(self, data = None, parent = None):
        """If the keys are tuples, convert them to tensor generators."""
        if data is None:
            data = {}
        data_processed = {}
        for term, coeff in dict(data).items():
            if isinstance(term, TensorGenerator):
                data_processed[term] = coeff
            else:
                data_processed[TensorGenerator(term, parent)] = coeff
        Element.__init__(self, data_processed)

    def fixLast(self, gen, parent = None):
        """Collect terms with ``gen`` as the last factor. Return the
        coefficient as either Element or TensorElement.

        """
        result = E0
        for term, coeff in self.items():
            if term[-1] == gen:
                if len(term) > 2:
                    result += coeff * TensorGenerator(term[:-1], parent)
                else: # len(term) == 2
                    result += coeff * term[0]
        return result

    def invertible(self):
        """Tests whether this element is invertible."""
        for term, coeff in self.items():
            for comp in term:
                if not (1*comp).invertible():
                    return False
        return True

    def inverse(self):
        """Returns the inverse of this element, if invertible. Undefined
        behavior if the element is not invertible.

        """
        return self

TensorGenerator.ELT_CLASS = TensorElement

def expandTensor(prod, parent = None):
    """Produces the tensor element formed by the tensor product of either
    generators or elements.

    ``prod`` is a tuple of either Generator or Element, corresponding to the
    components of the tensor product.

    For example, ((1*A+1*B),C) expands into 1*(A,C)+1*(B,C), and
    ((1*A-1*B),(1*C-1*D)) expands into 1*(A,C)-1*(B,C)-1*(A,D)+1*(B,D).

    ``parent`` specifies the Tensor module. If it is set to None, the default
    (with no additional operations defined) will be used (during the
    initialization of TensorElement).

    """
    assert isinstance(prod, tuple)
    num_part = len(prod)
    expanded = [(prod, 1)]
    for i in range(num_part):
        if len(expanded) == 0:
            return E0
        if isinstance(expanded[0][0][i], Generator):
            continue
        expanded2 = []
        for subterm, coeff in expanded:
            for gen, coeff2 in subterm[i].items():
                expanded2.append(((subterm[0:i]+(gen,)+subterm[i+1:]),
                                 coeff*coeff2))
        expanded = expanded2
    if isinstance(parent, TensorStar):
        return TensorStarElement(dict(expanded), parent)
    else:
        return TensorElement(dict(expanded), parent)

class TensorDGAlgebra(Tensor, DGAlgebra):
    """Tensor product of DGAlgebras is a DGAlgebra."""
    def diff(self, gen):
        return E0.accumulate([
                expandTensor(gen[:i]+(gen[i].diff(),)+gen[i+1:], self)
                for i in range(len(gen))])

    def multiply(self, gen1, gen2):
        if not isinstance(gen1, TensorGenerator) or gen1.parent != self:
            return NotImplemented
        if not isinstance(gen2, TensorGenerator) or gen2.parent != self:
            return NotImplemented

        return expandTensor(tuple([gen1[i]*gen2[i] for i in range(len(self))]),
                            self)

    def getGenerators(self):
        """Return the set of generators. Use product of sets of generators of
        the components. Currently only implemented for tensors of two algebras.

        """
        if len(self) != 2:
            return NotImplemented

        gens1 = self[0].getGenerators()
        gens2 = self[1].getGenerators()
        result = []
        for gen1 in gens1:
            for gen2 in gens2:
                result.append(TensorGenerator((gen1, gen2), self))
        return result

class TensorIdempotent(tuple):
    """Serves as idempotent to a tensor product of algebras."""
    def toAlgElt(self, parent):
        """Get the algebra element corresponding to this idempotent."""
        assert isinstance(parent, TensorDGAlgebra)
        return TensorGenerator(tuple([self[i].toAlgElt(parent[i])
                                      for i in range(len(self))]), parent)

class TensorStarGenerator(Generator, tuple):
    """Represents a generator of the tensor star algebra - a tuple (possibly
    with zero components) of elements in the same algebra.

    """
    def __new__(cls, data, parent = None, idem = None):
        return tuple.__new__(cls, tuple(data))

    def __init__(self, data, parent = None, idem = None):
        """Specifies the tuple of generators, and the algebra."""
        # Note tuple initialization is automatic
        if parent == None:
            assert len(data) > 0
            parent = TensorStar(data[0].parent)
        assert all([factor.parent == parent.baseModule for factor in data])
        Generator.__init__(self, parent)
        if len(data) > 0:
            self.left_idem = data[0].getLeftIdem()
            self.right_idem = data[-1].getRightIdem()
        else:
            assert idem != None
            self.left_idem = self.right_idem = idem

    def slice(self, start, end = None):
        """Returns the generator of TensorStar that contains factors in the
        range [start, end) (same convention as python slicing).

        If end is omitted, slice to the end of sequence.

        """
        if end is None:
            end = len(self)
        assert start <= end
        if start == end == 0:
            new_left_idem = new_right_idem = self.left_idem
        elif start == end == len(self):
            new_left_idem = new_right_idem = self[start-1].getRightIdem()
        else:
            new_left_idem = self[start].getLeftIdem()
            new_right_idem = self[end-1].getRightIdem()
        return TensorStarGenerator(self[start:end], self.parent, new_left_idem)

    def getLeftIdem(self):
        return self.left_idem

    def getRightIdem(self):
        return self.right_idem

class TensorStarElement(Element):
    """Represents an element of the tensor star algebra."""
    def __init__(self, data = None, parent = None):
        """If the keys are tuples, convert them to tensor star generators."""
        if data is None:
            data = {}
        data_processed = {}
        for term, coeff in dict(data).items():
            if isinstance(term, TensorStarGenerator):
                data_processed[term] = coeff
            else:
                data_processed[TensorStarGenerator(term, parent)] = coeff
        Element.__init__(self, data_processed)

class TensorStar(FreeModule):
    """Represents a free module that is the direct sum of n'th tensor product,
    over n >= 0, of some free module A. So each generator is a (possibly empty)
    sequence of generators of A.

    """
    def __init__(self, baseModule):
        """Specifies the base module A. All generators are then tuples of
        generators in A.

        """
        self.baseModule = baseModule
        FreeModule.__init__(self, baseModule.ring)

class CobarAlgebra(TensorStar, DGAlgebra):
    """The tensor star module over a dg-algebra can be given a dg-algebra
    structure. Multiplication is by joining the two sequences, and differential
    is the dual of either taking differential of one term of the sequence, or
    multiplying together two adjacent terms of the sequence.

    """
    def __init__(self, baseAlgebra):
        """Specifies the base algebra A. """
        assert isinstance(baseAlgebra, DGAlgebra)
        TensorStar.__init__(self, baseAlgebra)
        DGAlgebra.__init__(self, baseAlgebra.ring)

    def _singleDiff(self, gen):
        """Compute the differential, in the cobar-algebra, of the generator
        ((gen)), where gen is a generator of the base algebra. Sum of the dual
        of differential and multiplication.

        antiDiff() and factor() must be implemented for generators of the base
        algebra.

        """
        assert gen.parent == self.baseModule
        result = E0
        for term, coeff in gen.antiDiff().items():
            result += coeff * TensorStarGenerator((term,), self)
        for ((a, b), coeff) in gen.factor().items():
            if not (a.isIdempotent() or b.isIdempotent()):
                result += coeff * TensorStarGenerator((a, b), self)
        return result

    def diff(self, gen):
        """Compute the differential, in the cobar-algebra, of any generator.
        This is defined as either taking the anti-differential or factoring one
        term of the sequence.

        """
        return E0.accumulate(
            [gen.slice(0,i) * self._singleDiff(gen[i]) * gen.slice(i+1)
             for i in range(len(gen))])

    def multiply(self, gen1, gen2):
        """Multiplication is joining two sequences. """
        if not isinstance(gen1, TensorStarGenerator):
            return NotImplemented
        if not isinstance(gen2, TensorStarGenerator):
            return NotImplemented
        assert gen1.parent == self and gen2.parent == self, \
            "Algebra not compatible."
        return 1*TensorStarGenerator(tuple(gen1)+tuple(gen2), self)

def simplifyComplex(arrows, default_coeff = 0, find_homology_basis = False):
    """Simplify complex using the cancellation lemma.

    ``arrows`` specify the complex to be simplified. It is a dictionary whose
    keys are generators, and values are dictionaries mapping generators to
    coefficients. So ``arrows[x][y] = coeff`` means there is an arrow from x
    to y with coefficient ``coeff``.

    The simplification is done through the cancellation lemma: for any arrow
    from x to y, with invertible coefficient ``coeff``, the generators x and y
    can be cancelled as follows: remove x, y, and all arrows entering or
    leaving these two generators. For each arrow from a to y with coefficient
    c1, and each arrow from x to b with coefficient c2 in the previous complex,
    add an arrow from a to b with coefficient c_1*(coeff^-1)*c_2. The new
    coefficient is added onto the previous one if an arrow already exists from
    a to b, possibly cancelling the previous arrow.

    ``default_coeff`` specifies the value to be used for zero coefficients. It
    should be 0 if the coefficients are numbers, and E0 if they are of type
    Element.

    ``find_homology_basis`` specifies whether to keep track of the identity of
    the generators. Since very often we are only interested in finding a
    homotopy equivalent chain complex or module, this is not automatically
    done. However, if it is necessary to find the generators of the homology
    of a chain complex, in terms of the generators of the original complex,
    this parameter should be set to True. Then an attribute ``prev_meaning`` is
    set for each generator of ``arrows`` in the returned result, expressing
    each generator as a linear combination of generators of the original
    complex.

    This implementation uses two optimizations. First, a rev_arrows dictionary
    is generated and updated along with arrows, keeping for each generator y,
    the list of arrows going into y. This speeds up the query for the list of
    arrows going into y. Second, for each arrow from x to y, we compute its
    degree (|A|-1)(|B|-1), where |A| is the set of arrows going into y, and |B|
    is the set of arrows coming from x. This equals the number of arrows added
    in cancelling the arrow from x to y. We try to cancel arrows in increasing
    order of degree, using a priority queue (but not strictly so, as the degree
    of an arrow may have changed between when it is added to queue and when it
    is used for cancellation.

    """
    # Produce rev_arrows
    rev_arrows = dict()
    for x in arrows:
        rev_arrows[x] = dict()
    for x in arrows:
        for y in arrows[x]:
            if isinstance(arrows[x][y], Element):
                rev_arrows[y][x] = arrows[x][y].copy()
            else:
                rev_arrows[y][x] = arrows[x][y]

    cancel_list = []
    def tryAddEdge(x, y):
        """If the arrow from x to y is cancellable, then add it to cancel_list
        (a heap), along with its degree.

        """
        coeff = arrows[x][y]
        if coeff.invertible():
            cur_degree = (len(arrows[x])-1)*(len(rev_arrows[y])-1)
            heapq.heappush(cancel_list, (cur_degree, x, y))

    # Make a initial list of cancellable arrows
    for x in arrows:
        for y in arrows[x]:
            tryAddEdge(x, y)

    # Initialize prev_meaning attribute
    if find_homology_basis:
        for x in arrows:
            x.prev_meaning = 1*x

    def cancelEdge(x, y):
        """Cancel the edge from x to y."""
        coeff = arrows[x][y]
        assert coeff.invertible()
        inv_coeff = coeff.inverse()

        # List of edges going into y (other than that from x and y)
        alist = [(term, coeff) for term, coeff in rev_arrows[y].items()
                 if term not in (x, y)]
        # List of edges coming from x (other than that going to x and y)
        blist = [(term, coeff) for term, coeff in arrows[x].items()
                 if term not in (x, y)]

        # Remove all edges going into x or y
        for term in arrows[x]:
            if term not in (x, y):
                del rev_arrows[term][x]
        for term in arrows[y]:
            if term not in (x, y):
                del rev_arrows[term][y]
        for term in rev_arrows[x]:
            if term not in (x, y):
                del arrows[term][x]
        for term in rev_arrows[y]:
            if term not in (x, y):
                del arrows[term][y]
        # Remove x and y
        del arrows[x], arrows[y], rev_arrows[x], rev_arrows[y]

        # Add arrows from alist to blist
        for a, c1 in alist:
            c1_invc = c1 * inv_coeff
            for b, c2 in blist:
                new_coeff = c1_invc * c2
                if b not in arrows[a]:
                    arrows[a][b] = default_coeff
                    rev_arrows[b][a] = default_coeff
                arrows[a][b] += new_coeff
                rev_arrows[b][a] += new_coeff
                if arrows[a][b] == 0:
                    del arrows[a][b], rev_arrows[b][a]
                else:
                    tryAddEdge(a, b)

        # Update prev_meaning
        if find_homology_basis:
            for a, c1 in alist:
                a.prev_meaning += (c1 * inv_coeff) * x.prev_meaning

    # Main loop: try to cancel each edge in the queue
    while cancel_list:
        degree, x, y = heapq.heappop(cancel_list)
        if x in arrows and y in arrows[x] and arrows[x][y].invertible():
            cancelEdge(x, y)

    return arrows
