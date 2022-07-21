"""Various utilities useful for the project."""

from math import gcd
from numbers import Number

def memorize(function):
    """Function decorator: memorize returned values of this function.
    Based on: Daniel Lawrence, A simple example using a python cache decorator.

    """
    memo = {}
    class NoneSymbol(object):
        pass

    def wrapper(*args, **kwargs):
        key = args + tuple((k, v) for k, v in list(kwargs.items()))
        # Use memorization. Don't keep NotImplemented values.
        # Will throw TypeError if key is not hashable.
        val = memo.get(key, NoneSymbol)
        if val is not NoneSymbol:
            return val
        else:
            rv = function(*args, **kwargs)
            if rv is not NotImplemented:
                memo[key] = rv
            return rv
    return wrapper

def memorizeHash(function):
    """Decorator for hash function. Memorize hash value within the object
    attribute _hash_val.

    """
    def wrapper(*args):
        _self = args[0]
        if hasattr(_self, '_hash_val'):
            return _self._hash_val
        else:
            rv = function(*args)
            _self._hash_val = rv
            return rv
    return wrapper

def trace(function):
    """Decorator: print input and ouput for every invocation of the function."""
    def wrapper(*args):
        print("\nInputs:", args)
        rv = function(*args)
        print("Output:", rv)
        return rv
    return wrapper

def tolist(obj):
    """Force obj into a list."""
    if isinstance(obj, list):
        return obj
    else:
        return [obj]

def flatten(lst):
    """Flatten a once-nested list."""
    return sum(lst, [])

def fracToInt(frac):
    """Convert frac (which must have integer value regardless of type)
    to integers.

    """
    if isinstance(frac, int):
        return frac
    assert frac.denominator == 1
    return frac.numerator

def find(lst, item):
    """Find the location of the first occurrence of item in lst. Returns -1 if
    item does not exist in lst.

    """
    for i in range(len(lst)):
        if lst[i] == item:
            return i
    return -1

def sumColumns(matrix, num_col):
    """matrix is a list of lists of length num_col. Sum up by column."""
    result = [0]*num_col
    for row in matrix:
        assert len(row) == num_col
        for i in range(num_col):
            result[i] += row[i]
    return result

def subset(elts):
    """Returns the list of subsets of the given set (each element is a tuple).

    """
    if len(elts) == 0:
        return [()]
    all_but_zero = subset(elts[1:])
    return all_but_zero + [(elts[0],) + s for s in all_but_zero]

def _dictAddTo(dict1, dict2):
    """Add dict2 onto dict1 in place. If dict2 is a list, add each element of
    dict2 in place.

    """
    dict2 = [curdict.copy() for curdict in tolist(dict2) if curdict != 0]
    if dict1 == 0:
        if len(dict2) == 0:
            return dict1
        else:
            dict1 = dict2[0]
            dict2 = dict2[1:]
    for curdict in dict2:
        assert type(dict1) == type(curdict), "Incompatible types: %s, %s" % \
            (str(type(dict1)), str(type(curdict)))
        for k, v in list(curdict.items()):
            if k in dict1:
                dict1[k] += v
                if dict1[k] == 0:
                    del dict1[k]
            else:
                dict1[k] = v
    return dict1

def _dictMult(dict1, scalar):
    """Return a new dictionary with same type as self, the same keys, and
    each value multiplied by scalar.

    """
    if not isinstance(scalar, Number):
        return NotImplemented

    result = type(dict1)((k, scalar * v) for k, v in list(dict1.items()) if
                             scalar * v != 0)
    return result

class SummableDict(dict):
    """A dictionary type that supports sums and multiplication by scalar. Works
    in the same way as a free module with generators as keys and coefficients
    as values. Zeroes are automatically thrown away.

    """
    def __add__(self, other):
        return _dictAddTo(type(self)(), [self, other])

    def __iadd__(self, other):
        return _dictAddTo(self, other)

    def __sub__(self, other):
        return _dictAddTo(type(self)(), [self, -1*other])

    def __isub__(self, other):
        return _dictAddTo(self, -1*other)

    def accumulate(self, lst):
        """Similar to +=, except returns the sum."""
        return _dictAddTo(self, lst)

    def __mul__(self, other):
        return _dictMult(self, other)

    def __rmul__(self, other):
        return _dictMult(self, other)

    def __eq__(self, other):
        """Comparison to 0 is a test for empty dictionary."""
        if isinstance(other, int) and other == 0:
            return len(self) == 0
        else:
            return dict.__eq__(self, other)

    def __ne__(self, other):
        return not (self == other)

    def copy(self):
        """Copy function should preserve type."""
        return type(self)(self)

    def translateKey(self, key_map):
        """Translate keys of this dictionary using the dictionary key_map. All
        keys in self must appear in key_map and will be replaced by the
        corresponding value.

        """
        return type(self)([(key_map[k], v) for k, v in list(self.items())])

    def getElt(self):
        """Returns an arbitrary key from this dictionary. Must be non-empty."""
        return next(iter(self))

class Ring(object):
    def convert(self, data):
        """Try to convert data to an element of this ring."""
        raise NotImplementedError("convert function not specified for ring.")

class RingElement(Number):
    pass

class ModNRing(Ring):
    """The ring Z/nZ."""
    def __init__(self, n):
        self.n = n
        self.zero = self.convert(0)
        self.one = self.convert(1)

    def add(self, elt1, elt2):
        elt1, elt2 = self.convert(elt1), self.convert(elt2)
        if elt1 is NotImplemented or elt2 is NotImplemented:
            return NotImplemented
        return ModNElement(self, (elt1.val+elt2.val)%self.n)

    def multiply(self, elt1, elt2):
        elt1, elt2 = self.convert(elt1), self.convert(elt2)
        if elt1 is NotImplemented or elt2 is NotImplemented:
            return NotImplemented
        return ModNElement(self, (elt1.val*elt2.val)%self.n)

    def __eq__(self, other):
        return self.n == other.n

    def __ne__(self, other):
        return not (self == other)

    def __hash__(self):
        return hash((self.n, "ModNRing"))

    def convert(self, data):
        """Try to convert data to an element of this ring."""
        if isinstance(data, ModNElement) and data.parent == self:
            return data
        if isinstance(data, int):
            return ModNElement(self, data % self.n)
        return NotImplemented

class ModNElement(RingElement):
    """An element in a ring Z/nZ."""
    def __init__(self, parent, val):
        self.parent = parent
        self.val = val

    def __str__(self):
        return str(self.val)

    def __repr__(self):
        return str(self.val)

    def __add__(self, other):
        return self.parent.add(self, other)

    def __radd__(self, other):
        return self.parent.add(self, other)

    def __mul__(self, other):
        return self.parent.multiply(self, other)

    def __rmul__(self, other):
        return self.parent.multiply(self, other)

    def __eq__(self, other):
        """Can compare to integer 0 or 1."""
        if isinstance(other, int) and (other == 0 or other == 1):
            return self.val == other
        else:
            return self.val == other.val

    def invertible(self):
        """Returns whether this element is invertible in the ring."""
        return gcd(self.val, self.parent.n) == 1

    def inverse(self):
        """Returns the inverse of this element. Must be invertible"""
        # Currently only implemented for n = 2
        assert self.parent.n == 2
        return self

class Integer(Ring):
    """The ring Z."""
    def convert(self, data):
        """Try to convert data to an element of this ring."""
        assert isinstance(data, int)
        return data

class IntegerElement(RingElement, int):
    """An element in a ring Z."""
    def __new__(cls, parent, val):
        return int.__new__(cls, val)

    def __init__(self, parent, val):
        self.parent = parent

class NamedObject(object):
    """Provides functionality for an object to be described by name. If this is
    listed as a parent class, an object will use name in equality comparisons,
    hash functions, and string outputs.

    """
    def __init__(self, name):
        self.name = name

    def __eq__(self, other):
        return self.name == other.name

    def __ne__(self, other):
        return not (self == other)

    def __lt__(self, other):
        return self.name < other.name

    def __le__(self,other):
        return self.name <= other.name

    def __gt__(self, other):
        return self.name > other.name

    def __ge__(self, other):
        return self.name >= other.name
    
    @memorizeHash
    def __hash__(self):
        return hash(self.name)

    def __str__(self):
        return str(self.name)

    def __repr__(self):
        return str(self.name)

class MorObject(object):
    """If this is list as a parent class, an object be treated as a generator
    of some morphism complex. It will have source, coeff, and target fields.
    These will be used for equality comparisons, hash functions, and string
    outputs.

    """
    def __init__(self, source, coeff, target):
        self.source = source
        self.coeff = coeff
        self.target = target

    def __eq__(self, other):
        return self.source == other.source and self.coeff == other.coeff \
            and self.target == other.target

    def __ne__(self, other):
        return not (self == other)

    @memorizeHash
    def __hash__(self):
        return hash((self.source, self.coeff, self.target))

    def __str__(self):
        return "%s->%s*%s" % \
            (str(self.source), str(self.coeff), str(self.target))

    def __repr__(self):
        return str(self)

def safeMultiply(a, b):
    """Safely multiply the two sides using __mul__ and __rmul__. Return
    NotImplemented if both fails.

    """
    try:
        prod = a.__mul__(b)
    except TypeError:
        prod = NotImplemented
    if prod is NotImplemented:
        try:
            prod = b.__rmul__(a)
        except TypeError:
            prod = NotImplemented
    return prod

# Most commonly used rings: Z/2Z and Z
F2 = ModNRing(2)
ZZ = Integer()

# Constants for positive and negative orientation
POS, NEG = 1, -1

# Constants for left and right action
ACTION_LEFT, ACTION_RIGHT = 0, 1
def sideStr(side):
    if side == ACTION_LEFT: return "LEFT"
    else: return "RIGHT"
def oppSide(side):
    if side == ACTION_LEFT: return ACTION_RIGHT
    else: return ACTION_LEFT

# Constant for controlling amount of printing
PRINT_PROGRESS = 1

# How much assertions do you want to do? Higher value means more checks and
# slower program
ASSERT_LEVEL = 0

# For each grading group and grading set object, its type attribute is one of
# these
BIG_GRADING, SMALL_GRADING = 0, 1
DEFAULT_GRADING = SMALL_GRADING
def grTypeStr(gr_type):
    if gr_type == BIG_GRADING: return "big"
    else: return "small"

# Whether to use multiplicity-one algebra (or the full algebra)
MULT_ONE = True
