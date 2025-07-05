# Bordered Heegaard Floer Homology

This Python module does computations in Bordered Heegaard Floer
Homology.  To install, do
```
python -m pip install https://github.com/NathanDunfield/bfh_python/archive/master.zip
```
and then start Python and do
```
>>> from bfh_python.braid import BridgePresentation
>>> name = 'P(-2,3,5)'
>>> left_closure = (6,3,2,5,4,1)
>>> braid = [-1,-1,3,3,3,5,5,5,5,5]
>>> right_closure = (6,3,2,5,4,1)
>>> pretzel = BridgePresentation(name, left_closure, braid, right_closure)
>>> pretzel.getHF()
(compose Mor 28) 1,2,3,4,7,8,11,12,15,14,10,13,16,15,11,14,17,16,14,17,20,19,17,20,23,22,20,23,26
 Chain complex.
d(g1) = 0

```
For more on how to use `bfh_python`, see [the documentation by Robert
Lipshitz](https://pages.uoregon.edu/lipshitz/bfh/).

# Test suite

To makes sure everything is working, after installing `bfh_python` you
can do
```
python -m bfh_python.test
```
Please note that the testing of the final `braids` module can take
more than 20 minutes.

# Performance

The module `bfh_python` is pure Python and is largely self-contained
with no dependencies outside a few modules in Python's standard
library.  As such it can be run with [PyPy](https://pypy.org/), an
alternate Python implementation with a Just-in-Time (JIT) compiler.
On the test suite for `bfh_python`, PyPy bested Python 3.13 by a
factor of 3.5 under macOS with an Intel (x86) processor.


# Underlying theory

Based on a previous program (written in Sage) by Robert Lipshitz,
Peter Ozsvath and Dylan Thurston, which is [available on the former's
webpage](https://pages.uoregon.edu/lipshitz/Code/).

Implements various algorithms described in the following papers of
R. Lipshitz, P. Ozsvath, and D. Thurston

- arXiv:math/0810.0687 Bordered Heegaard Floer homology: Invariance and pairing
- arXiv:math/1003.0598 Bimodules in bordered Heegaard Floer homology
- arXiv:math/1005.1248 Heegaard Floer homology as morphism spaces
- arXiv:math/1010.2550 Computing HF^ by factoring mapping classes
- arXiv:math/1011.0499 Bordered Floer homology and the branched double cover I

and by B. Zhan

- arXiv:math/1403.6215 Explicit Koszul-dualizing bimodules in bordered Heegaard Floer homology
- arXiv:math/1405.3004 Combinatorial Proofs in Bordered Heegaard Floer homology

as well as other (not yet published) work.

# License

Copyright (C) 2013 by Bohua Zhan.

This program is free software; you can redistribute it and/or modify it under
the terms of the GNU General Public License as published by the Free Software
Foundation; either version 2 of the License, or (at your option) any later
version.
