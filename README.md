# Arithmat
Sage implementation of arithmetic matroids and toric arrangements

Authors: Giovanni Paolini and Roberto Pagaria

## Quick start
...

## Overview

Arithmat is a Sage package that implements arithmetic matroids.
At its core there is the class `ArithmeticMatroidMixin`, which is intended to be used in combination with any existing `Matroid` subclass of Sage (e.g. `RankMatroid`, `BasisMatroid`, `LinearMatroid`) via multiple inheritance.
The most common combinations are already defined, for example: `ArithmeticMatroid` (deriving from `RankMatroid`), `BasisArithmeticMatroid` (deriving from `BasisMatroid`), `LinearArithmeticMatroid` (deriving from `LinearMatroid`).
An additional class `ToricArithmeticMatroid` is implemented, for arithmetic matroids constructed from a fixed representation.

## Documentation

### Import

All defined classes can be imported at once:
```sage
from arithmat import *
```
Alternatively, it is possible to import only specific classes:
```sage
from arithmat import ArithmeticMatroid, ToricArithmeticMatroid
```


### Available classes for arithmetic matroids

All classes for arithmetic matroids derive from `ArithmeticMatroidMixin` and from some subclass of Sage's `Matroid`.
The class `ArithmeticMatroidMixin` is not intended to be used by itself, but it is possible to subclass it in order to create new classes for arithmetic matroids (see below).

A general way to construct an instance of some arithmetic matroid class `XxxArithmeticMatroid` (apart from `ToricArithmeticMatroid`, which is special) is the following. Suppose that `XxxArithmeticMatroid` derives from `XxxMatroid`. Then an instance of `XxxArithmeticMatroid` can be constructed with `XxxArithmeticMatroid(..., multiplicity_function=m)`, where the dots stand for arguments to construct an instance of `XxxMatroid`, and `m` is the multiplicity function.

The classes which are already provided in `arithmat` are the following.

* `ArithmeticMatroid` (derives from `ArithmeticMatroidMixin` and `RankMatroid`).



### Available methods

All classes for arithmetic matroids must also derive from some subclass of Sage's `Matroid`.
In particular, all `Matroid` methods are still available. For example:
```sage
M = ToricArithmeticMatroid(matrix(ZZ, [[1,2,3], [0,1, 1]]))
list(M.bases())
# [frozenset([0, 1]), frozenset([0, 2]), frozenset([1, 2])]
```

All subclasses of `ArithmeticMatroidMixin` also implement the following methods.

* `multiplicity(X)`
  Return the multiplicity of the set `X`.


### Creating new classes for arithmetic matroids
