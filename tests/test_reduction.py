"""
Tests for reduction.

Copyright (C) 2020 Roberto Pagaria
Copyright (C) 2020 Giovanni Paolini

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program. If not, see http://www.gnu.org/licenses/.
"""

from sage.combinat.subset import powerset
from sage.matrix.constructor import matrix
from sage.rings.integer_ring import ZZ

from arithmat.arithmetic_matroid import ArithmeticMatroid
from helpers import representation_to_matroid

def test_example():
    E = [1, 2, 3, 4, 5, 6]

    def rk(X):
        return min(2, len(X))

    def m(X):
        if len(X) == 0:
            return 1
        elif len(X) == 1:
            if X == frozenset([1]) or X == frozenset([2]):
                return 2
            else:
                return 1
        elif sum(1 for i in X if i not in [1, 2]) >= 2:
            return 1
        elif len(X) == 2:
            if X == frozenset([1, 2]):
                return 4
            else:
                return 2
        elif 3 in X:
            # X = {1,2,3}
            assert X == frozenset([1, 2, 3])
            return 1
        else:
            # X = {1,2,j} with j > 3
            return 2

    M = ArithmeticMatroid(E, rk, m)
    assert M.is_valid()

    M1 = M.reduction()
    assert M1.multiplicity([1, 2, 3]) == 2

    for X in powerset(E):
        if tuple(sorted(X)) != (1, 2, 3):
            assert M1.multiplicity(X) == M.multiplicity(X)

    assert not M1.is_valid()

    assert M1.is_torsion_free()
    assert M1.is_surjective()

def test_dual():
    A = matrix(ZZ, [[-1, 1, 0, -1, 2, 7], [6, 1, -1, -2, 2, 5]])
    M = representation_to_matroid(A)
    M2 = ArithmeticMatroid(M.groundset(), M.rank, lambda X: M._multiplicity(X) ** 2)

    assert M.is_valid()
    assert M2.is_valid()
    assert not M2.is_representable()

    assert not M.dual().reduction().is_isomorphic(M.reduction())
    assert M.dual().reduction().is_isomorphic(M.reduction().dual())
    assert M2.dual().reduction().is_isomorphic(M2.reduction().dual())
