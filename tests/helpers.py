"""
Helpers shared by several test modules.

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

import operator
from functools import reduce

from sage.matrix.constructor import matrix
from sage.matroids.basis_matroid import BasisMatroid
from sage.rings.integer_ring import ZZ

from arithmat.arithmetic_matroid import (
    ArithmeticMatroid,
    BasisArithmeticMatroid,
    LinearArithmeticMatroid,
    ToricArithmeticMatroid,
)

def representation_to_matroid(A):
    """
    Given an integer matrix A, return the associated ArithmeticMatroid.
    """
    E = range(A.ncols())

    def rk(X):
        return A[:, list(X)].rank()

    def m(X):
        return reduce(operator.mul, [d for d in A[:, list(X)].elementary_divisors() if d != 0], 1)

    return ArithmeticMatroid(E, rk, m)

# 0 and 2 are loops with torsion: m(empty) = 2, m({2}) = 1
T_A = matrix(ZZ, [[1, 0, 1, 1], [0, 1, 0, 1]], immutable=True)
T_Q = matrix(ZZ, [[2], [0]], immutable=True)
# 2 is a coloop with m({2}) = 2
A_LIN = matrix(ZZ, [[1, 1, 0], [0, 0, 2]], immutable=True)
def M_LIN(X):
    return reduce(operator.mul, [d for d in A_LIN[:, sorted(X)].elementary_divisors() if d != 0], 1) if X else 1
# U(2,5) on {1, ..., 5}, with m = 2 on pairs inside {3, 4, 5}
def RK_U25(X):
    return min(2, len(X))

def M_U25(X):
    return 2 if len(X) == 2 and all(x in [3, 4, 5] for x in X) else 1

# 2 is a loop without torsion
A_FREE_LOOP = matrix(ZZ, [[1, 0, 0], [0, 1, 0]], immutable=True)
# 2 is a coloop with m({2}) = 1
A_PLAIN_COLOOP = matrix(ZZ, [[1, 1, 0], [0, 0, 1]], immutable=True)

def make_matroid(name, T=None):
    # matroids built from the same T share its rank and multiplicity functions
    if T is None:
        T = ToricArithmeticMatroid(T_A, T_Q)
    return {
        "rank": lambda: ArithmeticMatroid(T.groundset(), T._rank, T._multiplicity),
        "toric": lambda: T,
        "basis": lambda: BasisArithmeticMatroid(BasisMatroid(groundset=T.groundset(), bases=T.bases()),
                                        multiplicity_function=T._multiplicity),
        "U(2,5)": lambda: ArithmeticMatroid([1, 2, 3, 4, 5], RK_U25, M_U25),
        "linear": lambda: LinearArithmeticMatroid(A_LIN, multiplicity_function=M_LIN),
        "free-loop": lambda: ToricArithmeticMatroid(A_FREE_LOOP),
        "plain-coloop": lambda: ToricArithmeticMatroid(A_PLAIN_COLOOP),
    }[name]()

TORIC_NAMES = ["toric", "free-loop", "plain-coloop"]
MATROID_NAMES = ["rank", "U(2,5)", "basis", "linear"] + TORIC_NAMES
