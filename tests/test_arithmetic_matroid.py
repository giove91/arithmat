"""
Tests for ArithmeticMatroid.

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

import copy

from sage.arith.misc import euler_phi, gcd
from sage.matrix.constructor import Matrix, matrix
from sage.matrix.special import random_matrix
from sage.matroids.linear_matroid import LinearMatroid
from sage.rings.integer_ring import ZZ
from sage.rings.rational_field import QQ

from arithmat.arithmetic_matroid import (
    ArithmeticMatroid,
    LinearArithmeticMatroid,
    ToricArithmeticMatroid,
)
from arithmat.shnf import signed_hermite_normal_form
from helpers import representation_to_matroid


def test_trivial():
    E = range(1, 10)

    def rk(X):
        return len(X)

    def m(X):
        return 1

    M = ArithmeticMatroid(E, rk, m)
    assert M.is_valid()
    assert M.full_rank() == 9
    assert M.is_representable()

def test_valid():
    # matroid with representation
    # [ 1, a, 1 ]
    # [ 0, b, 5 ]
    # with a|b and (b,5)=1

    E = [1, 2, 3]

    def rk(X):
        return min(2, len(X))

    def m(X):
        a = 3
        b = 6

        X = tuple(sorted(X))
        return {
            (): 1,
            (1,): 1,
            (2,): a,
            (3,): 1,
            (1, 2): b,
            (1, 3): 5,
            (2, 3): 5 * a - b,
            (1, 2, 3): 1,
        }[X]

    M = ArithmeticMatroid(E, rk, m)
    assert M.is_valid()
    assert M.is_representable()
    assert M.full_multiplicity() == 1

def test_not_valid():
    E = [1, 2, 3]

    def rk(X):
        return min(2, len(X))

    def m(X):
        a = 3
        b = 6

        X = tuple(sorted(X))
        return {
            (): 1,
            (1,): 1,
            (2,): a,
            (3,): 1,
            (1, 2): b,
            (1, 3): 5,
            (2, 3): 5 * a - b + 1,
            (1, 2, 3): 1,
        }[X]

    M = ArithmeticMatroid(E, rk, m)
    assert not M.is_valid()

def test_valid2():
    # valid, not representable, not orientable
    E = [1, 2, 3, 4, 5]

    def rk(X):
        return min(2, len(X))

    def m(X):
        if len(X) == 2 and all(x in [3, 4, 5] for x in X):
            return 2
        else:
            return 1

    M = ArithmeticMatroid(E, rk, m)
    assert M.is_valid()
    assert not M.is_representable()
    assert not M.is_orientable()
    assert M.representation() == None
    assert M.full_multiplicity() == 1

def test_realizable():
    # representable with a 2x2 matrix
    for a in range(6):
        for b in range(1, 6):
            A = matrix(ZZ, [[1, a], [0, b]])
            M = representation_to_matroid(A)

            assert M.is_valid()
            assert M.is_representable()
            assert M.is_orientable()

            R = ZZ['x, y']
            x, y = R._first_ngens(2)
            assert M.arithmetic_tutte_polynomial() == x ** 2 + x * (gcd(a, b) - 1) + b - gcd(a, b)

def test_representable2():
    A = matrix(ZZ, [[-1, -29, -1, 1], [1, -1, 0, 1]])
    M = representation_to_matroid(A)

    assert M.is_valid()
    assert M.is_representable()
    assert M.is_orientable()

def test_representable3():
    A = matrix(ZZ, [[-1, 1, 0, 0, -1], [6, 1, 1, -1, -1]])
    M = representation_to_matroid(A)

    assert M.is_valid()
    assert M.is_representable()
    assert M.is_orientable()

def test_representable4():
    A = matrix(ZZ, [[2, 2, 1, 0, 0], [1, 5, -1, 1, -2], [-2, 1, 0, -1, -1]])
    M = representation_to_matroid(A)

    assert M.is_valid()
    assert M.is_representable()
    assert M.is_orientable()

def test_representable_random():
    A = random_matrix(ZZ, 4, 6)
    M = representation_to_matroid(A)

    assert M.is_valid()
    assert M.is_representable()
    assert M.is_orientable()

def test_pseudo():
    E = [1, 2, 3]

    def rk(X):
        return min(2, len(X))

    def m(X):
        if len(X) == 2:
            return 3
        else:
            return 1

    M = ArithmeticMatroid(E, rk, m)
    assert M.is_valid()
    assert not M.is_representable()
    assert M.is_orientable()

def test_non_representable():
    A = matrix(ZZ, [[-1, 1, 0, -1], [6, 1, -1, -2]])
    M = representation_to_matroid(A)
    M2 = ArithmeticMatroid(M.groundset(), M.rank, lambda X: M._multiplicity(X) ** 2)

    assert M2.is_valid()
    assert M.is_representable()
    assert not M2.is_representable()

def test_non_representable2():
    A = matrix(ZZ, [[-1, 1, 0, -1, 2, 7], [6, 1, -1, -2, 2, 5]])
    M = representation_to_matroid(A)
    M2 = ArithmeticMatroid(M.groundset(), M.rank, lambda X: M._multiplicity(X) ** 2)

    assert M2.is_valid()
    assert M.is_representable()
    assert not M2.is_representable()

def test_non_orientable():
    E = range(1, 7)

    def rk(A):
        return min(len(A), 2)

    def m(X):
        A = tuple(sorted(X))
        if A == (1, 2, 3):
            return 1
        if len(A) == 0:
            return 1
        elif A == (1,) or A == (2,):
            return 2
        elif len(A) == 1:
            return 1
        elif len(A) == 2 and len([x for x in A if x <= 2]) == 1:
            return 2
        elif 1 in A and 2 in A and len(A) == 3:
            return 2
        elif A == (1, 2):
            return 4
        else:
            return 1

    M = ArithmeticMatroid(E, rk, m)
    assert M.is_valid()
    assert not M.is_representable()
    assert not M.is_orientable()

def test_non_essential():
    A = matrix(ZZ, [[-1, 1, 0], [6, 1, -1], [2, -3, 0], [1, 2, 3], [-1, 0, 0]])
    M = representation_to_matroid(A)

    assert M.is_valid()
    assert M.is_representable()

def test_hash():
    E = [1, 2, 3, 4, 5]

    def rk(X):
        return min(2, len(X))

    def m(X):
        if len(X) == 2 and all(x in [3, 4, 5] for x in X):
            return 2
        else:
            return 1

    M = ArithmeticMatroid(E, rk, m)
    s = {M}
    assert hash(M) == hash((frozenset(E), 2, 1, 1))
    assert M.groundset() == frozenset(E)

def test_isomorphism():
    E = [1, 2, 3, 4, 5]

    def rk(X):
        return min(2, len(X))

    def m(X):
        if len(X) == 2 and all(x in [3, 4, 5] for x in X):
            return 2
        else:
            return 1

    M = ArithmeticMatroid(E, rk, m)
    assert M.is_valid()
    assert not M.is_representable()
    assert not M.is_orientable()

    E1 = [4, 5, 6, 7, 8]

    def m1(X):
        if len(X) == 2 and all(x in [6, 7, 8] for x in X):
            return 2
        else:
            return 1

    M1 = ArithmeticMatroid(E1, rk, m1)
    assert M1.is_valid()
    assert not M1.is_representable()
    assert not M1.is_orientable()

    assert M.is_isomorphism(M1, {i: i + 3 for i in E})
    assert not M.equals(M1)

    M2 = ArithmeticMatroid(copy.copy(E), rk, m)
    assert M.equals(M2)
    assert M2.equals(M)

def test_isomorphic():
    A = matrix(ZZ, [[1]])
    B = matrix(ZZ, [[3, 1], [0, 1]])
    M = ToricArithmeticMatroid(A)
    DN = ToricArithmeticMatroid(B).minor(deletions=[0])
    CN = ToricArithmeticMatroid(B).minor(contractions=[0])
    assert DN.is_isomorphic(M)
    assert not CN.is_isomorphic(M)

def test_num_representations():
    r = 3

    for m in range(23, 26):  # m(E)
        E = range(r)
        A = matrix(ZZ, r, r, lambda i, j: 1 if i == j and i < r - 1 else 0 if j < r - 1 else m if i == r - 1 else 1)

        M = representation_to_matroid(A)
        assert M.num_representations() == (euler_phi(m) // 2) ** (r - 1)

def test_gcd():
    M = ToricArithmeticMatroid(matrix(ZZ, [[6, 3, -2, 2], [3, 21, 0, -9], [-1, -4, 3, -2]]))

    assert M.is_torsion_free()
    assert not M.is_surjective()
    assert M.is_representable()
    assert M.is_gcd()
    assert not M.is_strong_gcd()

    M1 = M.reduction()

    assert M1.is_torsion_free()
    assert M1.is_surjective()
    assert M1.is_representable()
    assert M1.is_gcd()

    assert M1.is_strong_gcd()


def test_linear_representation_is_sage():
    # arithmetic representation() hid LinearMatroid.representation(), so Sage's own calls to it with its arguments (Matrix(M), representation_vectors()) raised TypeError
    A = matrix(ZZ, [[1, 0, 1, 1], [0, 1, 1, 2]])
    L = LinearArithmeticMatroid(A, multiplicity_function=lambda X: representation_to_matroid(A).multiplicity(X))
    plain = LinearMatroid(A)

    assert L.representation() == plain.representation()
    assert L.representation() == A
    assert Matrix(L) == A
    assert L.representation_vectors() == plain.representation_vectors()
    B = L.basis()
    assert L.representation(B) == plain.representation(B)
    assert L.representation(B, reduced=True) == plain.representation(B, reduced=True)
    assert L.representation(labels=True) == plain.representation(labels=True)

    # Sage returns a new mutable matrix each call, and changing it leaves L as it was
    R = L.representation()
    assert R.is_mutable()
    assert L.representation() is not R
    R[0, 0] = 5
    assert L.representation() == A

    # the arithmetic search keeps its own name
    assert L.is_representable()
    H = L.arithmetic_representation()
    assert H == signed_hermite_normal_form(A)
    assert L.check_representation(H)
    assert L.arithmetic_representation(shnf=False).dimensions() == A.dimensions()

def test_linear_not_arithmetically_representable():
    # the stored matrix has a 2x2 minor of determinant 2, but every multiplicity is 1
    A = matrix(QQ, [[1, 0, 1, 1], [0, 1, 1, 2]])
    L = LinearArithmeticMatroid(A, multiplicity_function=lambda X: 1)
    assert L.representation() == A  # Sage's representation() is still the stored matrix
    assert L.arithmetic_representation() is None
    assert not L.is_representable()

def test_representation_alias():
    # outside LinearArithmeticMatroid, representation() is the arithmetic search
    E = [1, 2, 3, 4, 5]
    M = ArithmeticMatroid(E, lambda X: min(2, len(X)), lambda X: 1)
    assert M.representation() == M.arithmetic_representation()
    assert (M.representation(shnf=False) is None) == (M.arithmetic_representation(shnf=False) is None)
