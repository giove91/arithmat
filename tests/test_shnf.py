"""
Tests for the signed Hermite normal form.

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

import itertools

import pytest
from sage.matrix.constructor import matrix
from sage.matrix.special import diagonal_matrix, random_matrix
from sage.misc.prandom import choice
from sage.rings.integer_ring import ZZ

from arithmat.shnf import signed_hermite_normal_form as shnf


@pytest.mark.parametrize("x", range(-5, 5))
def test_1x1(x):
    A = matrix(ZZ, [[x]])
    assert shnf(A) == matrix(ZZ, [[abs(x)]])


def test_2x2():
    A = matrix(ZZ, [[5, 8], [0, 3]])
    assert shnf(A) == matrix(ZZ, [[5, 1], [0, 3]])


@pytest.mark.parametrize("r, n, all_signs", [
    pytest.param(1, 2, True, id="1x2"),
    pytest.param(1, 7, False, id="1x7"),
    pytest.param(2, 2, True, id="2x2"),
    pytest.param(2, 3, True, id="2x3"),
    pytest.param(2, 4, True, id="2x4"),
    pytest.param(3, 4, True, id="3x4"),
    pytest.param(4, 6, False, id="4x6"),
    pytest.param(6, 10, False, id="6x10"),
])
def test_random(r, n, all_signs):
    # shnf is invariant under row operations and column sign changes;
    # all_signs tries every sign pattern, otherwise one random pattern
    A = random_matrix(ZZ, r, n)
    U = random_matrix(ZZ, r, r, algorithm='unimodular')
    if all_signs:
        sign_patterns = itertools.product([1, -1], repeat=n)
    else:
        sign_patterns = [[choice([1, -1]) for i in range(n)]]
    for diag in sign_patterns:
        S = diagonal_matrix(diag)
        assert shnf(A) == shnf(U * A * S)


@pytest.mark.parametrize("rows", [
    pytest.param([[1, 1, 1], [0, 2, 3], [0, 0, 6]], id="2-3-6"),
    pytest.param([[1, 1, 1], [0, 2, 4], [0, 0, 8]], id="2-4-8"),
])
def test_cardinality_four_orbit(rows):
    # in these examples, an orbit of cardinality 4 occurs
    # see [PP19, Example 6.8]
    A = matrix(ZZ, rows)
    assert shnf(A) == A

    for i in range(3):
        U = random_matrix(ZZ, 3, 3, algorithm='unimodular')
        for diag in itertools.product([1, -1], repeat=3):
            S = diagonal_matrix(diag)
            assert shnf(U * A * S) == A


@pytest.mark.parametrize("u", [1, 2, 4, 5])
def test_cardinality_four_orbit_representative(u):
    # every matrix in the orbit of [PP19, Example 6.8] has the same shnf
    A = matrix(ZZ, [[1, 1, u], [0, 2, 3], [0, 0, 6]])
    assert shnf(A) == matrix(ZZ, [[1, 1, 1], [0, 2, 3], [0, 0, 6]])


def test_is_in_hermite_form():
    # check that the output of shnf is in Hermite normal form
    for i in range(10):
        A = random_matrix(ZZ, 3, 4)
        B = shnf(A)
        assert B == B.echelon_form()


def test_zero_pivot():
    # in this example, the second column has a zero pivot
    A = matrix(ZZ, [[1, -1, -1, 11], [1, -1, 0, -1], [3, -3, 7, -2]])
    U = random_matrix(ZZ, 3, 3, algorithm='unimodular')
    for diag in itertools.product([1, -1], repeat=4):
        S = diagonal_matrix(diag)
        assert shnf(U * A * S) == shnf(A)
        assert shnf(U * A) == shnf(A)


def test_not_reducible():
    # in this example, the matrix seems a block matrix but it isn't
    A = matrix(ZZ, [[2, 2, 0, 0, 0], [0, 4, 0, 0, 4], [0, 0, 4, 4, 0], [0, 0, 0, 8, 8]])
    for diag in itertools.product([1, -1], repeat=5):
        S = diagonal_matrix(diag)
        assert shnf(A * S) == shnf(A)
