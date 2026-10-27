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
import unittest

from sage.matrix.constructor import matrix
from sage.matrix.special import diagonal_matrix, random_matrix
from sage.misc.prandom import choice
from sage.rings.integer_ring import ZZ

from arithmat.shnf import signed_hermite_normal_form as shnf


class TestSHNF(unittest.TestCase):

    def test_1x1(self):
        for x in range(-5, 5):
            A = matrix(ZZ, [[x]])
            self.assertEqual(shnf(A), matrix(ZZ, [[abs(x)]]))

    def test_2x2(self):
        A = matrix(ZZ, [[5, 8], [0, 3]])
        self.assertEqual(shnf(A), matrix(ZZ, [[5, 1], [0, 3]]))

    def random_test(self, r, n, all=False):
        A = random_matrix(ZZ, r, n)
        U = random_matrix(ZZ, r, r, algorithm='unimodular')

        if all:
            for diag in itertools.product([1, -1], repeat=n):
                S = diagonal_matrix(diag)
                self.assertEqual(shnf(A), shnf(U * A * S))

        else:
            S = diagonal_matrix([choice([1, -1]) for i in range(n)])
            self.assertEqual(shnf(A), shnf(U * A * S))

    def test_random_1x2(self):
        self.random_test(1, 2, all=True)

    def test_random_1x7(self):
        self.random_test(1, 7)

    def test_random_2x2(self):
        self.random_test(2, 2, all=True)

    def test_random_2x3(self):
        self.random_test(2, 3, all=True)

    def test_random_2x4(self):
        self.random_test(2, 4, all=True)

    def test_random_3x4(self):
        self.random_test(3, 4, all=True)

    def test_random_4x6(self):
        self.random_test(4, 6)

    def test_random_6x10(self):
        self.random_test(6, 10)

    def test_cardinality_four_orbit(self):
        # in these examples, an orbit of cardinality 4 occurs
        # see [PP19, Example 6.8]
        for A in [
            matrix(ZZ, [[1, 1, 1], [0, 2, 3], [0, 0, 6]]),
            matrix(ZZ, [[1, 1, 1], [0, 2, 4], [0, 0, 8]]),
        ]:
            self.assertEqual(shnf(A), A)

            for i in range(3):
                U = random_matrix(ZZ, 3, 3, algorithm='unimodular')
                for diag in itertools.product([1, -1], repeat=3):
                    S = diagonal_matrix(diag)
                    self.assertEqual(shnf(U * A * S), A)

        for u in [1, 2, 4, 5]:
            A = matrix(ZZ, [[1, 1, u], [0, 2, 3], [0, 0, 6]])
            self.assertEqual(shnf(A), matrix(ZZ, [[1, 1, 1], [0, 2, 3], [0, 0, 6]]))

    def test_is_in_hermite_form(self):
        # check that the output of shnf is in Hermite normal form
        for i in range(10):
            A = random_matrix(ZZ, 3, 4)
            B = shnf(A)
            self.assertEqual(B, B.echelon_form())

    def test_zero_pivot(self):
        # in this example, the second column has a zero pivot
        A = matrix(ZZ, [[1, -1, -1, 11], [1, -1, 0, -1], [3, -3, 7, -2]])
        U = random_matrix(ZZ, 3, 3, algorithm='unimodular')
        for diag in itertools.product([1, -1], repeat=4):
            S = diagonal_matrix(diag)
            self.assertEqual(shnf(U * A * S), shnf(A))
            self.assertEqual(shnf(U * A), shnf(A))

    def test_not_reducible(self):
        # in this example, the matrix seems a block matrix but it isn't
        A = matrix(ZZ, [[2, 2, 0, 0, 0], [0, 4, 0, 0, 4], [0, 0, 4, 4, 0], [0, 0, 0, 8, 8]])
        for diag in itertools.product([1, -1], repeat=5):
            S = diagonal_matrix(diag)
            self.assertEqual(shnf(A * S), shnf(A))
