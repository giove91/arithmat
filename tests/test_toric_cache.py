"""
Tests for caching in ToricArithmeticMatroids.

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
import operator
import unittest
from functools import reduce

from sage.combinat.subset import powerset
from sage.matrix.constructor import matrix
from sage.matrix.special import block_matrix
from sage.rings.integer_ring import ZZ

from arithmat.arithmetic_matroid import ToricArithmeticMatroid


def cache_size(M):
    return len(M._rank_and_multiplicity.cache)

def reference_rank_and_multiplicity(M, X):
    """
    Pre-cache formula for (rank, multiplicity)
    """
    cols = [M._groundset_to_index[e] for e in X]
    T = block_matrix(ZZ, [[M._A[:, cols], M._Q]])
    rank = T.rank() - M._Q.ncols()
    multiplicity = reduce(operator.mul, [d for d in T.elementary_divisors() if d != 0], 1)
    return rank, multiplicity

class TestToricCache(unittest.TestCase):

    def setUp(self):
        self.A = matrix(ZZ, [[1, 0, 0, 2, 1, 3], [0, 1, 0, 1, 4, -1], [0, 0, 2, 1, -1, 5]])
        self.Q_kept = matrix(ZZ, [[4], [0], [0]])  # survives _normalize
        self.Q_stripped = matrix(ZZ, [[1], [0], [0]])  # removed by _normalize

    def matroids(self):
        return [
            ToricArithmeticMatroid(self.A),
            ToricArithmeticMatroid(self.A, self.Q_kept),
            ToricArithmeticMatroid(self.A, self.Q_stripped),
            ToricArithmeticMatroid(matrix(ZZ, [[-1, 1, 0, -1, 2, 7], [6, 1, -1, -2, 2, 5]])),
        ]

    def test_matches_reference_on_all_subsets(self):
        for M in self.matroids():
            for X in powerset(M.groundset()):
                X = frozenset(X)
                expected = reference_rank_and_multiplicity(M, X)
                self.assertEqual((M.rank(X), M.multiplicity(X)), expected)
                self.assertEqual((M._rank(X), M._multiplicity(X)), expected)

    def test_empty_set(self):
        M = ToricArithmeticMatroid(self.A)
        self.assertEqual(M._Q.ncols(), 0)
        self.assertEqual((M.rank([]), M.multiplicity([])), (0, 1))
        self.assertTrue(M.is_torsion_free())

        M = ToricArithmeticMatroid(self.A, self.Q_kept)
        self.assertEqual(M._Q.ncols(), 1)
        self.assertEqual((M.rank([]), M.multiplicity([])), (0, 4))
        self.assertFalse(M.is_torsion_free())

    def test_argument_types_share_one_entry(self):
        M = ToricArithmeticMatroid(self.A)
        expected = reference_rank_and_multiplicity(M, [1, 2])
        for X in ([1, 2], [2, 1], {1, 2}, (1, 2), frozenset([1, 2])):
            self.assertEqual((M._rank(X), M._multiplicity(X)), expected)
            self.assertEqual((M.rank(X), M.multiplicity(X)), expected)
        self.assertEqual(cache_size(M), 1)

    def test_hash(self):
        # __hash__ calls _multiplicity; must not recurse through the cache
        M1 = ToricArithmeticMatroid(self.A, self.Q_kept)
        M2 = ToricArithmeticMatroid(self.A, self.Q_kept)
        self.assertEqual(hash(M1), hash(M2))
        self.assertEqual(len({M1, M2}), 1)

    def test_cache_fills_once(self):
        M = ToricArithmeticMatroid(self.A)
        n = len(M.groundset())
        self.assertEqual(cache_size(M), 0)

        T1 = M.arithmetic_tutte_polynomial()
        self.assertEqual(cache_size(M), 2 ** n)

        T2 = M.arithmetic_tutte_polynomial()
        self.assertEqual(cache_size(M), 2 ** n)
        self.assertEqual(T1, T2)

    def test_subset_table_matches_reference(self):
        for M in self.matroids():
            table = M._subset_table()
            self.assertEqual(len(table), 2 ** len(M.groundset()))
            for X, value in table.items():
                self.assertEqual(value, reference_rank_and_multiplicity(M, X))

    def test_derived_matroids_have_own_cache(self):
        M = ToricArithmeticMatroid(self.A, self.Q_kept)
        M.arithmetic_tutte_polynomial()  # fill M's cache
        size = cache_size(M)

        derived = [
            M.minor(contractions=[1]),
            M.minor(deletions=[0, 3]),
            M.dual(),
            copy.copy(M),
            copy.deepcopy(M),
        ]
        for N in derived:
            self.assertEqual(cache_size(N), 0)
            for X in powerset(N.groundset()):
                X = frozenset(X)
                self.assertEqual((N.rank(X), N.multiplicity(X)), reference_rank_and_multiplicity(N, X))
            self.assertEqual(cache_size(M), size)

    def test_isomorphism_with_cache(self):
        U = matrix(ZZ, [[1, 1, 0], [0, 1, 0], [2, 3, 1]])  # unimodular
        perm = list(reversed(range(self.A.ncols())))
        M = ToricArithmeticMatroid(self.A)
        N = ToricArithmeticMatroid(U * self.A[:, perm])
        self.assertTrue(M._is_isomorphic(N))
        self.assertTrue(M._is_isomorphic(N))  # second call reads the filled caches

        B = copy.copy(self.A)
        B[2, 2] = 3
        self.assertFalse(M._is_isomorphic(ToricArithmeticMatroid(B)))
