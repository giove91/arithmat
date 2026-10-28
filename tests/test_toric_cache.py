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
from functools import reduce

import pytest
from sage.combinat.subset import powerset
from sage.matrix.constructor import matrix
from sage.matrix.special import block_matrix
from sage.rings.integer_ring import ZZ

from arithmat.arithmetic_matroid import ToricArithmeticMatroid

A = matrix(ZZ, [[1, 0, 0, 2, 1, 3], [0, 1, 0, 1, 4, -1], [0, 0, 2, 1, -1, 5]], immutable=True)
Q_KEPT = matrix(ZZ, [[4], [0], [0]], immutable=True)  # survives _normalize
Q_STRIPPED = matrix(ZZ, [[1], [0], [0]], immutable=True)  # removed by _normalize


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


@pytest.fixture(params=["no-Q", "Q-kept", "Q-stripped", "2x6"])
def toric(request):
    # A new matroid for each test, so no test sees a cache another test filled
    return {
        "no-Q": lambda: ToricArithmeticMatroid(A),
        "Q-kept": lambda: ToricArithmeticMatroid(A, Q_KEPT),
        "Q-stripped": lambda: ToricArithmeticMatroid(A, Q_STRIPPED),
        "2x6": lambda: ToricArithmeticMatroid(matrix(ZZ, [[-1, 1, 0, -1, 2, 7], [6, 1, -1, -2, 2, 5]])),
    }[request.param]()


def test_matches_reference_on_all_subsets(toric):
    for X in powerset(toric.groundset()):
        X = frozenset(X)
        expected = reference_rank_and_multiplicity(toric, X)
        assert (toric.rank(X), toric.multiplicity(X)) == expected
        assert (toric._rank(X), toric._multiplicity(X)) == expected


def test_subset_table_matches_reference(toric):
    table = toric._subset_table()
    assert len(table) == 2 ** len(toric.groundset())
    for X, value in table.items():
        assert value == reference_rank_and_multiplicity(toric, X)


def test_empty_set():
    M = ToricArithmeticMatroid(A)
    assert M._Q.ncols() == 0
    assert (M.rank([]), M.multiplicity([])) == (0, 1)
    assert M.is_torsion_free()

    M = ToricArithmeticMatroid(A, Q_KEPT)
    assert M._Q.ncols() == 1
    assert (M.rank([]), M.multiplicity([])) == (0, 4)
    assert not M.is_torsion_free()


def test_argument_types_share_one_entry():
    M = ToricArithmeticMatroid(A)
    expected = reference_rank_and_multiplicity(M, [1, 2])
    for X in ([1, 2], [2, 1], {1, 2}, (1, 2), frozenset([1, 2])):
        assert (M._rank(X), M._multiplicity(X)) == expected
        assert (M.rank(X), M.multiplicity(X)) == expected
    assert cache_size(M) == 1


def test_hash():
    # __hash__ calls _multiplicity; must not recurse through the cache
    M1 = ToricArithmeticMatroid(A, Q_KEPT)
    M2 = ToricArithmeticMatroid(A, Q_KEPT)
    assert hash(M1) == hash(M2)
    assert len({M1, M2}) == 1


def test_cache_fills_once():
    M = ToricArithmeticMatroid(A)
    n = len(M.groundset())
    assert cache_size(M) == 0

    T1 = M.arithmetic_tutte_polynomial()
    assert cache_size(M) == 2 ** n

    T2 = M.arithmetic_tutte_polynomial()
    assert cache_size(M) == 2 ** n
    assert T1 == T2


def test_derived_matroids_have_own_cache():
    M = ToricArithmeticMatroid(A, Q_KEPT)
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
        assert cache_size(N) == 0
        for X in powerset(N.groundset()):
            X = frozenset(X)
            assert (N.rank(X), N.multiplicity(X)) == reference_rank_and_multiplicity(N, X)
        assert cache_size(M) == size


def test_isomorphism_with_cache():
    U = matrix(ZZ, [[1, 1, 0], [0, 1, 0], [2, 3, 1]])  # unimodular
    perm = list(reversed(range(A.ncols())))
    M = ToricArithmeticMatroid(A)
    N = ToricArithmeticMatroid(U * A[:, perm])
    assert M._is_isomorphic(N)
    assert M._is_isomorphic(N)  # second call reads the filled caches

    B = copy.copy(A)
    B[2, 2] = 3
    assert not M._is_isomorphic(ToricArithmeticMatroid(B))
