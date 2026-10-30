"""
Tests for from_matroid.

Copyright (C) 2026 Hailey Jay Garcia

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

import pytest
from sage.combinat.subset import powerset
from sage.matrix.constructor import matrix
from sage.matroids.basis_matroid import BasisMatroid
from sage.matroids.linear_matroid import LinearMatroid
from sage.rings.integer_ring import ZZ
from sage.rings.rational_field import QQ

from arithmat.arithmetic_matroid import (
    ArithmeticMatroid,
    BasisArithmeticMatroid,
    DualArithmeticMatroid,
    LinearArithmeticMatroid,
    MinorArithmeticMatroid,
    ToricArithmeticMatroid,
)

from helpers import make_matroid, MATROID_NAMES, A_LIN


def assert_same_values(M, N):
    assert M.groundset() == N.groundset()
    for X in powerset(M.groundset()):
        assert M.rank(X) == N.rank(X), X
        assert M.multiplicity(X) == N.multiplicity(X), X


def seven(X): # methods compare equal
    return 7


@pytest.mark.parametrize("cls", [ArithmeticMatroid, BasisArithmeticMatroid])
@pytest.mark.parametrize("name", MATROID_NAMES)
def test_inherits_multiplicity(cls, name):
    N = make_matroid(name)
    M = cls.from_matroid(N)
    assert type(M) is cls
    assert M._multiplicity == N._multiplicity  # bound methods compare equal, not identical
    assert_same_values(M, N)


@pytest.mark.parametrize("cls", [ArithmeticMatroid, BasisArithmeticMatroid, LinearArithmeticMatroid])
def test_explicit_beats_source(cls):
    N = make_matroid("linear")
    M = cls.from_matroid(N, seven)
    assert M._multiplicity is seven


def test_plain_source_without_multiplicity():
    N = make_matroid("basis")
    plain = BasisMatroid(groundset=N.groundset(), bases=N.bases())
    assert BasisArithmeticMatroid.from_matroid(plain)._multiplicity is None
    assert ArithmeticMatroid.from_matroid(plain)._multiplicity is None


def test_linear_plain_source_without_multiplicity():
    L = LinearMatroid(matrix(QQ, A_LIN))
    assert LinearArithmeticMatroid.from_matroid(L)._multiplicity is None


def test_linear_inherits_from_arithmetic_source():
    N = make_matroid("linear")
    assert LinearArithmeticMatroid.from_matroid(N)._multiplicity is N._multiplicity


@pytest.mark.parametrize("name", ["toric", "free-loop", "plain-coloop"])
def test_toric_from_toric(name):
    T = make_matroid(name)
    M = ToricArithmeticMatroid.from_matroid(T)
    assert M is not T
    assert M == T
    assert_same_values(M, T)


def test_toric_from_linear():
    L = LinearMatroid(matrix(QQ, A_LIN))
    M = ToricArithmeticMatroid.from_matroid(L)
    assert_same_values(M, ToricArithmeticMatroid(A_LIN))


def test_toric_rejects_multiplicity():
    with pytest.raises(ValueError):
        ToricArithmeticMatroid.from_matroid(make_matroid("toric"), seven)


def test_toric_rejects_non_integral_representation():
    L = LinearMatroid(matrix(QQ, [[1, 1/2], [0, 1]]))
    with pytest.raises(TypeError):
        ToricArithmeticMatroid.from_matroid(L)


def test_basis_constructor_explicit_beats_source():
    N = make_matroid("basis")
    assert BasisArithmeticMatroid(N, multiplicity_function=seven)._multiplicity is seven
    assert BasisArithmeticMatroid(N)._multiplicity is N._multiplicity


class SubArithmetic(ArithmeticMatroid):
    pass


class SubBasis(BasisArithmeticMatroid):
    pass


class SubLinear(LinearArithmeticMatroid):
    pass


class SubToric(ToricArithmeticMatroid):
    pass


@pytest.mark.parametrize("cls, name", [(SubArithmetic, "rank"), (SubBasis, "basis"), (SubLinear, "linear"),
                                       (SubToric, "toric")])
def test_subclass_type(cls, name):
    assert type(cls.from_matroid(make_matroid(name))) is cls


@pytest.mark.parametrize("how", [copy.copy, copy.deepcopy])
def test_subclass_copies_keep_type(how):
    N = make_matroid("rank")
    for M in (SubArithmetic.from_matroid(N), SubToric.from_matroid(make_matroid("toric"))):
        M2 = how(M)
        assert type(M2) is type(M)
        assert M2 == M

    class SubMinor(MinorArithmeticMatroid):
        pass

    class SubDual(DualArithmeticMatroid):
        pass

    for M in (SubMinor(N, [0]), SubDual(N)):
        M2 = how(M)
        assert type(M2) is type(M)
        assert M2 == M
