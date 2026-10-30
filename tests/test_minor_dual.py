"""
Tests for duality and taking minors.

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

import pytest
from sage.combinat.subset import powerset
from sage.matrix.constructor import matrix
from sage.matroids.basis_matroid import BasisMatroid
from sage.matroids.linear_matroid import LinearMatroid
from sage.matroids.rank_matroid import RankMatroid
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

from helpers import make_matroid, MATROID_NAMES, TORIC_NAMES, T_A, T_Q


class _FreshCopyBasisMatroid(BasisMatroid):
    # copy by returning a new plain BasisMatroid instead of self, as Sage may
    def __copy__(self):
        return BasisMatroid(groundset=self.groundset(), bases=self.bases())

    def __deepcopy__(self, memo=None):
        return BasisMatroid(groundset=self.groundset(), bases=self.bases())


class FreshCopyBasisArithmeticMatroid(BasisArithmeticMatroid, _FreshCopyBasisMatroid):
    pass

def assert_minor_of(M, N, C, D):
    """
    Assert that M is N / C \\ D.
    """
    assert M.groundset() == N.groundset() - C - D
    for X in powerset(M.groundset()):
        assert M.rank(X) == N.rank(C.union(X)) - N.rank(C), (sorted(C), sorted(D), X)
        assert M.multiplicity(X) == N.multiplicity(C.union(X)), (sorted(C), sorted(D), X)

def assert_dual_of(M, N):
    """
    Assert that M is the dual of N.
    """
    E = N.groundset()
    assert M.groundset() == E
    for X in powerset(E):
        assert M.rank(X) == len(X) + N.rank(E.difference(X)) - N.full_rank(), X
        assert M.multiplicity(X) == N.multiplicity(E.difference(X)), X

def assert_same_arithmetic_matroid(M, N):
    assert M.groundset() == N.groundset()
    for X in powerset(N.groundset()):
        assert M.rank(X) == N.rank(X), X
        assert M.multiplicity(X) == N.multiplicity(X), X

def assert_valid(M):
    # Skip rank-1 linear matroids over a non-field ring (e.g. ZZ)
    if isinstance(M, LinearMatroid) and not M.base_ring().is_field() and M.full_rank() == 1:
        return
    assert M.is_valid()

def test_minor(matroid, subtests):
    N = matroid
    N_dual = N.dual()
    for C in map(frozenset, powerset(sorted(N.groundset()))):
        for D in map(frozenset, powerset(sorted(N.groundset() - C))):
            with subtests.test(contractions=sorted(C), deletions=sorted(D)):
                M = N.minor(C, D)
                assert_minor_of(M, N, C, D)
                assert_minor_of(N.minor(contractions=C).minor(deletions=D), N, C, D)
                assert_minor_of(M.dual().dual(), N, C, D)
                assert_valid(M)
                assert_valid(M.dual())
                assert_dual_of(M.dual(), M)
                # (N / C \ D)* = N* / D \ C
                assert_same_arithmetic_matroid(M.dual(), N_dual.minor(D, C))

# expected classes of (minor, dual, minor of dual)
EXPECTED_CLASSES = {
    "rank": (MinorArithmeticMatroid, DualArithmeticMatroid, DualArithmeticMatroid),
    "U(2,5)": (MinorArithmeticMatroid, DualArithmeticMatroid, DualArithmeticMatroid),
    "basis": (BasisArithmeticMatroid,) * 3,
    "linear": (LinearArithmeticMatroid,) * 3,
    "toric": (ToricArithmeticMatroid,) * 3,
    "free-loop": (ToricArithmeticMatroid,) * 3,
    "plain-coloop": (ToricArithmeticMatroid,) * 3,
}

@pytest.mark.parametrize("name", MATROID_NAMES)
def test_minor_and_dual_classes(name):
    N = make_matroid(name)
    minor_cls, dual_cls, dual_minor_cls = EXPECTED_CLASSES[name]
    e = min(N.groundset())
    assert isinstance(N.minor(contractions=[e]), minor_cls)
    assert isinstance(N.dual(), dual_cls)
    assert isinstance(N.dual().minor(contractions=[e]), dual_minor_cls)

def test_dual(matroid):
    assert_dual_of(matroid.dual(), matroid)

@pytest.mark.parametrize("form", [
    pytest.param(lambda e: e, id="element"),
    pytest.param(lambda e: [e], id="list"),
    pytest.param(lambda e: (e,), id="tuple"),
    pytest.param(lambda e: frozenset([e]), id="frozenset"),
    pytest.param(lambda e: (x for x in [e]), id="generator")
])
def test_minor_argument_forms(matroid, form):
    e = max(matroid.groundset())
    assert_minor_of(matroid.minor(contractions=form(e)), matroid, frozenset([e]), frozenset())
    assert_minor_of(matroid.minor(deletions=form(e)), matroid, frozenset(), frozenset([e]))

def test_minor_zero_element(matroid):
    # 0 is falsey but is an element
    if 0 not in matroid.groundset():
        pytest.skip("0 is not in the ground set")
    assert_minor_of(matroid.minor(contractions=0), matroid, frozenset([0]), frozenset())

def test_minor_not_in_groundset(matroid):
    e = max(matroid.groundset()) + 1
    with pytest.raises(ValueError):
        matroid.minor(contractions=[e])

def test_minor_not_disjoint(matroid):
    e = max(matroid.groundset())
    with pytest.raises(ValueError):
        matroid.minor(contractions=[e], deletions=[e])

def test_minor_repr_and_equality():
    N = make_matroid("rank")
    M1 = N.minor(contractions=[2]).minor(deletions=[3])
    M2 = N.minor(contractions=[2], deletions=[3])
    assert M1 == M2
    assert N.minor(contractions=[2]) != N.minor(deletions=[2])  # same matroid, not the same arithmetic matroid
    assert "/ {2}" in repr(M1)
    assert "\\ {3}" in repr(M1)

@pytest.mark.parametrize("how", [copy.copy, copy.deepcopy])
@pytest.mark.parametrize("name", MATROID_NAMES + ["minor", "dual"])
def test_copy(name, how):
    if name == "minor":
        M = make_matroid("rank").minor(contractions=[2], deletions=[3])
    elif name == "dual":
        M = make_matroid("rank").dual()
    else:
        M = make_matroid(name)
    M2 = how(M)
    assert_same_arithmetic_matroid(M2, M)
    if name in ("minor", "dual"):
        assert M2 == M
        assert repr(M2) == repr(M)
        assert '_multiplicity' not in vars(M2)  # uses its own method, not M's

@pytest.mark.parametrize("name", ["rank", "U(2,5)"])
@pytest.mark.parametrize("how", [copy.copy, copy.deepcopy])
def test_copy_rank_matroid_equality(name, how):
    # Regression: copies of rank-based matroids could compare unequal to the original.
    # "rank" has a bound method as rank function, and a deep copy of it is a different method.
    M = make_matroid(name)
    M2 = how(M)
    assert M2 is not M  # our copy ran, not Sage returning M itself
    assert isinstance(M2, ArithmeticMatroid)
    assert M2 == M
    assert hash(M2) == hash(M)
    assert M2._rank_function is M._rank_function
    assert M2._multiplicity is M._multiplicity

@pytest.mark.parametrize("name", ["rank", "U(2,5)"])
def test_with_multiplicity_of_plain_rank_matroid(name):
    # when the parent's copy returns a plain RankMatroid with the same rank function
    M = make_matroid(name)
    plain = RankMatroid(M.groundset(), M._rank_function)
    assert M._with_multiplicity(plain, M._multiplicity) == M

@pytest.mark.parametrize("how", [copy.copy, copy.deepcopy])
@pytest.mark.parametrize("wrap", ["plain", "minor", "dual"])
def test_copy_keeps_multiplicity_function(how, wrap):
    # Regression: deepcopy replaced the multiplicity function with a deep copy
    T = ToricArithmeticMatroid(T_A, T_Q)
    M, twin = make_matroid("basis", T), make_matroid("basis", T)
    if wrap == "minor":
        M, twin = MinorArithmeticMatroid(M, [2]), MinorArithmeticMatroid(twin, [2])
    elif wrap == "dual":
        M, twin = DualArithmeticMatroid(M), DualArithmeticMatroid(twin)
    M2 = how(M)
    assert M2 == M
    assert M == twin  # the original is unchanged
    assert_same_arithmetic_matroid(M2, M)

@pytest.mark.parametrize("how", [copy.copy, copy.deepcopy])
def test_fresh_copy_multiplicity(how):
    N = make_matroid("basis")
    M = FreshCopyBasisArithmeticMatroid(N, multiplicity_function=N._multiplicity)
    M2 = how(M)
    assert M2 is not M
    assert M2 == M
    assert_same_arithmetic_matroid(M2, M)

@pytest.mark.parametrize("name", [x for x in MATROID_NAMES if x not in TORIC_NAMES])
def test_with_multiplicity_of_plain_matroid(name):
    # Sage currently returns arithmetic classes here, so call the fallback directly on plain matroids
    N = make_matroid(name)
    plain = {"rank": lambda: RankMatroid(N.groundset(), N.rank),
             "U(2,5)": lambda: RankMatroid(N.groundset(), N.rank),
             "basis": lambda: BasisMatroid(groundset=N.groundset(), bases=N.bases()),
             "linear": lambda: LinearMatroid(matrix=N.representation(), groundset=N.groundset_list()),
             }[name]()
    M = N._with_multiplicity(plain, N._multiplicity)
    assert isinstance(M, type(N))
    assert_same_arithmetic_matroid(M, N)

def test_dual_arithmetic_matroid():
    # test of DualArithmeticMatroid

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

    M1 = M.dual()
    assert isinstance(M1, DualArithmeticMatroid)
    assert M1.is_valid()
    assert M1 != M

    M2 = M1.dual()
    assert isinstance(M2, ArithmeticMatroid)
    assert not isinstance(M2, DualArithmeticMatroid)
    assert M == M2

    M3 = M1.minor()
    assert isinstance(M3, DualArithmeticMatroid)
    # assert M1 == M3 # this is not an equality for (non-arithmetic) matroids
    assert M3.is_valid()

def test_contract_delete():
    A = matrix(ZZ, [[2, 0, 4], [0, 3, 0]])
    M = ToricArithmeticMatroid(A)
    D = M.delete([1, 2])
    C = M.contract([1, 2])
    assert C.full_multiplicity() == 6
    assert D.full_multiplicity() == 2
    assert C.multiplicity([]) == 12
    assert D.multiplicity([]) == 1
    assert C == M.minor(contractions=[1, 2], deletions=[])
    assert D == M.minor(contractions=[], deletions=[1, 2])

def test_dual_linear_matroid():
    A = matrix(QQ, [[-1, 1, 0, -1, 2, 7], [6, 1, -1, -2, 2, 5]])

    def m(X):
        return 1

    M = LinearArithmeticMatroid(A, multiplicity_function=m)
    assert M.is_valid()

    M1 = M.dual()
    assert M1.is_valid()

    M2 = M1.dual()
    assert M2.is_valid()
    assert M != M2  # multiplicity functions are not equal in the sense of ==
    assert M.equals(M2)
    assert M != M1
    assert M1 != M2

def test_minor_linear_matroid():
    A = matrix(QQ, [[-1, 1, 0, -1, 2, 7], [6, 1, -1, -2, 2, 5]])

    def m(X):
        return 1

    M = LinearArithmeticMatroid(A, multiplicity_function=m)
    assert M.is_valid()

    M1 = M.minor(contractions=frozenset([1]), deletions=frozenset([2]))
    assert isinstance(M1, LinearArithmeticMatroid)
    assert M1.is_valid()

    M1a = M.contract([1]).delete([2])
    assert M1 != M1a  # multiplicity functions are not equal in the sense of ==
    assert M1.equals(M1a)

    M2 = M1.dual()
    assert M2.is_valid()
    assert isinstance(M2, LinearArithmeticMatroid)
    assert M != M2
    assert not M.equals(M2)

    N1 = M.dual()
    assert isinstance(N1, LinearArithmeticMatroid)
    assert N1.is_valid()

    N2 = N1.minor(contractions=frozenset([2]), deletions=frozenset([1]))
    assert N2.is_valid()
    assert isinstance(N2, LinearArithmeticMatroid)
    assert M2 != N2  # multiplicity functions are not equal in the sense of ==
    assert M2.equals(N2)

    N3 = N1.minor(contractions=frozenset([1]), deletions=frozenset([2]))
    assert N3.is_valid()
    assert isinstance(N3, LinearArithmeticMatroid)
    assert M2 != N3
    assert not M2.equals(N3)
