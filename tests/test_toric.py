"""
Tests for ToricArithmeticMatroid.

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
import itertools

import pytest
from sage.combinat.posets.posets import Poset
from sage.combinat.set_partition import SetPartition
from sage.combinat.subset import powerset
from sage.matrix.constructor import matrix
from sage.rings.integer_ring import ZZ
from sage.homology.homology_group import HomologyGroup

from arithmat.arithmetic_matroid import ToricArithmeticMatroid
from arithmat.shnf import signed_hermite_normal_form

def assert_mutable_copy(R, M):
    """
    Assert that R is a new mutable matrix, and that changing it leaves M as it was.
    """
    twin = copy.deepcopy(M)
    assert R.is_mutable()
    assert R is not M._A
    assert R is not M._Q
    if R.nrows() and R.ncols():
        R[0, 0] += 100
    assert not M._A.is_mutable()
    assert not M._Q.is_mutable()
    assert M == twin
    assert M.arrangement_matrix() == twin.arrangement_matrix()


def test_without_Q():
    A = matrix(ZZ, [[-1, 1, 0, -1, 2, 7], [6, 1, -1, -2, 2, 5]])
    M = ToricArithmeticMatroid(A)

    assert M.full_rank() == M.rank(M.groundset())

    assert M._Q == matrix(ZZ, 2, 0)
    assert M.rank([0]) == 1
    assert M.rank([0, 1]) == 2
    assert M.rank([0, 1, 2]) == 2

    assert M._multiplicity([1, 2]) == 1

    # minor
    M2 = M.minor(contractions=[], deletions=[1])
    assert M2.groundset() == frozenset([0, 2, 3, 4, 5])
    assert M2._A == matrix(ZZ, [[-1, 0, -1, 2, 7], [6, -1, -2, 2, 5]])

    M2 = M.minor(contractions=[1], deletions=[])
    assert M2.groundset() == frozenset([0, 2, 3, 4, 5])
    assert M2._Q == matrix(ZZ, [[]])
    assert M2._multiplicity([0, 2]) == 1

    # dual
    M1 = M.dual()
    assert M1.rank(M1.groundset()) == 4
    M2 = M1.dual()
    assert M.equals(M2)

    # check representation
    assert M.check_representation(A)
    assert M.is_representable()
    assert M.representation(shnf=False) == A
    assert M.representation(ordered_groundset=[5, 4, 3, 2, 1, 0], shnf=False) == \
        matrix(ZZ, [[7, 2, -1, 0, 1, -1], [5, 2, -2, -1, 1, 6]])
    assert M.representation() == signed_hermite_normal_form(A)

    # check orientability
    assert M.is_orientable()

@pytest.mark.parametrize("A,Q,ncols,mult_checks,rank_checks",[
pytest.param(matrix(ZZ, [[-1, 1, 0, -1, 2, 7], [6, 1, -1, -2, 2, 5]], immutable=True),
            matrix(ZZ, [[5, 9, 1], [-3, -2, -1]], immutable=True), 0, (([0], 1),), ()),
pytest.param(matrix(ZZ, [[-1, 1, 0, -1, 2, 7], [0, 1, -1, -2, 2, 5]], immutable=True),
            matrix(ZZ, [[3, 9, 6], [0, 0, 0]]), 1, (([0], 1), ([1], 3), ([1,2], 1)),
            (([0], 0), ([1, 2], 1))),
pytest.param(matrix(ZZ, [[-1, 1, 0, 7], [6, 1, -1, -2]]),
            matrix(ZZ, [[5], [1]]), 0, (), ())
])
def test_with_Q(A, Q, ncols, mult_checks, rank_checks):
    M = ToricArithmeticMatroid(A, torus_matrix=Q)
    assert M.full_rank() == M.rank(M.groundset())
    assert M.full_multiplicity() == M.multiplicity()
    assert M._Q.ncols() == ncols
    for (ls, mult) in mult_checks:
        assert M._multiplicity(ls) == mult
    for (ls, rk) in rank_checks:
        assert M.rank(ls) == rk

def test_minor():
    A = matrix(ZZ, [[-1, 1, 0, -1, 2, 7], [6, 1, -1, -2, 2, 5]])
    M = ToricArithmeticMatroid(A)

    # minor
    M2 = M.minor(contractions=[], deletions=[1])
    assert M2.groundset() == frozenset([0, 2, 3, 4, 5])
    assert M2._A == matrix(ZZ, [[-1, 0, -1, 2, 7], [6, -1, -2, 2, 5]])

    M2 = M.minor(contractions=[1], deletions=[])
    assert M2.groundset() == frozenset([0, 2, 3, 4, 5])
    assert M2._Q == matrix(ZZ, [[]])
    assert M2.multiplicity([0, 2]) == 1

def test_minor2():
    A = matrix(ZZ, [[2, 1], [0, 3]])
    M = ToricArithmeticMatroid(A).minor(contractions=[0], deletions=[])
    N = ToricArithmeticMatroid(A).minor(contractions=[], deletions=[0])
    assert not M.is_torsion_free()
    assert N.is_torsion_free()
    assert not N.is_isomorphic(M)

def test_ordered_groundset():
    A = matrix(ZZ, [[-1, 1, 0, -1, 2, 7], [6, 1, -1, -2, 2, 5]])
    M = ToricArithmeticMatroid(A, ordered_groundset=['a', 'b', 'f', 'e', 'd', 'c'])

    assert M.groundset() == frozenset(['a', 'b', 'f', 'e', 'd', 'c'])
    assert M.full_rank() == M.rank(M.groundset())

    assert M._Q == matrix(ZZ, 2, 0)
    assert M.rank(['a']) == 1
    assert M.rank(['a', 'b']) == 2
    assert M.rank(['f', 'a', 'b']) == 2

    assert M._multiplicity(['b', 'f']) == 1

    # test failure if there are repetition in the ordered groundset
    with pytest.raises(AssertionError):
        M = ToricArithmeticMatroid(A, ordered_groundset=['a', 'b', 'f', 'e', 'd', 'b'])

@pytest.mark.parametrize("A, expect_equiv", [
    (matrix(ZZ, [[1, 1, 2], [0, 7, 7]], immutable=True), True),
    (matrix(ZZ, [[1, 2, 3], [0, 7, 7]], immutable=True), False)
])
@pytest.mark.parametrize("ordered_groundset, morphism", [
    pytest.param(None, None, id="same-groundset"),
    pytest.param(["a", "b", "c"], {0:"a", 1:"b", 2:"c"}, id="with-morphism")
])
def test_equivalence(A, expect_equiv, ordered_groundset, morphism):
    # without Q

    # invert
    msihprom = None if morphism is None else {val:key for key,val in morphism.items()}

    M = ToricArithmeticMatroid(A)
    N = ToricArithmeticMatroid(matrix(ZZ, [[-1, -1, 2], [-1, 6, -5]]), ordered_groundset=ordered_groundset)

    assert M.is_equivalent(N, morphism) == expect_equiv
    assert N.is_equivalent(M, msihprom) == expect_equiv

    assert M.is_isomorphic(N)

def test_non_equivalent_representations():
    M = ToricArithmeticMatroid(matrix(ZZ, [[6, 3, -2, 2], [3, 21, 0, -9], [-1, -4, 3, -2]]))

    assert M.num_representations() == 4
    matroids = [ToricArithmeticMatroid(A) for A in M.all_representations()]
    for N, O in itertools.combinations(matroids, 2):
        assert not N.is_equivalent(O)
        assert N.is_isomorphic(O)

def test_all_representations():
    A = matrix(ZZ, [[1, 1, 4, 1, 0], [0, 2, 3, 1, 1], [0, 0, 6, 0, 1]])
    M = ToricArithmeticMatroid(A)
    representations = list(M.all_representations(shnf=False))
    [B] = representations

    shnf_representations = list(M.all_representations())
    [C] = shnf_representations
    assert C == M.representation()

    assert signed_hermite_normal_form(B) == C
    assert signed_hermite_normal_form(A) == C

def test_ordered_groundset_is_list_copy():
    A = matrix(ZZ, [[-1, 1, 0, -1, 2, 7], [6, 1, -1, -2, 2, 5]])
    E = ['a', 'b', 'f', 'e', 'd', 'c']
    M = ToricArithmeticMatroid(A, ordered_groundset=E)
    L = M.ordered_groundset()
    assert isinstance(L, list)
    assert L == E
    assert M.ordered_groundset() is not L  # new list on each call, like Sage's groundset_list()
    L.reverse()
    L.append('z')
    assert M.ordered_groundset() == E
    assert M.representation(shnf=False) == A

A_ACCESSORS = matrix(ZZ, [[-1, 1, 0, -1, 2, 7], [0, 1, -1, -2, 2, 5]], immutable=True)

@pytest.mark.parametrize("Q", [
    pytest.param(None, id="without-Q"),
    pytest.param(matrix(ZZ, [[3, 9, 6], [0, 0, 0]], immutable=True), id="with-Q")
])
@pytest.mark.parametrize("accessor, internal", [
    ("arrangement_matrix", "_A"),
    ("torus_matrix", "_Q")
])
def test_matrix_accessors_return_mutable_copies(Q, accessor, internal):
    M = ToricArithmeticMatroid(A_ACCESSORS) if Q is None else ToricArithmeticMatroid(A_ACCESSORS, Q)
    rank_and_multiplicity = [(M.rank(X), M.multiplicity(X)) for X in powerset(M.groundset())]
    R = getattr(M, accessor)()
    assert R == getattr(M, internal)
    assert_mutable_copy(R, M)
    # the cache of rank and multiplicity is still valid
    assert [(M.rank(X), M.multiplicity(X)) for X in powerset(M.groundset())] == rank_and_multiplicity

@pytest.mark.parametrize("shnf", [True, False])
@pytest.mark.parametrize("E", [
    pytest.param(None, id="default-order"),
    pytest.param([9, 8, 7, 6, 5], id="reversed")
])
def test_representations_are_mutable_copies(shnf, E):
    A = matrix(ZZ, [[1, 1, 4, 1, 0], [0, 2, 3, 1, 1], [0, 0, 6, 0, 1]])
    M = ToricArithmeticMatroid(A, ordered_groundset=[5, 6, 7, 8, 9])
    R = M.representation(ordered_groundset=E, shnf=shnf)
    assert R == M.arithmetic_representation(ordered_groundset=E, shnf=shnf)
    assert_mutable_copy(R, M)

@pytest.mark.parametrize("shnf", [True, False])
@pytest.mark.parametrize("rows, E", [
    pytest.param([[1, 1, 4, 1, 0], [0, 2, 3, 1, 1], [0, 0, 6, 0, 1]], [5, 6, 7, 8, 9], id="full-mult-1"),
    pytest.param([[1, 1, 0], [0, 2, 2]], None, id="full-mult-2")
])
def test_all_representations_are_mutable_copies(shnf, rows, E):
    N = ToricArithmeticMatroid(matrix(ZZ, rows), ordered_groundset=E)
    # full multiplicity 1 and 2 take different branches of all_representations()
    expected = list(N.all_representations(shnf=shnf))
    assert expected
    for R in N.all_representations(shnf=shnf):
        assert_mutable_copy(R, N)
        assert list(N.all_representations(shnf=shnf)) == expected

def test_representation_uses_toric_override():
    # representation() goes through arithmetic_representation(), so Toric's version (which keeps A) is used
    A = matrix(ZZ, [[1, 1, 4, 1, 0], [0, 2, 3, 1, 1], [0, 0, 6, 0, 1]])
    M = ToricArithmeticMatroid(A)
    assert M.representation(shnf=False) == A
    assert M.arithmetic_representation(shnf=False) == A
    assert M.representation() == signed_hermite_normal_form(A)

@pytest.mark.parametrize("rows, blocks, decomposable", [
    pytest.param([[1, 5, -1, 5], [1, 7, 2, 6], [0, -2, -2, -1]], [[2], [0, 1, 3]], True, id="two-blocks"),
    pytest.param([[1, 5, -10, 50], [0, 0, 2, 6], [0, 0, -2, -1]], [[0, 1], [2, 3]], True, id="block-diagonal"),
    pytest.param([[1, 1], [0, 2]], [[0, 1]], False, id="indecomposable"),
])
def test_decomposition(rows, blocks, decomposable):
    M = ToricArithmeticMatroid(matrix(ZZ, rows))
    assert M.decomposition() == SetPartition(blocks)
    assert M.is_decomposable() == decomposable

@pytest.mark.parametrize("ordered_groundset,expected_labels", [
    pytest.param(None, None, id="no-groundset"),
    pytest.param(['a', 'c', 'b'], {(), ('a',), ('b',), ('c',), ('a', 'c', 'b'), ('a', 'b')}, id="alpha-groundset")
])
def test_poset_of_layers(ordered_groundset, expected_labels):
    M = ToricArithmeticMatroid(matrix(ZZ, [[1, 0, 1], [0, 1, 3]]), ordered_groundset=ordered_groundset)
    P = M.poset_of_layers()
    assert len(P) == 7
    assert P.has_bottom()
    assert not P.has_top()
    assert len(P.cover_relations()) == 10
    assert P.is_isomorphic(
        Poset(data=(range(7), [(0, 1), (0, 2), (0, 3), (1, 4), (2, 4), (2, 5), (2, 6), (3, 4), (3, 5), (3, 6)])))
    if expected_labels is not None:
        assert set(S for (S, x) in P) == expected_labels

# see [DGP17]
@pytest.mark.parametrize("rows, betti", [
    pytest.param([[2, 0, 1, 1], [0, 2, 1, -1]], [0, 3], id="C2"),
    pytest.param([[2, 0, 0, 1, 1, 1, 1, 0, 0], [0, 2, 0, 1, -1, 0, 0, 1, 1], [0, 0, 2, 0, 0, 1, -1, 1, -1]],
                 [0, 0, 15], id="C3"),
    pytest.param([[1, 0, 1, 1], [0, 1, 1, -1]], [0, 1], id="B2"),
    pytest.param([[1, 0, 0, 1, 1, 1, 1, 0, 0], [0, 1, 0, 1, -1, 0, 0, 1, 1], [0, 0, 1, 0, 0, 1, -1, 1, -1]],
                 [0, 0, 6], id="B3"),
    pytest.param([[1, 1, 1, 1, 0, 0], [1, -1, 0, 0, 1, 1], [0, 0, 1, -1, 1, -1]], [0, 0, 3], id="D3"),
])
def test_poset_of_layers_root_system(rows, betti):
    M = ToricArithmeticMatroid(matrix(ZZ, rows))
    P = M.poset_of_layers()
    homology = P.subposet([a for a in P if a != P.bottom()]).order_complex(on_ints=True).homology()
    assert [homology[i].ngens() for i in range(len(betti))] == betti
    assert homology[len(betti) - 1].ngens() == M.arithmetic_tutte_polynomial()(x=0, y=0)


def test_poset_of_layers_root_system_C2_size():
    # see [DGP17]
    P = ToricArithmeticMatroid(matrix(ZZ, [[2, 0, 1, 1], [0, 2, 1, -1]])).poset_of_layers()
    assert len(P) == 11
    assert len(P.cover_relations()) == 18


def test_different_posets_of_layers():
    # see [PP19, Section 8]
    A = matrix(ZZ, [[1, 1, 1, 3], [0, 5, 0, 5], [0, 0, 5, 5]])
    B = matrix(ZZ, [[1, 4, 1, 6], [0, 5, 0, 5], [0, 0, 5, 5]])
    C = matrix(ZZ, [[1, 2, 2, 1], [0, 5, 0, 5], [0, 0, 5, -5]])

    # check signed Hermite normal form
    assert signed_hermite_normal_form(A) == matrix(ZZ, [[1, 1, 1, -3], [0, 5, 0, -5], [0, 0, 5, -5]])
    assert signed_hermite_normal_form(B) == matrix(ZZ, [[1, 1, 1, -1], [0, 5, 0, 5], [0, 0, 5, -5]])
    assert signed_hermite_normal_form(C) == matrix(ZZ, [[1, 2, 2, 1], [0, 5, 0, 5], [0, 0, 5, -5]])

    MA = ToricArithmeticMatroid(A)
    MB = ToricArithmeticMatroid(B)
    MC = ToricArithmeticMatroid(C)

    assert MA.is_isomorphic(MB)
    assert MA.is_isomorphic(MC)
    assert MB.is_isomorphic(MC)

    assert not MA.poset_of_layers().is_isomorphic(MB.poset_of_layers())
    assert not MA.poset_of_layers().is_isomorphic(MC.poset_of_layers())
    assert not MB.poset_of_layers().is_isomorphic(MC.poset_of_layers())

@pytest.mark.parametrize("ordered_groundset,expected_labels", [
    pytest.param(None, None, id="no-groundset"),
    pytest.param(['a', 'c', 'b'], {(), ('a',), ('b',), ('c',), ('a', 'c'), ('a', 'b'), ('c', 'b')}, id="alpha-groundset")
])
def test_independence_poset(ordered_groundset, expected_labels):
    A = matrix(ZZ, [[1, 0, 1], [0, 1, 3]])
    M = ToricArithmeticMatroid(A, ordered_groundset=ordered_groundset)
    P = M.arithmetic_independence_poset()

    assert len(P) == 9
    assert len(P.cover_relations()) == 13
    assert P.rank() == 2
    assert P.has_bottom()
    if expected_labels is not None:
        assert set(S for (S, x) in P) == expected_labels

@pytest.mark.parametrize("poset, expected_homology", [
    pytest.param
    ("poset_of_layers", {
            0: HomologyGroup(0, ZZ),  # trivial group
            1: HomologyGroup(1, ZZ, [5]),  # Z/5Z
            2: HomologyGroup(48, ZZ),  # Z^48
        }, id="poset_of_layers"),
    pytest.param("arithmetic_independence_poset", {
            0: HomologyGroup(0, ZZ),  # trivial group
            1: HomologyGroup(1, ZZ, [5]),  # Z/5Z
            2: HomologyGroup(73, ZZ),  # Z^73
        }, id="arithmetic_independence_poset")
])
def test_nonshellable_posets(poset, expected_homology):
    # see [PP19, Section 8]
    A = matrix(ZZ, [[1, 1, 1, -3], [0, 5, 0, -5], [0, 0, 5, -5]])
    B = matrix(ZZ, [[1, 1, 1, -1], [0, 5, 0, 5], [0, 0, 5, -5]])
    C = matrix(ZZ, [[1, 2, 2, 1], [0, 5, 0, 5], [0, 0, 5, -5]])

    for X in [A, B, C]:
        M = ToricArithmeticMatroid(X)
        P = getattr(M, poset)()
        homology = P.subposet([a for a in P if a != P.bottom()]).order_complex(on_ints=True).homology()

        assert homology == expected_homology

# Disconnected matroids: signs of whole components must also be considered
@pytest.mark.parametrize("A, B, expect_equiv", [
    pytest.param(matrix(ZZ, [[-2, 2, 2, 2], [1, -1, -1, 3]]),
                 matrix(ZZ, [[3, -3, -3, -65], [2, -2, -2, -46]]), True, id="with-coloop"),
    pytest.param(matrix(ZZ, [[1, -1, -1, 3, 1], [-1, 2, -2, -3, -1], [2, -2, -2, -3, -1]]),
                 matrix(ZZ, [[5, 5, -5, 3, 1], [15, 18, -6, 0, 0], [-9, -14, -6, 9, 3]]), True, id="two-components"),
    pytest.param(matrix(ZZ, [[1, -1, 1], [0, 0, 5]]),
                 matrix(ZZ, [[1, -1, 2], [0, 0, 5]]), False, id="not-equivalent"),
])
def test_equivalence_disconnected(A, B, expect_equiv):
    M = ToricArithmeticMatroid(A)
    N = ToricArithmeticMatroid(B)

    assert len(M.components()) > 1
    assert M.is_isomorphic(N)
    assert M.is_equivalent(N) == expect_equiv
    assert N.is_equivalent(M) == expect_equiv
