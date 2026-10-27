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
import unittest

from sage.combinat.posets.posets import Poset
from sage.combinat.set_partition import SetPartition
from sage.combinat.subset import powerset
from sage.matrix.constructor import matrix
from sage.rings.integer_ring import ZZ

from arithmat.arithmetic_matroid import ToricArithmeticMatroid
from arithmat.shnf import signed_hermite_normal_form


class TestToric(unittest.TestCase):

    def test_without_Q(self):
        A = matrix(ZZ, [[-1, 1, 0, -1, 2, 7], [6, 1, -1, -2, 2, 5]])
        M = ToricArithmeticMatroid(A)

        self.assertEqual(M.full_rank(), M.rank(M.groundset()))

        self.assertEqual(M._Q, matrix(ZZ, 2, 0))
        self.assertEqual(M.rank([0]), 1)
        self.assertEqual(M.rank([0, 1]), 2)
        self.assertEqual(M.rank([0, 1, 2]), 2)

        self.assertEqual(M._multiplicity([1, 2]), 1)

        # minor
        M2 = M.minor(contractions=[], deletions=[1])
        self.assertEqual(M2.groundset(), frozenset([0, 2, 3, 4, 5]))
        self.assertEqual(M2._A, matrix(ZZ, [[-1, 0, -1, 2, 7], [6, -1, -2, 2, 5]]))

        M2 = M.minor(contractions=[1], deletions=[])
        self.assertEqual(M2.groundset(), frozenset([0, 2, 3, 4, 5]))
        self.assertEqual(M2._Q, matrix(ZZ, [[]]))
        self.assertEqual(M2._multiplicity([0, 2]), 1)

        # dual
        M1 = M.dual()
        self.assertEqual(M1.rank(M1.groundset()), 4)
        M2 = M1.dual()
        self.assertTrue(M.equals(M2))

        # check representation
        self.assertTrue(M.check_representation(A))
        self.assertTrue(M.is_representable())
        self.assertEqual(M.representation(shnf=False), A)
        self.assertEqual(
            M.representation(ordered_groundset=[5, 4, 3, 2, 1, 0], shnf=False),
            matrix(ZZ, [[7, 2, -1, 0, 1, -1], [5, 2, -2, -1, 1, 6]])
        )
        self.assertEqual(M.representation(), signed_hermite_normal_form(A))

        # check orientability
        self.assertTrue(M.is_orientable())

    def test_with_Q(self):
        A = matrix(ZZ, [[-1, 1, 0, -1, 2, 7], [6, 1, -1, -2, 2, 5]])
        Q = matrix(ZZ, [[5, 9, 1], [-3, -2, -1]])
        M = ToricArithmeticMatroid(A, torus_matrix=Q)

        self.assertEqual(M.full_rank(), M.rank(M.groundset()))

        self.assertEqual(M._Q.ncols(), 0)
        self.assertEqual(M._multiplicity([0]), 1)

    def test_with_Q_2(self):
        A = matrix(ZZ, [[-1, 1, 0, -1, 2, 7], [0, 1, -1, -2, 2, 5]])
        Q = matrix(ZZ, [[3, 9, 6], [0, 0, 0]])
        M = ToricArithmeticMatroid(A, torus_matrix=Q)

        self.assertEqual(M.full_rank(), M.rank(M.groundset()))

        self.assertEqual(M._Q, matrix(ZZ, [[3], [0]]))
        self.assertEqual(M.rank([0]), 0)
        self.assertEqual(M.rank([1, 2]), 1)

        self.assertEqual(M._multiplicity([0]), 1)
        self.assertEqual(M._multiplicity([1]), 3)
        self.assertEqual(M._multiplicity([1, 2]), 1)

    def test_with_Q_3(self):
        A = matrix(ZZ, [[-1, 1, 0, 7], [6, 1, -1, -2]])
        Q = matrix(ZZ, [[5], [1]])
        M = ToricArithmeticMatroid(A, Q)

        self.assertEqual(M.full_rank(), M.rank(M.groundset()))
        self.assertEqual(M._Q.ncols(), 0)

    def test_minor(self):
        A = matrix(ZZ, [[-1, 1, 0, -1, 2, 7], [6, 1, -1, -2, 2, 5]])
        M = ToricArithmeticMatroid(A)

        self.assertEqual(M.full_rank(), M.rank(M.groundset()))

        self.assertEqual(M._Q, matrix(ZZ, 2, 0))
        self.assertEqual(M.rank([0]), 1)
        self.assertEqual(M.rank([0, 1]), 2)
        self.assertEqual(M.rank([0, 1, 2]), 2)

        self.assertEqual(M.multiplicity([1, 2]), 1)

        # minor
        M2 = M.minor(contractions=[], deletions=[1])
        self.assertEqual(M2.groundset(), frozenset([0, 2, 3, 4, 5]))
        self.assertEqual(M2._A, matrix(ZZ, [[-1, 0, -1, 2, 7], [6, -1, -2, 2, 5]]))

        M2 = M.minor(contractions=[1], deletions=[])
        self.assertEqual(M2.groundset(), frozenset([0, 2, 3, 4, 5]))
        self.assertEqual(M2._Q, matrix(ZZ, [[]]))
        self.assertEqual(M2.multiplicity([0, 2]), 1)

    def test_minor2(self):
        A = matrix(ZZ, [[2, 1], [0, 3]])
        M = ToricArithmeticMatroid(A).minor(contractions=[0], deletions=[])
        N = ToricArithmeticMatroid(A).minor(contractions=[], deletions=[0])
        self.assertFalse(M.is_torsion_free())
        self.assertTrue(N.is_torsion_free())
        self.assertFalse(N.is_isomorphic(M))

    def test_ordered_groundset(self):
        A = matrix(ZZ, [[-1, 1, 0, -1, 2, 7], [6, 1, -1, -2, 2, 5]])
        M = ToricArithmeticMatroid(A, ordered_groundset=['a', 'b', 'f', 'e', 'd', 'c'])

        self.assertEqual(M.groundset(), frozenset(['a', 'b', 'f', 'e', 'd', 'c']))
        self.assertEqual(M.full_rank(), M.rank(M.groundset()))

        self.assertEqual(M._Q, matrix(ZZ, 2, 0))
        self.assertEqual(M.rank(['a']), 1)
        self.assertEqual(M.rank(['a', 'b']), 2)
        self.assertEqual(M.rank(['f', 'a', 'b']), 2)

        self.assertEqual(M._multiplicity(['b', 'f']), 1)

        # test failure if there are repetition in the ordered groundset
        with self.assertRaises(AssertionError):
            M = ToricArithmeticMatroid(A, ordered_groundset=['a', 'b', 'f', 'e', 'd', 'b'])

    def test_equivalence(self):
        # without Q
        A = matrix(ZZ, [[1, 1, 2], [0, 7, 7]])
        B = matrix(ZZ, [[-1, -1, 2], [-1, 6, -5]])

        M = ToricArithmeticMatroid(A)
        N = ToricArithmeticMatroid(B)

        self.assertTrue(M.is_equivalent(N))
        self.assertTrue(N.is_equivalent(M))

        self.assertTrue(M.is_isomorphic(N))

        C = matrix(ZZ, [[1, 2, 3], [0, 7, 7]])
        O = ToricArithmeticMatroid(C)
        self.assertFalse(M.is_equivalent(O))
        self.assertFalse(N.is_equivalent(O))
        self.assertFalse(O.is_equivalent(M))
        self.assertFalse(O.is_equivalent(N))

        self.assertTrue(M.is_isomorphic(O))

    def test_non_equivalent_representations(self):
        M = ToricArithmeticMatroid(matrix(ZZ, [[6, 3, -2, 2], [3, 21, 0, -9], [-1, -4, 3, -2]]))

        self.assertEqual(M.num_representations(), 4)
        matroids = [ToricArithmeticMatroid(A) for A in M.all_representations()]
        for N, O in itertools.combinations(matroids, 2):
            self.assertFalse(N.is_equivalent(O))
            self.assertTrue(N.is_isomorphic(O))

    def test_all_representations(self):
        A = matrix(ZZ, [[1, 1, 4, 1, 0], [0, 2, 3, 1, 1], [0, 0, 6, 0, 1]])
        M = ToricArithmeticMatroid(A)
        representations = list(M.all_representations(shnf=False))
        [B] = representations

        shnf_representations = list(M.all_representations())
        [C] = shnf_representations
        self.assertEqual(C, M.representation())

        self.assertEqual(signed_hermite_normal_form(B), C)
        self.assertEqual(signed_hermite_normal_form(A), C)

    def assertMutableCopy(self, R, M):
        """
        Assert that R is a new mutable matrix, and that changing it leaves M as it was.
        """
        twin = copy.deepcopy(M)
        self.assertTrue(R.is_mutable())
        self.assertIsNot(R, M._A)
        self.assertIsNot(R, M._Q)
        if R.nrows() and R.ncols():
            R[0, 0] += 100
        self.assertFalse(M._A.is_mutable())
        self.assertFalse(M._Q.is_mutable())
        self.assertEqual(M, twin)
        self.assertEqual(M.arrangement_matrix(), twin.arrangement_matrix())

    def test_ordered_groundset_is_list_copy(self):
        A = matrix(ZZ, [[-1, 1, 0, -1, 2, 7], [6, 1, -1, -2, 2, 5]])
        E = ['a', 'b', 'f', 'e', 'd', 'c']
        M = ToricArithmeticMatroid(A, ordered_groundset=E)
        L = M.ordered_groundset()
        self.assertIsInstance(L, list)
        self.assertEqual(L, E)
        self.assertIsNot(M.ordered_groundset(), L)  # new list on each call, like Sage's groundset_list()
        L.reverse()
        L.append('z')
        self.assertEqual(M.ordered_groundset(), E)
        self.assertEqual(M.representation(shnf=False), A)

    def test_matrix_accessors_return_mutable_copies(self):
        A = matrix(ZZ, [[-1, 1, 0, -1, 2, 7], [0, 1, -1, -2, 2, 5]])
        Q = matrix(ZZ, [[3, 9, 6], [0, 0, 0]])
        for name, M in (("without Q", ToricArithmeticMatroid(A)), ("with Q", ToricArithmeticMatroid(A, Q))):
            with self.subTest(matroid=name):
                rank_and_multiplicity = [(M.rank(X), M.multiplicity(X)) for X in powerset(M.groundset())]
                for accessor, internal in (("arrangement_matrix", M._A), ("torus_matrix", M._Q)):
                    with self.subTest(accessor=accessor):
                        R = getattr(M, accessor)()
                        self.assertEqual(R, internal)
                        self.assertMutableCopy(R, M)
                # the cache of rank and multiplicity is still valid
                self.assertEqual([(M.rank(X), M.multiplicity(X)) for X in powerset(M.groundset())],
                                 rank_and_multiplicity)

    def test_representations_are_mutable_copies(self):
        A = matrix(ZZ, [[1, 1, 4, 1, 0], [0, 2, 3, 1, 1], [0, 0, 6, 0, 1]])
        M = ToricArithmeticMatroid(A, ordered_groundset=[5, 6, 7, 8, 9])
        for shnf in (True, False):
            for E in (None, [9, 8, 7, 6, 5]):
                with self.subTest(shnf=shnf, ordered_groundset=E):
                    R = M.representation(ordered_groundset=E, shnf=shnf)
                    self.assertEqual(R, M.arithmetic_representation(ordered_groundset=E, shnf=shnf))
                    self.assertMutableCopy(R, M)
        # full multiplicity 1 and 2 take different branches of all_representations()
        for N in (M, ToricArithmeticMatroid(matrix(ZZ, [[1, 1, 0], [0, 2, 2]]))):
            for shnf in (True, False):
                with self.subTest(shnf=shnf, method="all_representations", full_multiplicity=N.full_multiplicity()):
                    expected = list(N.all_representations(shnf=shnf))
                    self.assertTrue(expected)
                    for R in N.all_representations(shnf=shnf):
                        self.assertMutableCopy(R, N)
                    self.assertEqual(list(N.all_representations(shnf=shnf)), expected)

    def test_representation_uses_toric_override(self):
        # representation() goes through arithmetic_representation(), so Toric's version (which keeps A) is used
        A = matrix(ZZ, [[1, 1, 4, 1, 0], [0, 2, 3, 1, 1], [0, 0, 6, 0, 1]])
        M = ToricArithmeticMatroid(A)
        self.assertEqual(M.representation(shnf=False), A)
        self.assertEqual(M.arithmetic_representation(shnf=False), A)
        self.assertEqual(M.representation(), signed_hermite_normal_form(A))

    def test_decomposition(self):
        M = ToricArithmeticMatroid(matrix(ZZ, [[1, 5, -1, 5], [1, 7, 2, 6], [0, -2, -2, -1]]))
        self.assertEqual(M.decomposition(), SetPartition([[2], [0, 1, 3]]))
        self.assertTrue(M.is_decomposable())

    def test_decomposition2(self):
        M = ToricArithmeticMatroid(matrix(ZZ, [[1, 5, -10, 50], [0, 0, 2, 6], [0, 0, -2, -1]]))
        self.assertEqual(M.decomposition(), SetPartition([[0, 1], [2, 3]]))
        self.assertTrue(M.is_decomposable())

    def test_decomposition3(self):
        M = ToricArithmeticMatroid(matrix(ZZ, [[1, 1], [0, 2]]))
        self.assertEqual(M.decomposition(), SetPartition([[0, 1]]))
        self.assertFalse(M.is_decomposable())

    def test_poset_of_layers(self):
        M = ToricArithmeticMatroid(matrix(ZZ, [[1, 0, 1], [0, 1, 3]]))
        P = M.poset_of_layers()
        self.assertEqual(len(P), 7)
        self.assertTrue(P.has_bottom())
        self.assertFalse(P.has_top())
        self.assertEqual(len(P.cover_relations()), 10)
        self.assertTrue(P.is_isomorphic(
            Poset(data=(range(7), [(0, 1), (0, 2), (0, 3), (1, 4), (2, 4), (2, 5), (2, 6), (3, 4), (3, 5), (3, 6)]))))

    def test_poset_of_layers_relabeled(self):
        M = ToricArithmeticMatroid(matrix(ZZ, [[1, 0, 1], [0, 1, 3]]), ordered_groundset=['a', 'c', 'b'])
        P = M.poset_of_layers()
        self.assertEqual(len(P), 7)
        self.assertTrue(P.has_bottom())
        self.assertFalse(P.has_top())
        self.assertEqual(len(P.cover_relations()), 10)
        self.assertTrue(P.is_isomorphic(
            Poset(data=(range(7), [(0, 1), (0, 2), (0, 3), (1, 4), (2, 4), (2, 5), (2, 6), (3, 4), (3, 5), (3, 6)]))))
        self.assertEqual(set(S for (S, x) in P), {(), ('a',), ('b',), ('c',), ('a', 'c', 'b'), ('a', 'b')})

    def test_poset_of_layers_root_system_C2(self):
        # see [DGP17]
        A = matrix(ZZ, [[2, 0, 1, 1, ], [0, 2, 1, -1]])
        M = ToricArithmeticMatroid(A)
        P = M.poset_of_layers()
        self.assertEqual(len(P), 11)
        self.assertEqual(len(P.cover_relations()), 18)
        homology = P.subposet([a for a in P if a != P.bottom()]).order_complex(on_ints=True).homology()
        self.assertEqual([homology[i].ngens() for i in range(2)], [0, 3])
        self.assertEqual(homology[1].ngens(), M.arithmetic_tutte_polynomial()(x=0, y=0))

    def test_poset_of_layers_root_system_C3(self):
        # see [DGP17]
        A = matrix(ZZ, [[2, 0, 0, 1, 1, 1, 1, 0, 0], [0, 2, 0, 1, -1, 0, 0, 1, 1], [0, 0, 2, 0, 0, 1, -1, 1, -1]])
        M = ToricArithmeticMatroid(A)
        P = M.poset_of_layers()
        homology = P.subposet([a for a in P if a != P.bottom()]).order_complex(on_ints=True).homology()
        self.assertEqual([homology[i].ngens() for i in range(3)], [0, 0, 15])
        self.assertEqual(homology[2].ngens(), M.arithmetic_tutte_polynomial()(x=0, y=0))

    def test_poset_of_layers_root_system_B2(self):
        # see [DGP17]
        A = matrix(ZZ, [[1, 0, 1, 1, ], [0, 1, 1, -1]])
        M = ToricArithmeticMatroid(A)
        P = M.poset_of_layers()
        homology = P.subposet([a for a in P if a != P.bottom()]).order_complex(on_ints=True).homology()
        self.assertEqual([homology[i].ngens() for i in range(2)], [0, 1])
        self.assertEqual(homology[1].ngens(), M.arithmetic_tutte_polynomial()(x=0, y=0))

    def test_poset_of_layers_root_system_B3(self):
        # see [DGP17]
        A = matrix(ZZ, [[1, 0, 0, 1, 1, 1, 1, 0, 0], [0, 1, 0, 1, -1, 0, 0, 1, 1], [0, 0, 1, 0, 0, 1, -1, 1, -1]])
        M = ToricArithmeticMatroid(A)
        P = M.poset_of_layers()
        homology = P.subposet([a for a in P if a != P.bottom()]).order_complex(on_ints=True).homology()
        self.assertEqual([homology[i].ngens() for i in range(3)], [0, 0, 6])
        self.assertEqual(homology[2].ngens(), M.arithmetic_tutte_polynomial()(x=0, y=0))

    def test_poset_of_layers_root_system_D3(self):
        # see [DGP17]
        A = matrix(ZZ, [[1, 1, 1, 1, 0, 0], [1, -1, 0, 0, 1, 1], [0, 0, 1, -1, 1, -1]])
        M = ToricArithmeticMatroid(A)
        P = M.poset_of_layers()
        homology = P.subposet([a for a in P if a != P.bottom()]).order_complex(on_ints=True).homology()
        self.assertEqual([homology[i].ngens() for i in range(3)], [0, 0, 3])
        self.assertEqual(homology[2].ngens(), M.arithmetic_tutte_polynomial()(x=0, y=0))

    def test_different_posets_of_layers(self):
        # see [PP19, Section 8]
        A = matrix(ZZ, [[1, 1, 1, 3], [0, 5, 0, 5], [0, 0, 5, 5]])
        B = matrix(ZZ, [[1, 4, 1, 6], [0, 5, 0, 5], [0, 0, 5, 5]])
        C = matrix(ZZ, [[1, 2, 2, 1], [0, 5, 0, 5], [0, 0, 5, -5]])

        # check signed Hermite normal form
        self.assertEqual(signed_hermite_normal_form(A), matrix(ZZ, [[1, 1, 1, -3], [0, 5, 0, -5], [0, 0, 5, -5]]))
        self.assertEqual(signed_hermite_normal_form(B), matrix(ZZ, [[1, 1, 1, -1], [0, 5, 0, 5], [0, 0, 5, -5]]))
        self.assertEqual(signed_hermite_normal_form(C), matrix(ZZ, [[1, 2, 2, 1], [0, 5, 0, 5], [0, 0, 5, -5]]))

        MA = ToricArithmeticMatroid(A)
        MB = ToricArithmeticMatroid(B)
        MC = ToricArithmeticMatroid(C)

        self.assertTrue(MA.is_isomorphic(MB))
        self.assertTrue(MA.is_isomorphic(MC))
        self.assertTrue(MB.is_isomorphic(MC))

        self.assertFalse(MA.poset_of_layers().is_isomorphic(MB.poset_of_layers()))
        self.assertFalse(MA.poset_of_layers().is_isomorphic(MC.poset_of_layers()))
        self.assertFalse(MB.poset_of_layers().is_isomorphic(MC.poset_of_layers()))

    def test_independence_poset(self):
        A = matrix(ZZ, [[1, 0, 1], [0, 1, 3]])
        M = ToricArithmeticMatroid(A)
        P = M.arithmetic_independence_poset()

        self.assertEqual(len(P), 9)
        self.assertEqual(len(P.cover_relations()), 13)
        self.assertEqual(P.rank(), 2)
        self.assertTrue(P.has_bottom())

    def test_independence_poset_relabeled(self):
        M = ToricArithmeticMatroid(matrix(ZZ, [[1, 0, 1], [0, 1, 3]]), ordered_groundset=['a', 'c', 'b'])
        P = M.arithmetic_independence_poset()

        self.assertEqual(len(P), 9)
        self.assertEqual(len(P.cover_relations()), 13)
        self.assertEqual(P.rank(), 2)
        self.assertTrue(P.has_bottom())
        self.assertEqual(set(S for (S, x) in P), {(), ('a',), ('b',), ('c',), ('a', 'c'), ('a', 'b'), ('c', 'b')})

    def test_nonshellable_poset_of_layers(self):
        # see [PP19, Section 8]
        A = matrix(ZZ, [[1, 1, 1, -3], [0, 5, 0, -5], [0, 0, 5, -5]])
        B = matrix(ZZ, [[1, 1, 1, -1], [0, 5, 0, 5], [0, 0, 5, -5]])
        C = matrix(ZZ, [[1, 2, 2, 1], [0, 5, 0, 5], [0, 0, 5, -5]])

        for X in [A, B, C]:
            M = ToricArithmeticMatroid(X)
            P = M.poset_of_layers()
            homology = P.subposet([a for a in P if a != P.bottom()]).order_complex(on_ints=True).homology()

            from sage.homology.homology_group import HomologyGroup
            self.assertEqual(homology, {
                0: HomologyGroup(0, ZZ),  # trivial group
                1: HomologyGroup(1, ZZ, [5]),  # Z/5Z
                2: HomologyGroup(48, ZZ),  # Z^48
            })

    def test_nonshellable_arithmetic_independence_poset(self):
        # see [PP19, Section 8]
        A = matrix(ZZ, [[1, 1, 1, -3], [0, 5, 0, -5], [0, 0, 5, -5]])
        B = matrix(ZZ, [[1, 1, 1, -1], [0, 5, 0, 5], [0, 0, 5, -5]])
        C = matrix(ZZ, [[1, 2, 2, 1], [0, 5, 0, 5], [0, 0, 5, -5]])

        for X in [A, B, C]:
            M = ToricArithmeticMatroid(X)
            P = M.arithmetic_independence_poset()
            homology = P.subposet([a for a in P if a != P.bottom()]).order_complex(on_ints=True).homology()

            from sage.homology.homology_group import HomologyGroup
            self.assertEqual(homology, {
                0: HomologyGroup(0, ZZ),  # trivial group
                1: HomologyGroup(1, ZZ, [5]),  # Z/5Z
                2: HomologyGroup(73, ZZ),  # Z^73
            })
