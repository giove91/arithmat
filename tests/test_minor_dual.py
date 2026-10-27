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
import operator
import unittest
from functools import reduce

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


class _FreshCopyBasisMatroid(BasisMatroid):
    # copy by returning a new plain BasisMatroid instead of self, as Sage may
    def __copy__(self):
        return BasisMatroid(groundset=self.groundset(), bases=self.bases())

    def __deepcopy__(self, memo=None):
        return BasisMatroid(groundset=self.groundset(), bases=self.bases())


class FreshCopyBasisArithmeticMatroid(BasisArithmeticMatroid, _FreshCopyBasisMatroid):
    def _arithmetic_copy(self, matroid):
        return FreshCopyBasisArithmeticMatroid(matroid)


class TestDualAndMinor(unittest.TestCase):

    def setUp(self):
        # 0 and 2 are loops with torsion: m(empty) = 2, m({2}) = 1
        self.T = ToricArithmeticMatroid(matrix(ZZ, [[1, 0, 1, 1], [0, 1, 0, 1]]), matrix(ZZ, [[2], [0]]))
        # 2 is a coloop with m({2}) = 2
        self.A_lin = matrix(ZZ, [[1, 1, 0], [0, 0, 2]])

    def matroids(self):
        T = self.T
        A = self.A_lin

        def m_lin(X):
            return reduce(operator.mul, [d for d in A[:, sorted(X)].elementary_divisors() if d != 0], 1) if X else 1

        # U(2,5) on {1, ..., 5}, with m = 2 on pairs inside {3, 4, 5}
        def rk_u25(X):
            return min(2, len(X))

        def m_u25(X):
            return 2 if len(X) == 2 and all(x in [3, 4, 5] for x in X) else 1

        return {
            "rank": ArithmeticMatroid(T.groundset(), T._rank, T._multiplicity),
            "U(2,5)": ArithmeticMatroid([1, 2, 3, 4, 5], rk_u25, m_u25),
            "basis": BasisArithmeticMatroid(BasisMatroid(groundset=T.groundset(), bases=T.bases()),
                                            multiplicity_function=T._multiplicity),
            "linear": LinearArithmeticMatroid(A, multiplicity_function=m_lin),
            "toric": T,
        }

    def assertMinorOf(self, M, N, C, D):
        """
        Assert that M is N / C \\ D.
        """
        self.assertEqual(M.groundset(), N.groundset() - C - D)
        for X in powerset(M.groundset()):
            self.assertEqual(M.rank(X), N.rank(C.union(X)) - N.rank(C), (sorted(C), sorted(D), X))
            self.assertEqual(M.multiplicity(X), N.multiplicity(C.union(X)), (sorted(C), sorted(D), X))

    def assertDualOf(self, M, N):
        """
        Assert that M is the dual of N.
        """
        E = N.groundset()
        self.assertEqual(M.groundset(), E)
        for X in powerset(E):
            self.assertEqual(M.rank(X), len(X) + N.rank(E.difference(X)) - N.full_rank(), X)
            self.assertEqual(M.multiplicity(X), N.multiplicity(E.difference(X)), X)

    def assertSameArithmeticMatroid(self, M, N):
        self.assertEqual(M.groundset(), N.groundset())
        for X in powerset(N.groundset()):
            self.assertEqual(M.rank(X), N.rank(X), X)
            self.assertEqual(M.multiplicity(X), N.multiplicity(X), X)

    def assertValid(self, M):
        # Skip rank-1 linear matroids over a non-field ring (e.g. ZZ). Sage's LinearMatroid.is_valid() checks
        # these via cross_ratios(), which raises on rank 1 (min() of an empty set in _line_ratios, Sage 10.10),
        # so is_valid() wrongly returns False. TODO: restore once fixed upstream.
        if isinstance(M, LinearMatroid) and not M.base_ring().is_field() and M.full_rank() == 1:
            return
        self.assertTrue(M.is_valid())

    def test_minor(self):
        for name, N in self.matroids().items():
            # if N.dual() fails, skip this matroid rather than failing every (C, D) below
            with self.subTest(matroid=name):
                N_dual = N.dual()
                for C in map(frozenset, powerset(sorted(N.groundset()))):
                    for D in map(frozenset, powerset(sorted(N.groundset() - C))):
                        with self.subTest(contractions=sorted(C), deletions=sorted(D)):
                            M = N.minor(C, D)
                            self.assertMinorOf(M, N, C, D)
                            self.assertMinorOf(N.minor(contractions=C).minor(deletions=D), N, C, D)
                            self.assertMinorOf(M.dual().dual(), N, C, D)
                            self.assertValid(M)
                            self.assertValid(M.dual())
                            self.assertDualOf(M.dual(), M)
                            # (N / C \ D)* = N* / D \ C
                            self.assertSameArithmeticMatroid(M.dual(), N_dual.minor(D, C))

    def test_minor_and_dual_classes(self):
        Minor, Dual = MinorArithmeticMatroid, DualArithmeticMatroid
        # (minor, dual, minor of dual)
        expected = {
            "rank": (Minor, Dual, Dual),
            "U(2,5)": (Minor, Dual, Dual),
            "basis": (BasisArithmeticMatroid,) * 3,
            "linear": (LinearArithmeticMatroid,) * 3,
            "toric": (ToricArithmeticMatroid,) * 3,
        }
        for name, N in self.matroids().items():
            e = min(N.groundset())
            minor_cls, dual_cls, dual_minor_cls = expected[name]
            with self.subTest(matroid=name):
                self.assertIsInstance(N.minor(contractions=[e]), minor_cls)
                self.assertIsInstance(N.dual(), dual_cls)
                self.assertIsInstance(N.dual().minor(contractions=[e]), dual_minor_cls)

    def test_dual(self):
        for name, N in self.matroids().items():
            with self.subTest(matroid=name):
                self.assertDualOf(N.dual(), N)

    def test_minor_argument_forms(self):
        for name, N in self.matroids().items():
            e = max(N.groundset())
            forms = {"element": lambda: e, "list": lambda: [e], "tuple": lambda: (e,),
                     "frozenset": lambda: frozenset([e]), "generator": lambda: (x for x in [e])}
            for form, arg in forms.items():
                with self.subTest(matroid=name, arg=form):
                    self.assertMinorOf(N.minor(contractions=arg()), N, frozenset([e]), frozenset())
                    self.assertMinorOf(N.minor(deletions=arg()), N, frozenset(), frozenset([e]))
            if 0 in N.groundset():
                with self.subTest(matroid=name, arg=0):
                    # 0 is falsy but is an element
                    self.assertMinorOf(N.minor(contractions=0), N, frozenset([0]), frozenset())

    def test_minor_invalid_arguments(self):
        for name, N in self.matroids().items():
            e = max(N.groundset())
            with self.subTest(matroid=name, arg="not in groundset"):
                self.assertRaises(ValueError, N.minor, contractions=[9])
            with self.subTest(matroid=name, arg="not disjoint"):
                self.assertRaises(ValueError, N.minor, contractions=[e], deletions=[e])

    def test_minor_repr_and_equality(self):
        N = self.matroids()["rank"]
        M1 = N.minor(contractions=[2]).minor(deletions=[3])
        M2 = N.minor(contractions=[2], deletions=[3])
        self.assertEqual(M1, M2)
        self.assertNotEqual(N.minor(contractions=[2]), N.minor(deletions=[2]))  # same matroid, not the same arithmetic matroid
        self.assertIn("/ {2}", repr(M1))
        self.assertIn("\\ {3}", repr(M1))

    def test_copy(self):
        N = self.matroids()["rank"]
        cases = dict(self.matroids())
        cases["minor"] = N.minor(contractions=[2], deletions=[3])
        cases["dual"] = N.dual()
        for name, M in cases.items():
            for how in (copy.copy, copy.deepcopy):
                with self.subTest(matroid=name, copy=how.__name__):
                    # Sage matroids are immutable, so Sage's copy may return M itself
                    M2 = how(M)
                    self.assertSameArithmeticMatroid(M2, M)
                    if name in ("minor", "dual"):
                        self.assertEqual(M2, M)
                        self.assertEqual(repr(M2), repr(M))
                        self.assertNotIn('_multiplicity', vars(M2))  # uses its own method, not M's

    def test_copy_rank_matroid_equality(self):
        # Regression: copies of rank-based matroids could compare unequal to the original.
        # "rank" has a bound method as rank function, and a deep copy of it is a different method.
        for name in ("rank", "U(2,5)"):
            M = self.matroids()[name]
            for how in (copy.copy, copy.deepcopy):
                with self.subTest(matroid=name, copy=how.__name__):
                    M2 = how(M)
                    self.assertIsNot(M2, M)  # our copy ran, not Sage returning M itself
                    self.assertIsInstance(M2, ArithmeticMatroid)
                    self.assertEqual(M2, M)
                    self.assertEqual(hash(M2), hash(M))
                    self.assertIs(M2._rank_function, M._rank_function)
                    self.assertIs(M2._multiplicity, M._multiplicity)
            with self.subTest(matroid=name, copy="fresh plain copy"):
                # when the parent's copy returns a plain RankMatroid with the same rank function
                plain = RankMatroid(M.groundset(), M._rank_function)
                self.assertEqual(M._with_multiplicity(plain, M._multiplicity), M)

    def test_copy_keeps_multiplicity_function(self):
        # Regression: deepcopy replaced the multiplicity function with a deep copy
        def make():
            return self.matroids()["basis"]
        for how in (copy.copy, copy.deepcopy):
            for wrap in ("plain", "minor", "dual"):
                with self.subTest(copy=how.__name__, wrap=wrap):
                    M, twin = make(), make()
                    if wrap == "minor":
                        M, twin = MinorArithmeticMatroid(M, [2]), MinorArithmeticMatroid(twin, [2])
                    elif wrap == "dual":
                        M, twin = DualArithmeticMatroid(M), DualArithmeticMatroid(twin)
                    M2 = how(M)
                    self.assertEqual(M2, M)
                    self.assertEqual(M, twin)  # the original is unchanged
                    self.assertSameArithmeticMatroid(M2, M)
            with self.subTest(copy=how.__name__, wrap="fresh plain copy"):
                N = make()
                M = FreshCopyBasisArithmeticMatroid(N, multiplicity_function=N._multiplicity)
                M2 = how(M)
                self.assertIsNot(M2, M)
                self.assertEqual(M2, M)
                self.assertSameArithmeticMatroid(M2, M)

    def test_arithmetic_copy(self):
        # Sage currently returns arithmetic classes here, so call the fallback directly on plain matroids
        for name, N in self.matroids().items():
            if name == "toric":
                continue
            with self.subTest(matroid=name):
                plain = {
                    "rank": lambda: RankMatroid(N.groundset(), N.rank),
                    "U(2,5)": lambda: RankMatroid(N.groundset(), N.rank),
                    "basis": lambda: BasisMatroid(groundset=N.groundset(), bases=N.bases()),
                    "linear": lambda: LinearMatroid(matrix=N.representation(), groundset=N.groundset_list()),
                }[name]()
                M = N._with_multiplicity(plain, N._multiplicity)
                self.assertIsInstance(M, type(N))
                self.assertSameArithmeticMatroid(M, N)

    def test_dual_arithmetic_matroid(self):
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
        self.assertTrue(M.is_valid())

        M1 = M.dual()
        self.assertIsInstance(M1, DualArithmeticMatroid)
        self.assertTrue(M1.is_valid())
        self.assertNotEqual(M1, M)

        M2 = M1.dual()
        self.assertIsInstance(M2, ArithmeticMatroid)
        self.assertNotIsInstance(M2, DualArithmeticMatroid)
        self.assertEqual(M, M2)

        M3 = M1.minor()
        self.assertIsInstance(M3, DualArithmeticMatroid)
        # self.assertEqual(M1, M3) # this is not an equality for (non-arithmetic) matroids
        self.assertTrue(M3.is_valid())

    def test_contract_delete(self):
        A = matrix(ZZ, [[2, 0, 4], [0, 3, 0]])
        M = ToricArithmeticMatroid(A)
        D = M.delete([1, 2])
        C = M.contract([1, 2])
        self.assertEqual(C.full_multiplicity(), 6)
        self.assertEqual(D.full_multiplicity(), 2)
        self.assertEqual(C.multiplicity([]), 12)
        self.assertEqual(D.multiplicity([]), 1)
        self.assertEqual(C, M.minor(contractions=[1, 2], deletions=[]))
        self.assertEqual(D, M.minor(contractions=[], deletions=[1, 2]))

    def test_dual_linear_matroid(self):
        A = matrix(QQ, [[-1, 1, 0, -1, 2, 7], [6, 1, -1, -2, 2, 5]])

        def m(X):
            return 1

        M = LinearArithmeticMatroid(A, multiplicity_function=m)
        self.assertTrue(M.is_valid())

        M1 = M.dual()
        self.assertTrue(M1.is_valid())

        M2 = M1.dual()
        self.assertTrue(M2.is_valid())
        self.assertNotEqual(M, M2)  # multiplicity functions are not equal in the sense of ==
        self.assertTrue(M.equals(M2))
        self.assertNotEqual(M, M1)
        self.assertNotEqual(M1, M2)

    def test_minor_linear_matroid(self):
        A = matrix(QQ, [[-1, 1, 0, -1, 2, 7], [6, 1, -1, -2, 2, 5]])

        def m(X):
            return 1

        M = LinearArithmeticMatroid(A, multiplicity_function=m)
        self.assertTrue(M.is_valid())

        M1 = M.minor(contractions=frozenset([1]), deletions=frozenset([2]))
        self.assertIsInstance(M1, LinearArithmeticMatroid)
        self.assertTrue(M1.is_valid())

        M1a = M.contract([1]).delete([2])
        self.assertNotEqual(M1, M1a)  # multiplicity functions are not equal in the sense of ==
        self.assertTrue(M1.equals(M1a))

        M2 = M1.dual()
        self.assertTrue(M2.is_valid())
        self.assertIsInstance(M2, LinearArithmeticMatroid)
        self.assertNotEqual(M, M2)
        self.assertFalse(M.equals(M2))

        N1 = M.dual()
        self.assertIsInstance(N1, LinearArithmeticMatroid)
        self.assertTrue(N1.is_valid())

        N2 = N1.minor(contractions=frozenset([2]), deletions=frozenset([1]))
        self.assertTrue(N2.is_valid())
        self.assertIsInstance(N2, LinearArithmeticMatroid)
        self.assertNotEqual(M2, N2)  # multiplicity functions are not equal in the sense of ==
        self.assertTrue(M2.equals(N2))

        N3 = N1.minor(contractions=frozenset([1]), deletions=frozenset([2]))
        self.assertTrue(N3.is_valid())
        self.assertIsInstance(N3, LinearArithmeticMatroid)
        self.assertNotEqual(M2, N3)
        self.assertFalse(M2.equals(N3))
