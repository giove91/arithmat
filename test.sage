# -*- coding: utf-8 -*-
"""
Tests for Arithmat.

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

import unittest
import itertools
from arithmat.arithmetic_matroid import *
from arithmat.shnf import signed_hermite_normal_form as shnf


def representation_to_matroid(A):
    """
    Given an integer matrix A, return the associated ArithmeticMatroid.
    """
    E = range(A.ncols())

    def rk(X):
        return A[:,list(X)].rank()

    def m(X):
        return reduce(operator.mul, [d for d in A[:,list(X)].elementary_divisors() if d != 0], 1)

    return ArithmeticMatroid(E, rk, m)



class TestArithmeticMatroid(unittest.TestCase):

    def test_trivial(self):
        E = range(1, 10)

        def rk(X):
            return len(X)

        def m(X):
            return 1

        M = ArithmeticMatroid(E, rk, m)
        self.assertTrue(M.is_valid())
        self.assertEqual(M.full_rank(), 9)
        self.assertTrue(M.is_representable())


    def test_valid(self):
        # matroid with representation
        # [ 1, a, 1 ]
        # [ 0, b, 5 ]
        # with a|b and (b,5)=1

        E = [1,2,3]

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
                (1,2): b,
                (1,3): 5,
                (2,3): 5*a-b,
                (1,2,3): 1,
            }[X]

        M = ArithmeticMatroid(E, rk, m)
        self.assertTrue(M.is_valid())
        self.assertTrue(M.is_representable())
        self.assertEqual(M.full_multiplicity(), 1)


    def test_not_valid(self):
        E = [1,2,3]

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
                (1,2): b,
                (1,3): 5,
                (2,3): 5*a-b+1,
                (1,2,3): 1,
            }[X]

        M = ArithmeticMatroid(E, rk, m)
        self.assertFalse(M.is_valid())


    def test_valid2(self):
        # valid, not representable, not orientable
        E = [1,2,3,4,5]

        def rk(X):
            return min(2, len(X))

        def m(X):
            if len(X) == 2 and all(x in [3,4,5] for x in X):
                return 2
            else:
                return 1

        M = ArithmeticMatroid(E, rk, m)
        self.assertTrue(M.is_valid())
        self.assertFalse(M.is_representable())
        self.assertFalse(M.is_orientable())
        self.assertEqual(M.representation(), None)
        self.assertEqual(M.full_multiplicity(), 1)


    def test_relizable(self):
        # representable with a 2x2 matrix
        for a in range(6):
            for b in range(1, 6):
                A = matrix(ZZ, [[1, a], [0, b]])
                M = representation_to_matroid(A)

                self.assertTrue(M.is_valid())
                self.assertTrue(M.is_representable())
                self.assertTrue(M.is_orientable())

                R = ZZ['x, y']
                x, y = R._first_ngens(2)
                self.assertEqual(M.arithmetic_tutte_polynomial(), x**2 + x*(gcd(a,b)-1) + b-gcd(a,b))


    def test_representable2(self):
        A = matrix(ZZ, [[-1, -29, -1, 1], [1, -1, 0, 1]])
        M = representation_to_matroid(A)

        self.assertTrue(M.is_valid())
        self.assertTrue(M.is_representable())
        self.assertTrue(M.is_orientable())


    def test_representable3(self):
        A = matrix(ZZ, [[-1, 1, 0, 0, -1], [6, 1, 1, -1, -1]])
        M = representation_to_matroid(A)

        self.assertTrue(M.is_valid())
        self.assertTrue(M.is_representable())
        self.assertTrue(M.is_orientable())


    def test_representable4(self):
        A = matrix(ZZ, [[2, 2, 1, 0, 0], [1, 5, -1, 1, -2], [-2, 1, 0, -1, -1]])
        M = representation_to_matroid(A)

        self.assertTrue(M.is_valid())
        self.assertTrue(M.is_representable())
        self.assertTrue(M.is_orientable())


    def test_representable_random(self):
        A = random_matrix(ZZ,4,6)
        M = representation_to_matroid(A)

        self.assertTrue(M.is_valid())
        self.assertTrue(M.is_representable())
        self.assertTrue(M.is_orientable())


    def test_pseudo(self):
        E = [1,2,3]

        def rk(X):
            return min(2, len(X))

        def m(X):
            if len(X) == 2:
                return 3
            else:
                return 1

        M = ArithmeticMatroid(E, rk, m)
        self.assertTrue(M.is_valid())
        self.assertFalse(M.is_representable())
        self.assertTrue(M.is_orientable())


    def test_non_representable(self):
        A = matrix(ZZ, [[-1, 1, 0, -1], [6, 1, -1, -2]])
        M = representation_to_matroid(A)
        M2 = ArithmeticMatroid(M.groundset(), M.rank, lambda X: M._multiplicity(X)**2)

        self.assertTrue(M2.is_valid())
        self.assertTrue(M.is_representable())
        self.assertFalse(M2.is_representable())


    def test_non_representable2(self):
        A = matrix(ZZ, [[-1, 1, 0, -1, 2, 7], [6, 1, -1, -2, 2, 5]])
        M = representation_to_matroid(A)
        M2 = ArithmeticMatroid(M.groundset(), M.rank, lambda X: M._multiplicity(X)**2)

        self.assertTrue(M2.is_valid())
        self.assertTrue(M.is_representable())
        self.assertFalse(M2.is_representable())


    def test_non_orientable(self):
        E = range(1, 7)

        def rk(A):
            return min(len(A), 2)

        def m(X):
            A = tuple(sorted(X))
            if A == (1,2,3):
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
            elif A == (1,2):
                return 4
            else:
                return 1

        M = ArithmeticMatroid(E, rk, m)
        self.assertTrue(M.is_valid())
        self.assertFalse(M.is_representable())
        self.assertFalse(M.is_orientable())


    def test_non_essential(self):
        A = matrix(ZZ, [[-1,  1,  0], [ 6,  1, -1], [2, -3, 0], [1, 2, 3], [-1, 0, 0]])
        M = representation_to_matroid(A)

        self.assertTrue(M.is_valid())
        self.assertTrue(M.is_representable())


    def test_hash(self):
        E = [1,2,3,4,5]

        def rk(X):
            return min(2, len(X))

        def m(X):
            if len(X) == 2 and all(x in [3,4,5] for x in X):
                return 2
            else:
                return 1

        M = ArithmeticMatroid(E, rk, m)
        s = set([M])
        self.assertEqual(hash(M), hash((frozenset(E), 2, 1, 1)))
        self.assertEqual(M.groundset(), frozenset(E))


    def test_isomorphism(self):
        E = [1,2,3,4,5]

        def rk(X):
            return min(2, len(X))

        def m(X):
            if len(X) == 2 and all(x in [3,4,5] for x in X):
                return 2
            else:
                return 1

        M = ArithmeticMatroid(E, rk, m)
        self.assertTrue(M.is_valid())
        self.assertFalse(M.is_representable())
        self.assertFalse(M.is_orientable())

        E1 = [4,5,6,7,8]

        def m1(X):
            if len(X) == 2 and all(x in [6,7,8] for x in X):
                return 2
            else:
                return 1

        M1 = ArithmeticMatroid(E1, rk, m1)
        self.assertTrue(M1.is_valid())
        self.assertFalse(M1.is_representable())
        self.assertFalse(M1.is_orientable())

        self.assertTrue(M.is_isomorphism(M1, {i: i+3 for i in E}))
        self.assertFalse(M.equals(M1))

        M2 = ArithmeticMatroid(copy.copy(E), rk, m)
        self.assertTrue(M.equals(M2))
        self.assertTrue(M2.equals(M))


    def test_isomorphic(self):
        A = matrix(ZZ, [[1]])
        B = matrix(ZZ, [[3,1],[0,1]])
        M = ToricArithmeticMatroid(A)
        DN = ToricArithmeticMatroid(B).minor(deletions = [0])
        CN = ToricArithmeticMatroid(B).minor(contractions = [0])
        self.assertTrue(DN.is_isomorphic(M))
        self.assertFalse(CN.is_isomorphic(M))

    def test_num_representations(self):
        r = 3

        for m in range(23, 26): # m(E)
            E = range(r)
            A = matrix(ZZ, r, r, lambda i, j: 1 if i == j and i < r-1 else 0 if j < r-1 else m if i == r-1 else 1)

            M = representation_to_matroid(A)
            self.assertEqual(M.num_representations(), (euler_phi(m)//2)**(r-1))


    def test_gcd(self):
        M = ToricArithmeticMatroid(matrix(ZZ, [[6, 3, -2, 2], [3, 21, 0, -9], [-1, -4, 3, -2]]))

        self.assertTrue(M.is_torsion_free())
        self.assertFalse(M.is_surjective())
        self.assertTrue(M.is_representable())
        self.assertTrue(M.is_gcd())
        self.assertFalse(M.is_strong_gcd())

        M1 = M.reduction()

        self.assertTrue(M1.is_torsion_free())
        self.assertTrue(M1.is_surjective())
        self.assertTrue(M1.is_representable())
        self.assertTrue(M1.is_gcd())

        self.assertTrue(M1.is_strong_gcd())



class TestDualAndMinor(unittest.TestCase):
    # TODO: test delete and contract, test multiplicity

    def test_dual_arithmetic_matroid(self):
        # test of DualArithmeticMatroid

        E = [1,2,3,4,5]

        def rk(X):
            return min(2, len(X))

        def m(X):
            if len(X) == 2 and all(x in [3,4,5] for x in X):
                return 2
            else:
                return 1

        M = ArithmeticMatroid(E, rk, m)
        self.assertTrue(M.is_valid())

        M1 = M.dual()
        self.assertIsInstance(M1, DualArithmeticMatroid)
        self.assertTrue(M1.is_valid())

        M2 = M1.dual()
        self.assertIsInstance(M2, ArithmeticMatroid)
        self.assertNotIsInstance(M2, DualArithmeticMatroid)
        self.assertEqual(M, M2)

        M3 = M1.minor()
        self.assertIsInstance(M3, DualArithmeticMatroid)
        # self.assertEqual(M1, M3) # this is not an equality for (non-arithmetic) matroids
        self.assertTrue(M3.is_valid())


    def test_minor_arithmetic_matroid(self):
        # test of MinorArithmeticMatroid

        E = [1,2,3,4,5]

        def rk(X):
            return min(2, len(X))

        def m(X):
            if len(X) == 2 and all(x in [3,4,5] for x in X):
                return 2
            else:
                return 1

        M = ArithmeticMatroid(E, rk, m)
        self.assertTrue(M.is_valid())

        M1 = M.minor(contractions=[1], deletions=[2])
        self.assertIsInstance(M1, MinorArithmeticMatroid)
        self.assertTrue(M1.is_valid())

        M2 = M1.dual()
        self.assertTrue(M2.is_valid())
        self.assertIsInstance(M2, DualArithmeticMatroid)
        self.assertNotEqual(M, M2)


        N1 = M.dual()
        self.assertIsInstance(N1, DualArithmeticMatroid)
        self.assertTrue(N1.is_valid())

        N2 = N1.minor(contractions=[2], deletions=[1])
        self.assertTrue(N2.is_valid())
        self.assertIsInstance(N2, DualArithmeticMatroid)
        self.assertEqual(M2, N2)

        N3 = N1.minor(contractions=[1], deletions=[2])
        self.assertTrue(N3.is_valid())
        self.assertIsInstance(N3, DualArithmeticMatroid)
        self.assertNotEqual(M2, N3)


    def test_contract_delete(self):
        A = matrix(ZZ, [[2,0,4],[0,3,0]])
        M = ToricArithmeticMatroid(A)
        D = M.delete([1,2])
        C = M.contract([1,2])
        self.assertEqual(C.full_multiplicity(), 6)
        self.assertEqual(D.full_multiplicity(), 2)
        self.assertEqual(C.multiplicity([]), 12)
        self.assertEqual(D.multiplicity([]), 1)
        self.assertEqual(C,M.minor(contractions=[1,2], deletions=[]))
        self.assertEqual(D,M.minor(contractions=[], deletions=[1,2]))


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
        self.assertNotEqual(M, M2) # multiplicity functions are not equal in the sense of ==
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
        self.assertNotEqual(M1, M1a) # multiplicity functions are not equal in the sense of ==
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
        self.assertNotEqual(M2, N2) # multiplicity functions are not equal in the sense of ==
        self.assertTrue(M2.equals(N2))

        N3 = N1.minor(contractions=frozenset([1]), deletions=frozenset([2]))
        self.assertTrue(N3.is_valid())
        self.assertIsInstance(N3, LinearArithmeticMatroid)
        self.assertNotEqual(M2, N3)
        self.assertFalse(M2.equals(N3))


    def test_delete(self):
        E = [1,2,3,4,5]

        def rk(X):
            return min(2, len(X))

        def m(X):
            if len(X) == 2 and all(x in [3,4,5] for x in X):
                return 2
            else:
                return 1

        M = ArithmeticMatroid(E, rk, m)
        N = M.delete([2,3])
        self.assertEqual(N.groundset(), frozenset([1,4,5]))

        N = N.delete([5])
        self.assertEqual(N.groundset(), frozenset([1,4]))
        self.assertFalse(repr(N) is None)


    def test_contract(self):
        E = [1,2,3,4,5]

        def rk(X):
            return min(2, len(X))

        def m(X):
            if len(X) == 2 and all(x in [3,4,5] for x in X):
                return 2
            else:
                return 1

        M = ArithmeticMatroid(E, rk, m)
        N = M.contract([2,3])
        self.assertEqual(N.groundset(), frozenset([1,4,5]))

        N = N.delete([5])
        self.assertEqual(N.groundset(), frozenset([1,4]))




class TestToric(unittest.TestCase):

    def test_without_Q(self):
        A = matrix(ZZ, [[-1, 1, 0, -1, 2, 7], [6, 1, -1, -2, 2, 5]])
        M = ToricArithmeticMatroid(A)

        self.assertEqual(M.full_rank(), M.rank(M.groundset()))

        self.assertEqual(M._Q, matrix(ZZ, 2, 0))
        self.assertEqual(M.rank([0]), 1)
        self.assertEqual(M.rank([0,1]), 2)
        self.assertEqual(M.rank([0,1,2]), 2)

        self.assertEqual(M._multiplicity([1,2]), 1)

        # minor
        M2 = M.minor(contractions=[], deletions=[1])
        self.assertEqual(M2.groundset(), frozenset([0,2,3,4,5]))
        self.assertEqual(M2._A, matrix(ZZ, [[-1, 0, -1, 2, 7], [6, -1, -2, 2, 5]]))

        M2 = M.minor(contractions=[1], deletions=[])
        self.assertEqual(M2.groundset(), frozenset([0,2,3,4,5]))
        self.assertEqual(M2._Q, matrix(ZZ, [[]]))
        self.assertEqual(M2._multiplicity([0,2]), 1)

        # dual
        M1 = M.dual()
        self.assertEqual(M1.rank(M1.groundset()), 4)
        M2 = M1.dual()
        self.assertTrue(M.equals(M2))

        # check representation
        self.assertTrue(M.check_representation(A))
        self.assertTrue(M.is_representable())
        self.assertEqual(M.representation(), A)
        self.assertEqual(M.representation(ordered_groundset=[5,4,3,2,1,0]), matrix(ZZ, [[7, 2, -1, 0, 1, -1], [5, 2, -2, -1, 1, 6]]))

        # orientability
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
        self.assertEqual(M.rank([1,2]), 1)

        self.assertEqual(M._multiplicity([0]), 1)
        self.assertEqual(M._multiplicity([1]), 3)
        self.assertEqual(M._multiplicity([1,2]), 1)


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
        self.assertEqual(M.rank([0,1]), 2)
        self.assertEqual(M.rank([0,1,2]), 2)

        self.assertEqual(M.multiplicity([1,2]), 1)

        # minor
        M2 = M.minor(contractions=[], deletions=[1])
        self.assertEqual(M2.groundset(), frozenset([0,2,3,4,5]))
        self.assertEqual(M2._A, matrix(ZZ, [[-1, 0, -1, 2, 7], [6, -1, -2, 2, 5]]))

        M2 = M.minor(contractions=[1], deletions=[])
        self.assertEqual(M2.groundset(), frozenset([0,2,3,4,5]))
        self.assertEqual(M2._Q, matrix(ZZ, [[]]))
        self.assertEqual(M2.multiplicity([0,2]), 1)


    def test_minor2(self):
        A = matrix(ZZ, [[2,1],[0,3]])
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


    def test_decomposition(self):
        M = ToricArithmeticMatroid(matrix(ZZ, [[1, 5, -1, 5], [1, 7, 2, 6], [0, -2, -2, -1]]))
        self.assertEqual(M.decomposition(), SetPartition([[2], [0,1,3]]))
        self.assertTrue(M.is_decomposable())


    def test_decomposition2(self):
        M = ToricArithmeticMatroid(matrix(ZZ, [[1, 5, -10, 50], [0, 0, 2, 6], [0, 0, -2, -1]]))
        self.assertEqual(M.decomposition(), SetPartition([[0,1], [2,3]]))
        self.assertTrue(M.is_decomposable())


    def test_decomposition3(self):
        M = ToricArithmeticMatroid(matrix(ZZ, [[1, 1], [0, 2]]))
        self.assertEqual(M.decomposition(), SetPartition([[0,1]]))
        self.assertFalse(M.is_decomposable())


    def test_poset_of_layers(self):
        M = ToricArithmeticMatroid(matrix(ZZ, [[1,0,1], [0,1,3]]))
        P = M.poset_of_layers()
        self.assertEqual(len(P), 7)
        self.assertTrue(P.has_bottom())
        self.assertFalse(P.has_top())
        self.assertEqual(len(P.cover_relations()), 10)
        self.assertTrue(P.is_isomorphic(Poset(data=(range(7), [(0,1), (0,2), (0,3), (1,4), (2,4), (2,5), (2,6), (3,4), (3,5), (3,6)]))))


    def test_poset_of_layers_relabeled(self):
        M = ToricArithmeticMatroid(matrix(ZZ, [[1,0,1], [0,1,3]]), ordered_groundset=['a', 'c', 'b'])
        P = M.poset_of_layers()
        self.assertEqual(len(P), 7)
        self.assertTrue(P.has_bottom())
        self.assertFalse(P.has_top())
        self.assertEqual(len(P.cover_relations()), 10)
        self.assertTrue(P.is_isomorphic(Poset(data=(range(7), [(0,1), (0,2), (0,3), (1,4), (2,4), (2,5), (2,6), (3,4), (3,5), (3,6)]))))
        self.assertEqual(set(S for (S, x) in P), set([(), ('a',), ('b',), ('c',), ('a', 'c', 'b'), ('a', 'b')]))


    def test_poset_of_layers_root_system_C2(self):
        # see [DGP17]
        A = matrix(ZZ, [[2,0,1,1,], [0,2,1,-1]])
        M = ToricArithmeticMatroid(A)
        P = M.poset_of_layers()
        self.assertEqual(len(P), 11)
        self.assertEqual(len(P.cover_relations()), 18)
        homology = P.subposet([a for a in P if a != P.bottom()]).order_complex(on_ints=True).homology()
        self.assertEqual([homology[i].ngens() for i in range(2)], [0,3])
        self.assertEqual(homology[1].ngens(), M.arithmetic_tutte_polynomial()(x=0, y=0))


    def test_poset_of_layers_root_system_C3(self):
        # see [DGP17]
        A = matrix(ZZ, [[2,0,0,1,1,1,1,0,0], [0,2,0,1,-1,0,0,1,1], [0,0,2,0,0,1,-1,1,-1]])
        M = ToricArithmeticMatroid(A)
        P = M.poset_of_layers()
        homology = P.subposet([a for a in P if a != P.bottom()]).order_complex(on_ints=True).homology()
        self.assertEqual([homology[i].ngens() for i in range(3)], [0,0,15])
        self.assertEqual(homology[2].ngens(), M.arithmetic_tutte_polynomial()(x=0, y=0))


    def test_poset_of_layers_root_system_B2(self):
        # see [DGP17]
        A = matrix(ZZ, [[1,0,1,1,], [0,1,1,-1]])
        M = ToricArithmeticMatroid(A)
        P = M.poset_of_layers()
        homology = P.subposet([a for a in P if a != P.bottom()]).order_complex(on_ints=True).homology()
        self.assertEqual([homology[i].ngens() for i in range(2)], [0,1])
        self.assertEqual(homology[1].ngens(), M.arithmetic_tutte_polynomial()(x=0, y=0))


    def test_poset_of_layers_root_system_B3(self):
        # see [DGP17]
        A = matrix(ZZ, [[1,0,0,1,1,1,1,0,0], [0,1,0,1,-1,0,0,1,1], [0,0,1,0,0,1,-1,1,-1]])
        M = ToricArithmeticMatroid(A)
        P = M.poset_of_layers()
        homology = P.subposet([a for a in P if a != P.bottom()]).order_complex(on_ints=True).homology()
        self.assertEqual([homology[i].ngens() for i in range(3)], [0,0,6])
        self.assertEqual(homology[2].ngens(), M.arithmetic_tutte_polynomial()(x=0, y=0))


    def test_poset_of_layers_root_system_D3(self):
        # see [DGP17]
        A = matrix(ZZ, [[1,1,1,1,0,0], [1,-1,0,0,1,1], [0,0,1,-1,1,-1]])
        M = ToricArithmeticMatroid(A)
        P = M.poset_of_layers()
        homology = P.subposet([a for a in P if a != P.bottom()]).order_complex(on_ints=True).homology()
        self.assertEqual([homology[i].ngens() for i in range(3)], [0,0,3])
        self.assertEqual(homology[2].ngens(), M.arithmetic_tutte_polynomial()(x=0, y=0))


    def test_different_posets_of_layers(self):
        # see [PP19, Section 8]
        A = matrix(ZZ, [[1,1,1,3], [0,5,0,5], [0,0,5,5]])
        B = matrix(ZZ, [[1,4,1,6], [0,5,0,5], [0,0,5,5]])
        C = matrix(ZZ, [[1,2,2,1], [0,5,0,5], [0,0,5,-5]])

        # check signed Hermite normal form
        self.assertEqual(shnf(A), matrix(ZZ, [[1,1,1,-3], [0,5,0,-5], [0,0,5,-5]]))
        self.assertEqual(shnf(B), matrix(ZZ, [[1,1,1,-1], [0,5,0,5], [0,0,5,-5]]))
        self.assertEqual(shnf(C), matrix(ZZ, [[1,2,2,1], [0,5,0,5], [0,0,5,-5]]))

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
        A = matrix(ZZ, [[1,0,1], [0,1,3]])
        M = ToricArithmeticMatroid(A)
        P = M.arithmetic_independence_poset()

        self.assertEqual(len(P), 9)
        self.assertEqual(len(P.cover_relations()), 13)
        self.assertEqual(P.rank(), 2)
        self.assertTrue(P.has_bottom())


    def test_independence_poset_relabeled(self):
        M = ToricArithmeticMatroid(matrix(ZZ, [[1,0,1], [0,1,3]]), ordered_groundset=['a', 'c', 'b'])
        P = M.arithmetic_independence_poset()

        self.assertEqual(len(P), 9)
        self.assertEqual(len(P.cover_relations()), 13)
        self.assertEqual(P.rank(), 2)
        self.assertTrue(P.has_bottom())
        self.assertEqual(set(S for (S, x) in P), set([(), ('a',), ('b',), ('c',), ('a', 'c'), ('a', 'b'), ('c', 'b')]))


    def test_nonshellable_poset_of_layers(self):
        # see [PP19, Section 8]
        A = matrix(ZZ, [[1,1,1,-3], [0,5,0,-5], [0,0,5,-5]])
        B = matrix(ZZ, [[1,1,1,-1], [0,5,0,5], [0,0,5,-5]])
        C = matrix(ZZ, [[1,2,2,1], [0,5,0,5], [0,0,5,-5]])

        for X in [A, B, C]:
            M = ToricArithmeticMatroid(X)
            P = M.poset_of_layers()
            homology = P.subposet([a for a in P if a != P.bottom()]).order_complex(on_ints=True).homology()

            from sage.homology.homology_group import HomologyGroup
            self.assertEqual(homology, {
                0: HomologyGroup(0, ZZ),    # trivial group
                1: HomologyGroup(1, ZZ, [5]),   # Z/5Z
                2: HomologyGroup(48, ZZ),   # Z^48
                })


    def test_nonshellable_arithmetic_independence_poset(self):
        # see [PP19, Section 8]
        A = matrix(ZZ, [[1,1,1,-3], [0,5,0,-5], [0,0,5,-5]])
        B = matrix(ZZ, [[1,1,1,-1], [0,5,0,5], [0,0,5,-5]])
        C = matrix(ZZ, [[1,2,2,1], [0,5,0,5], [0,0,5,-5]])

        for X in [A, B, C]:
            M = ToricArithmeticMatroid(X)
            P = M.arithmetic_independence_poset()
            homology = P.subposet([a for a in P if a != P.bottom()]).order_complex(on_ints=True).homology()

            from sage.homology.homology_group import HomologyGroup
            self.assertEqual(homology, {
                0: HomologyGroup(0, ZZ),    # trivial group
                1: HomologyGroup(1, ZZ, [5]),   # Z/5Z
                2: HomologyGroup(73, ZZ),   # Z^73
                })



class TestReduction(unittest.TestCase):

    def test_example(self):
        E = [1,2,3,4,5,6]

        def rk(X):
            return min(2, len(X))

        def m(X):
            if len(X) == 0:
                return 1
            elif len(X) == 1:
                if X == frozenset([1]) or X == frozenset([2]):
                    return 2
                else:
                    return 1
            elif sum(1 for i in X if i not in [1,2]) >= 2:
                return 1
            elif len(X) == 2:
                if X == frozenset([1,2]):
                    return 4
                else:
                    return 2
            elif 3 in X:
                # X = {1,2,3}
                assert X == frozenset([1,2,3])
                return 1
            else:
                # X = {1,2,j} with j > 3
                return 2

        M = ArithmeticMatroid(E, rk, m)
        self.assertTrue(M.is_valid())

        M1 = M.reduction()
        self.assertEqual(M1.multiplicity([1,2,3]), 2)

        for X in powerset(E):
            if tuple(sorted(X)) != (1,2,3):
                self.assertEqual(M1.multiplicity(X), M.multiplicity(X))

        self.assertFalse(M1.is_valid())

        self.assertTrue(M1.is_torsion_free())
        self.assertTrue(M1.is_surjective())


    def test_dual(self):
        A = matrix(ZZ, [[-1, 1, 0, -1, 2, 7], [6, 1, -1, -2, 2, 5]])
        M = representation_to_matroid(A)
        M2 = ArithmeticMatroid(M.groundset(), M.rank, lambda X: M._multiplicity(X)**2)

        self.assertTrue(M.is_valid())
        self.assertTrue(M2.is_valid())
        self.assertFalse(M2.is_representable())

        self.assertFalse(M.dual().reduction().is_isomorphic(M.reduction()))
        self.assertTrue(M.dual().reduction().is_isomorphic(M.reduction().dual()))
        self.assertTrue(M2.dual().reduction().is_isomorphic(M2.reduction().dual()))



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
                self.assertEqual(shnf(A), shnf(U*A*S))

        else:
            S = diagonal_matrix([sage.misc.prandom.choice([1, -1]) for i in range(n)])
            self.assertEqual(shnf(A), shnf(U*A*S))


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
                    self.assertEqual(shnf(U*A*S), A)

        for u in [1,2,4,5]:
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
            self.assertEqual(shnf(U*A*S), shnf(A))
            self.assertEqual(shnf(U*A), shnf(A))


    def test_not_reducible(self):
        # in this example, the matrix seems a block matrix but it isn't
        A = matrix(ZZ, [[2,2,0,0,0],[0,4,0,0,4],[0,0,4,4,0],[0,0,0,8,8]])
        for diag in itertools.product([1, -1], repeat=5):
            S = diagonal_matrix(diag)
            self.assertEqual(shnf(A*S), shnf(A))



if __name__ == '__main__':
    unittest.main()
