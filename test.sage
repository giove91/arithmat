import unittest
from arithmetic_matroid import *


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
        self.assertTrue(M.is_realizable())
        
    
    def test_valid(self):
        # matroid with realization
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
        self.assertTrue(M.is_realizable())
        # print M.realization()


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
        # valid, not realizable, not orientable
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
        self.assertFalse(M.is_realizable())
        self.assertFalse(M.is_orientable())
        self.assertEqual(M.realization(), None)
    
    
    def test_realization_to_matroid(self):
        # realizable with a 2x2 matrix
        for a in xrange(6):
            for b in xrange(1, 6):
                A = matrix(ZZ, [[1, a], [0, b]])
                M = realization_to_matroid(A)
                
                self.assertTrue(M.is_valid())
                self.assertTrue(M.is_realizable())
                self.assertTrue(M.is_orientable())
                
                R = ZZ['x, y']
                x, y = R._first_ngens(2)
                self.assertEqual(M.arithmetic_tutte_polynomial(), x**2 + x*(gcd(a,b)-1) + b-gcd(a,b))
    
    
    def test_realization_to_matroid2(self):
        A = matrix(ZZ, [[-1, -29, -1, 1], [1, -1, 0, 1]])
        M = realization_to_matroid(A)
        
        self.assertTrue(M.is_valid())
        self.assertTrue(M.is_realizable())
        self.assertTrue(M.is_orientable())


    def test_realization_to_matroid3(self):
        A = matrix(ZZ, [[-1,  1,  0,  0, -1], [ 6,  1,  1, -1, -1]])
        M = realization_to_matroid(A)
        
        self.assertTrue(M.is_valid())
        self.assertTrue(M.is_realizable())
        self.assertTrue(M.is_orientable())
    
    
    def test_realization_to_matroid4(self):
        A = matrix(ZZ, [[ 2,  2,  1,  0,  0], [ 1,  5, -1,  1, -2], [-2,  1,  0, -1, -1]])
        M = realization_to_matroid(A)
        
        self.assertTrue(M.is_valid())
        self.assertTrue(M.is_realizable())
        self.assertTrue(M.is_orientable())


    def test_realization_to_matroid_random(self):
        A = random_matrix(ZZ,4,6)
        M = realization_to_matroid(A)
        
        self.assertTrue(M.is_valid())
        self.assertTrue(M.is_realizable())
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
        self.assertFalse(M.is_realizable())
        self.assertTrue(M.is_orientable())
    
    
    def test_non_realizable(self):
        A = matrix(ZZ, [[-1,  1,  0, -1], [ 6,  1, -1, -2]])
        M = realization_to_matroid(A)
        M2 = ArithmeticMatroid(M.groundset(), M._rank, lambda X: M._multiplicity(X)**2)
        
        self.assertTrue(M2.is_valid())
        self.assertTrue(M.is_realizable())
        self.assertFalse(M2.is_realizable())
    
    def test_non_realizable2(self):
        A = matrix(ZZ, [[-1,  1,  0, -1, 2, 7], [ 6,  1, -1, -2, 2, 5]])
        M = realization_to_matroid(A)
        M2 = ArithmeticMatroid(M.groundset(), M._rank, lambda X: M._multiplicity(X)**2)
        
        self.assertTrue(M2.is_valid())
        self.assertTrue(M.is_realizable())
        self.assertFalse(M2.is_realizable())
    
    
    def test_not_orientable(self):
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
        self.assertFalse(M.is_realizable())
        self.assertFalse(M.is_orientable())
    
    
    def test_non_essential(self):
        A = matrix(ZZ, [[-1,  1,  0], [ 6,  1, -1], [2, -3, 0], [1, 2, 3], [-1, 0, 0]])
        M = realization_to_matroid(A)
        
        self.assertTrue(M.is_valid())
        self.assertTrue(M.is_realizable())



class TestDual(unittest.TestCase):
    
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
        
        M3 = M1._minor()
        self.assertIsInstance(M3, DualArithmeticMatroid)
        # self.assertEqual(M1, M3) # this is not an equality for (non-arithmetic) matroids
        self.assertTrue(M3.is_valid())
    
    
    def test_linear_matroid(self):
        A = matrix(QQ, [[-1,  1,  0, -1, 2, 7], [ 6,  1, -1, -2, 2, 5]])
        
        def m(X):
            return 1
        
        M = LinearArithmeticMatroid(A, multiplicity_function=m)
        self.assertTrue(M.is_valid())
        
        M1 = M.dual()
        self.assertTrue(M1.is_valid())
        
        M2 = M1.dual()
        self.assertTrue(M2.is_valid())
        self.assertEqual(M, M2)
        self.assertNotEqual(M, M1)
        self.assertNotEqual(M1, M2)



if __name__ == '__main__':
    unittest.main()

