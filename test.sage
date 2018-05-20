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
        self.assertEqual(M.r, 9)
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
        # valid, not realizable
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
        self.assertEqual(M.realization(), None)
    
    
    def test_realization_to_matroid(self):
        # realizable with a 2x2 matrix
        for a in xrange(6):
            for b in xrange(1, 6):
                A = matrix(ZZ, [[1, a], [0, b]])
                M = realization_to_matroid(A)
                
                self.assertTrue(M.is_valid())
                self.assertTrue(M.is_realizable())
    
    def test_realization_to_matroid2(self):
        A = matrix(ZZ, [[-1, -29, -1, 1], [1, -1, 0, 1]])
        M = realization_to_matroid(A)
        
        self.assertTrue(M.is_valid())
        self.assertTrue(M.is_realizable())


    def test_realization_to_matroid3(self):
        A = matrix(ZZ, [[-1,  1,  0,  0, -1], [ 6,  1,  1, -1, -1]])
        M = realization_to_matroid(A)
        
        self.assertTrue(M.is_valid())
        self.assertTrue(M.is_realizable())
    
    
    def test_realization_to_matroid4(self):
        A = matrix(ZZ, [[ 2,  2,  1,  0,  0], [ 1,  5, -1,  1, -2], [-2,  1,  0, -1, -1]])
        M = realization_to_matroid(A)
        
        self.assertTrue(M.is_valid())
        self.assertTrue(M.is_realizable())


    def test_realization_to_matroid_random(self):
        A = random_matrix(ZZ,4,6)
        M = realization_to_matroid(A)
        
        self.assertTrue(M.is_valid())
        self.assertTrue(M.is_realizable())


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

if __name__ == '__main__':
    unittest.main()
