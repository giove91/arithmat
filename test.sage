import unittest
from arithmetic_matroid import ArithmeticMatroid


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
    
    
    def test_valid(self):
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


if __name__ == '__main__':
    unittest.main()
