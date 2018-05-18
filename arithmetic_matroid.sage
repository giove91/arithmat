import sage.matroids.matroid
import itertools
import networkx as nx
import operator
import sys
from fractions import gcd

"""
TODO

Comparison:

    def __hash__(self)
    def __eq__(self, other)
    def __ne__(self, other)

In Cythonized classes, use __richcmp__() instead of __eq__(), __ne__().

Copying, loading, saving:

    def __copy__(self)
    def __deepcopy__(self, memo={})
    def __reduce__(self)
"""


class ArithmeticMatroid(sage.matroids.matroid.Matroid):
    def __init__(self, E, rk, m):
        self.E = E
        self.rk = rk
        self.m = m
        
        self.r = self._rank(self.E)

    def __repr__(self):
        return "Arithmetic matroid of rank %d on %d elements" % (self.r, len(self.E))

    def groundset(self):
        return self.E

    def _rank(self, X):
        return self.rk(X)
    
    def _multiplicity(self, X):
        return self.m(X)
    
    def is_independent_from(self, v, X):
        return self._rank(X+[v]) != self._rank(X)
    
    def is_dependent_from(self, v, X):
        return not self.is_independent_from(v, X)
    
    
    def is_valid(self):
        if not super(ArithmeticMatroid, self).is_valid():
            return False
        
        # check axioms for arithmetic matroids
        for X in powerset(self.E):
            for v in self.E:
                if v not in X:
                    if self.is_dependent_from(v, X):
                        # check axiom 1
                        if self._multiplicity(X) % self._multiplicity(X+[v]) != 0:
                            return False
                    
                    else:
                        # check axiom 2
                        if self._multiplicity(X+[v]) % self._multiplicity(X) != 0:
                            return False
        
        for Y in powerset(self.E):
            for X in powerset(Y):
                T = []
                F = []
                for y in Y:
                    if y not in X:
                        if self.is_dependent_from(y, X):
                            T.append(y)
                        else:
                            F.append(y)
                
                if len(F) + self._rank(X) == self._rank(Y):
                    # (X, Y) is a molecule
                    
                    # check axiom 3
                    if self._multiplicity(X) * self._multiplicity(Y) != self._multiplicity(X+F) * self._multiplicity(X+T):
                        return False
                    
                    # check axiom P
                    if (-1)**(len(T)) * sum((-1)**(len(Y)-len(Z)-len(X)) * self._multiplicity(X+Z) for Z in powerset([z for z in Y if z not in X])) < 0:
                        return False
        
        return True
    
    
    def check_realization(self, A):
        """
        Check if the given matrix is a realization for the matroid.
        """
        E = self.E
        r = self.r
        n = len(E)
        
        if A.ncols() != n:
            return False
        
        for S in powerset(range(n)):
            T = [E[i] for i in S]   # corresponding subset of E
            if A[:,S].rank() != self._rank(T):
                # print >> sys.stderr, "Not realizable, rank of %r is incorrect" % T
                return False

            if reduce(operator.mul, [d for d in A[:,S].elementary_divisors() if d != 0], 1) != self._multiplicity(T):
                # print >> sys.stderr, "Not realizable, multiplicity of %r is incorrect" % T
                return False
        
        return True
        
    
    def realization_surjective(self):
        """
        Find a realization (if it exists) for a surjective matroid (m(E)=1).
        """
        E = self.E
        r = self.r
        n = len(E)
        assert self._multiplicity(E) == 1
        
        B = list(sorted(self.basis()))
        # print "Basis:", B
                
        # find bipartite graph
        edges = [(x,y) for x in B for y in E if y not in B and self._rank([z for z in B if z != x]+[y]) == r]
        
        spanning_forest = nx.Graph()
        
        # find spanning forest
        uf = DisjointSet(self.E)
        for (x,y) in edges:
            if uf.find(x) != uf.find(y):
                spanning_forest.add_edge(x,y)
                uf.union(x,y)
        
        # print "Graph:", edges
        # print "Spanning forest:", spanning_forest.edges()
        
        # compute entries of matrix A
        def entry(i,j):
            x = B[i]
            y = E[j]
            
            if y in B:
                return self._multiplicity(B) if x == y else 0
            
            elif (x,y) in edges:
                return self._multiplicity([z for z in B if z != x]+[y])
            
            else:
                return 0
        
        A = matrix(ZZ, r, n, entry)
        # print A
        
        B_to_index = {B[i]: i for i in xrange(r)}
        E_to_index = {E[j]: j for j in xrange(n)}
        
        
        graph = spanning_forest
        while graph.number_of_edges() < len(edges):
            # find all paths in the graph
            paths = nx.all_pairs_dijkstra_path(graph)
            for (x,y) in sorted(edges, key = lambda (x,y): len(paths[x][y])):
                if len(paths[x][y]) == 2:
                    # (x,y) is in the graph
                    assert (x,y) in graph.edges()
                    continue
            
                i = B_to_index[x]
                j = E_to_index[y]
                
                rows = [B_to_index[z] for z in paths[x][y][::2]]
                columns = [E_to_index[z] for z in paths[x][y][1::2]]
                
                # print x, y
                # print "rows:", rows
                # print "columns:", columns
                
                new_tuple = [z for z in B + paths[x][y] if z not in B or z not in paths[x][y]]
                # print "new_tuple:", new_tuple
                expected_mult = self._multiplicity(new_tuple) * self._multiplicity(B)**(len(rows)-1) if self._rank(new_tuple) == r else 0
                if abs(A[rows,columns].determinant()) != expected_mult:
                    # change sign
                    # print "change sign!"
                    # print A[rows,columns].determinant()
                    A[i,j] = -A[i,j]
                    
                    if abs(A[rows,columns].determinant()) != expected_mult:
                        # print A
                        # print A[rows,columns].determinant(), expected_mult
                        return None
                
                graph.add_edge(x,y)
                break
        
        D, P, Q = A.smith_form()
        res = Q.inverse()[:r,:]
        res = matrix(ZZ, res)
        
        # print >> sys.stderr, "Candidate realization:"
        # print >> sys.stderr, res
        
        # check if this is indeed a realization
        if not self.check_realization(res):
            return None
        
        return res


    def all_realizations(self):
        """
        Generator of all non-equivalent essential realizations.
        """
        E = self.E
        r = self.r
        n = len(E)
        if self._multiplicity(E) == 1:
            res = self.realization_surjective()
            if res is not None:
                yield res
            return
        
        # construct "reduced" matroid
        denominator = reduce(gcd, [self._multiplicity(B) for B in self.bases()], 0)
        
        def m_bar(X):
            return reduce(gcd, [self._multiplicity(B) for B in self.bases() if self._rank(X) == self._rank([x for x in X if x in B])], 0) // denominator
        
        M = ArithmeticMatroid(E, self._rank, m_bar)
        
        if not M.is_valid():
            return
        
        # get realization of "reduced" matroid
        A = M.realization_surjective()
        
        if A is None:
            return
        
        # try all left Hermite normal forms
        for H in hermite_normal_forms(r, self._multiplicity(E)):
            if self.check_realization(H*A):
                yield H*A

    
    def realization(self):
        """
        Compute any essential realization.
        Returns None if the matroid is not realizable.
        """
        for A in self.all_realizations():
            return A
        return None

    
    def is_realizable(self):
        return self.realization() is not None



def hermite_normal_forms(r, det):
    """
    Generate all r x r integer matrices in (left) Hermite normal form
    with the given determinant.
    """
    if r == 0:
        if det == 1:
            yield matrix(ZZ, 0, [])
    
    else:
        for d in divisors(det):
            for A in hermite_normal_forms(r-1, det//d):
                for column in itertools.product(range(d), repeat=r-1):
                    yield matrix(ZZ, r, r, lambda i, j:
                        A[i,j] if i < r-1 and j < r-1
                        else d if i == r-1 and j == r-1
                        else column[i] if j == r-1
                        else 0
                    )
    
    


def realization_to_matroid(A):
    """
    Given an integer matrix A, return the associated ArithmeticMatroid.
    """
    E = range(A.ncols())
    
    def rk(X):
        return A[:,list(X)].rank()
    
    def m(X):
        return reduce(operator.mul, [d for d in A[:,list(X)].elementary_divisors() if d != 0], 1)
    
    return ArithmeticMatroid(E, rk, m)



if __name__ == '__main__':
    E = [1,2,3,4,5]

    def rk(X):
        return min(2, len(X))


    def m(X):
        if len(X) == 2 and all(x in [3,4,5] for x in X):
            return 2
        else:
            return 1


    M = ArithmeticMatroid(E, rk, m)

    print M
    print "Valid:", M.is_valid()

    print M.realization()
    
    
