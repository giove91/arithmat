import sage.matroids.matroid
import itertools
import networkx as nx
import operator
import sys

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
    
    
    def realization(self):
        # FIXME only works for m(E)=1
        
        E = self.E
        r = self.r
        n = len(E)
        
        B = list(sorted(self.basis()))
        # print "Basis:", B
        
        # C = matrix(r, n-r, lambda i,j: self._rank(B[:i]+B[i+1:]+[J[j]]) == r)
        # print C
        
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
            
            expected_mult = self._multiplicity([z for z in B + paths[x][y] if z not in B or z not in paths[x][y]])
            if abs(A[rows,columns].determinant()) != expected_mult:
                # change sign
                # print "change sign!"
                A[i,j] = -A[i,j]
                
                if abs(A[rows,columns].determinant()) != expected_mult:
                    return None
            
            graph.add_edge(x,y)
        
        D, P, Q = A.smith_form()
        res = Q.inverse()[:r,:]
        res = matrix(ZZ, res)
        
        print >> sys.stderr, "Candidate realization:"
        print >> sys.stderr, res
        
        # check if this is indeed a realization
        for S in powerset(range(n)):
            T = [E[i] for i in S]   # corresponding subset of E
            if res[:,S].rank() != self._rank(T):
                print >> sys.stderr, "Not realizable, rank of %r is incorrect" % T
                return None

            if reduce(operator.mul, [d for d in res[:,S].elementary_divisors() if d != 0], 1) != self._multiplicity(T):
                print >> sys.stderr, "Not realizable, multiplicity of %r is incorrect" % T
                return None
        
        return res

    
    def is_realizable(self):
        return self.realization() is not None


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
    
    print realization_to_matroid(matrix(ZZ, 2, [1,2,3, 4,5,6]))

    
