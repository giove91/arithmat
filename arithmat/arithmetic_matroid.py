# -*- coding: utf-8 -*-
"""
Arithmetic matroids classes.

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

from __future__ import absolute_import

import itertools
import operator
import copy
import networkx as nx
from functools import reduce

from sage.structure.sage_object import SageObject
from sage.matroids.matroid import Matroid
from sage.matroids.advanced import *
from sage.combinat.set_partition import SetPartition
from sage.sets.disjoint_set import DisjointSet
from sage.combinat.posets.posets import Poset

from sage.rings.integer_ring import ZZ
from sage.matrix.constructor import matrix, vector
from sage.matrix.special import identity_matrix, block_matrix

from sage.misc.misc import powerset
from sage.arith.misc import divisors, gcd

from .shnf import signed_hermite_normal_form


class ArithmeticMatroidMixin(SageObject):
    """
    A general mixin for arithmetic matroids, that can be added to any Matroid subclass of Sage.
    """

    def __init__(self, *args, **kwargs):
        # get multiplicity function
        try:
            multiplicity = kwargs.pop('multiplicity_function')
        except KeyError:
            multiplicity = None # multiplicity function must be set later

        super(ArithmeticMatroidMixin, self).__init__(*args, **kwargs)
        self._multiplicity = multiplicity

    def _repr_(self):
        return "Arithmetic matroid of rank %d on %d elements" % (self.full_rank(), len(self.groundset()))


    def __hash__(self):
        return hash((self.groundset(), self.full_rank(), self._multiplicity(frozenset()), self._multiplicity(self.groundset())))

    def __eq__(self, other):
        if not isinstance(other, ArithmeticMatroidMixin):
            return False

        return super(ArithmeticMatroidMixin, self).__eq__(other) and (self._multiplicity == other._multiplicity)

    def __ne__(self, other):
        return not self == other


    def __copy__(self):
        N = super(ArithmeticMatroidMixin, self).__copy__()
        N._multiplicity = self._multiplicity
        N.__class__ = type(self)
        return N

    def __deepcopy__(self, *args, **kwargs):
        N = super(ArithmeticMatroidMixin, self).__deepcopy__(*args, **kwargs)
        N._multiplicity = copy.deepcopy(self._multiplicity)
        N.__class__ = type(self)
        return N

    def __reduce__(self):
        raise TypeError("unfortunately, functions cannot be saved reliably, so this class doesn't have load/save support.")



    def multiplicity(self, X=None):
        """
        Return the multiplicity of X.
        """
        if X is None:
            return self.full_multiplicity()
        if not isinstance(X, frozenset):
            X = frozenset(X)
        return self._multiplicity(X)


    def full_multiplicity(self):
        """
        Return the multiplicity of the groundset.
        """
        return self._multiplicity(self.groundset())


    def _is_isomorphism(self, other, morphism):
        """
        Version of is_isomorphism() that does no type checking.
        (see Matroid.is_isomorphism)
        """
        for X in powerset(self.groundset()):
            Y = frozenset(morphism[e] for e in X)
            if self.rank(X) != other.rank(Y) or self.multiplicity(X) != other.multiplicity(Y):
                return False
        return True


    def _is_isomorphic(self, other, certificate=False):
        """
        Test if self is isomorphic to other.
        Internal version that performs no checks on input (see Matroid.is_isomorphic).
        """
        if certificate:
            # TODO
            raise NotImplementedError

        if len(self.groundset()) != len(other.groundset()) or self.full_rank() != other.full_rank() or self.full_multiplicity() != other.full_multiplicity():
            return False

        # TODO: try to make more efficient
        E = list(self.groundset())
        for perm in itertools.permutations(other.groundset()):
            morphism = {e: perm[i] for i, e in enumerate(E)}
            if self.is_isomorphism(other, morphism):
                return True

        return False


    def _is_independent_from(self, v, X):
        return self.rank(X.union([v])) != self.rank(X)

    def _is_dependent_from(self, v, X):
        return not self._is_independent_from(v, X)


    def is_valid(self):
        """
        Check if the arithmetic matroid axioms are satisfied.
        """
        if not super(ArithmeticMatroidMixin, self).is_valid():
            # check validity of the underlying matroid
            return False

        E = self.groundset()

        # check axioms for arithmetic matroids
        for X in powerset(E):
            X = frozenset(X)
            for v in E:
                if v not in X:
                    if self._is_dependent_from(v, X):
                        # check axiom 1
                        if self.multiplicity(X) % self.multiplicity(X.union([v])) != 0:
                            return False

                    else:
                        # check axiom 2
                        if self.multiplicity(X.union([v])) % self.multiplicity(X) != 0:
                            return False

        for Y in powerset(E):
            Y = frozenset(Y)

            for X in powerset(Y):
                X = frozenset(X)

                T = []
                F = []

                for y in Y:
                    if y not in X:
                        if self._is_dependent_from(y, X):
                            T.append(y)
                        else:
                            F.append(y)

                if len(F) + self.rank(X) == self.rank(Y):
                    # (X, Y) is a molecule

                    T = frozenset(T)
                    F = frozenset(F)

                    # check axiom 3
                    if self.multiplicity(X) * self.multiplicity(Y) != self.multiplicity(X.union(F)) * self.multiplicity(X.union(T)):
                        return False

                    # check axiom P
                    if (-1)**(len(T)) * sum((-1)**(len(Y)-len(Z)-len(X)) * self.multiplicity(X.union(Z)) for Z in powerset(Y.difference(X))) < 0:
                        return False

        return True


    def is_torsion_free(self):
        """
        Check if the matroid is torsion-free, i.e. m({}) = 1.
        """
        return self.multiplicity(frozenset()) == 1


    def is_surjective(self):
        """
        Check if the matroid is surjective, i.e. m(E) = 1.
        """
        return self.full_multiplicity() == 1


    def is_gcd(self):
        """
        Check if the matroid satisfies the gcd property, defined in [DM13, Section 3].
        """
        return all(self.multiplicity(X) == reduce(gcd, (self.multiplicity(I) for I in powerset(X) if self.rank(X) == self.rank(I) and self.is_independent(I)), 0) for X in powerset(self.groundset()))


    def is_strong_gcd(self):
        """
        Check if the matroid satisfies the gcd property, defined in [PP19, Section 3].
        """
        return all(self.multiplicity(X) == reduce(gcd, (self.multiplicity(B) for B in self.bases() if len(B.intersection(X)) == self.rank(X)), 0) for X in powerset(self.groundset()))


    def minor(self, contractions=None, deletions=None):
        # get minor as a (non-arithmetic) matroid
        matroid = super(ArithmeticMatroidMixin, self).minor(contractions, deletions)

        contractions = list(contractions) if contractions else []
        deletions = list(deletions) if deletions else []

        if isinstance(matroid, MinorMatroid):
            # return an instance of MinorArithmeticMatroid
            return MinorArithmeticMatroid(self, contractions, deletions)

        else:
            # we use the same (arithmetic) class here
            matroid.__class__ = type(self)

            # add multiplicity function
            matroid._multiplicity = lambda X : self._multiplicity(frozenset(contractions).union(X))

            return matroid


    def dual(self):
        # get dual as a (non-arithmetic) matroid
        matroid = super(ArithmeticMatroidMixin, self).dual()

        if isinstance(matroid, DualMatroid):
            # return an instance of DualArithmeticMatroid
            return DualArithmeticMatroid(self)

        else:
            # we use the same (arithmetic) class here
            matroid.__class__ = type(self)

            # add multiplicity function
            matroid._multiplicity = lambda X : self._multiplicity(self.groundset().difference(X))

            return matroid


    def reduction(self):
        """
        Return reduction of the matroid, as defined in [PP19, Section 4].
        """
        d = reduce(gcd, [self.multiplicity(B) for B in self.bases()], 0)

        def m_bar(X):
            return reduce(gcd, [self.multiplicity(B) for B in self.bases() if self.rank(X) == self.rank(X.intersection(B))], 0) // d

        return ArithmeticMatroid(self.groundset(), self._rank, multiplicity_function=m_bar)


    def check_representation(self, A, ordered_groundset=None, check_bases=False):
        """
        Check if the given matrix is a representation for the matroid.
        If check_bases==True, check that the multiplicity is correct only on the bases.
        """
        r = self.full_rank()
        n = len(self.groundset())

        if ordered_groundset is not None:
            # use the groundset in the given order
            E = ordered_groundset
            assert frozenset(E) == self.groundset()
            assert len(E) == len(self.groundset())
        else:
            # to sort the groundset
            E = list(sorted(self.groundset()))

        if A.ncols() != n:
            return False

        for S in powerset(range(n)):
            T = frozenset(E[i] for i in S)   # corresponding subset of E
            if A[:,S].rank() != self.rank(T):
                return False

            if check_bases and len(T) != r and self.rank(T) < r:
                # skip multiplicity check
                continue

            if reduce(operator.mul, [d for d in A[:,S].elementary_divisors() if d != 0], 1) != self.multiplicity(T):
                return False

        return True


    def _representation_surjective(self, ordered_groundset=None, check_bases=False):
        """
        Find a representation (if it exists) for a surjective matroid (m(emptyset)=m(E)=1).
        If check_bases==True, find a representation of a matroid (E,rk,m')
        such that m'(B)=m(B) for every basis B.
        Return None if no representation exists.
        """
        assert self.full_multiplicity() == 1

        r = self.full_rank()
        n = len(self.groundset())

        if ordered_groundset is not None:
            # use the groundset in the given order
            E = ordered_groundset
            assert frozenset(E) == self.groundset()
            assert len(E) == len(self.groundset())
        else:
            # sort the groundset
            E = list(sorted(self.groundset()))

        B = self.basis()

        # find bipartite graph
        edges = [(x,y) for x in B for y in E if y not in B and self.is_basis(B.difference([x]).union([y]))]

        spanning_forest = nx.Graph()

        # find spanning forest
        uf = DisjointSet(E)
        for (x,y) in edges:
            if uf.find(x) != uf.find(y):
                spanning_forest.add_edge(x,y)
                uf.union(x,y)

        # fix an order of B
        B_ordered = list(sorted(B))

        # compute entries of matrix A
        def entry(i,j):
            x = B_ordered[i]
            y = E[j]

            if y in B:
                return self.multiplicity(B) if x == y else 0

            elif (x,y) in edges:
                return self.multiplicity(B.difference([x]).union([y]))

            else:
                return 0

        A = matrix(ZZ, r, n, entry)

        B_to_index = {B_ordered[i]: i for i in range(r)}
        E_to_index = {E[j]: j for j in range(n)}


        graph = spanning_forest
        while graph.number_of_edges() < len(edges):
            # find all paths in the graph
            paths = dict(nx.all_pairs_dijkstra_path(graph))
            for (x,y) in sorted(edges, key=lambda edge: len(paths[edge[0]][edge[1]])):
                if len(paths[x][y]) == 2:
                    # (x,y) is in the graph
                    assert (x,y) in graph.edges()
                    continue

                i = B_to_index[x]
                j = E_to_index[y]

                rows = [B_to_index[z] for z in paths[x][y][::2]]
                columns = [E_to_index[z] for z in paths[x][y][1::2]]

                new_tuple = [z for z in B_ordered + paths[x][y] if z not in B or z not in paths[x][y]]
                expected_mult = self.multiplicity(new_tuple) * self.multiplicity(B)**(len(rows)-1) if self.rank(new_tuple) == r else 0
                if abs(A[rows,columns].determinant()) != expected_mult:
                    # change sign
                    A[i,j] = -A[i,j]

                    if abs(A[rows,columns].determinant()) != expected_mult:
                        return None

                graph.add_edge(x,y)
                break

        D, U, V = A.smith_form()
        res = V.inverse()[:r,:]
        res = matrix(ZZ, res)

        # check if this is indeed a representation
        if not self.check_representation(res, ordered_groundset=ordered_groundset, check_bases=check_bases):
            return None

        return res


    def all_representations(self, ordered_groundset=None):
        """
        Generator of all non-equivalent essential representations.
        Uses the algorithm of [PP19, Section 5].
        """
        if self.multiplicity([]) > 1:
            raise NotImplementedError
            # TODO implement m({}) > 1

        r = self.full_rank()

        if self.full_multiplicity() == 1:
            res = self._representation_surjective(ordered_groundset=ordered_groundset)
            if res is not None:
                yield res
            return

        M = self.reduction()

        if not M.is_valid():
            return

        # get representation of reduced matroid
        A = M._representation_surjective(ordered_groundset=ordered_groundset)

        if A is None:
            return

        found_representations = set()

        # try all left Hermite normal forms
        for H in _hermite_normal_forms(r, self.multiplicity(self.groundset())):
            if self.check_representation(H*A):
                B = signed_hermite_normal_form(H*A)
                if B not in found_representations:
                    found_representations.add(B)
                    yield B


    def num_representations(self):
        """
        Compute the number of non-equivalent essential representations.
        """
        return sum(1 for _ in self.all_representations())


    def representation(self, ordered_groundset=None):
        """
        Compute any essential representation.
        Return None if the matroid is not representable.
        """
        for A in self.all_representations(ordered_groundset=ordered_groundset):
            return A
        return None


    def is_representable(self):
        """
        Determine if the matroid is a representable arithmetic matroid.
        """
        return self.representation() is not None


    def is_orientable(self):
        """
        Determine if the matroid is an orientable arithmetic matroid, as defined in [Pag18].
        """
        M = self.reduction() # note that this might not be an arithmetic matroid

        return M._representation_surjective(check_bases=True) is not None
        # TODO maybe it is not necessary to check (on the bases) that the result is a representation


    def arithmetic_tutte_polynomial(self, x=None, y=None):
        """
        Return the arithmetic Tutte polynomial of the matroid.
        """
        r = self.full_rank()

        a = x
        b = y
        R = ZZ['x, y']
        x, y = R._first_ngens(2)
        T = R(0)
        for X in powerset(self.groundset()):
            T += self.multiplicity(X) * (x-1) ** (r - self.rank(X)) * (y-1) ** (len(X) - self.rank(X))
        if a is not None and b is not None:
            T = T(a, b)
        return T



class ArithmeticMatroid(ArithmeticMatroidMixin, RankMatroid):
    def __init__(self, groundset, rank_function, multiplicity_function):
        # take multiplicity function as third positional argument
        return super(ArithmeticMatroid, self).__init__(groundset, rank_function, multiplicity_function=multiplicity_function)


class LinearArithmeticMatroid(ArithmeticMatroidMixin, LinearMatroid):
    def _repr_(self):
        return "Linear arithmetic matroid of rank %d on %d elements" % (self.full_rank(), len(self.groundset()))


class BasisArithmeticMatroid(ArithmeticMatroidMixin, BasisMatroid):
    def __init__(self, M=None, *args, **kwargs):
        if isinstance(M, ArithmeticMatroidMixin):
            # extract multiplicity function from the given arithmetic matroid
            kwargs['multiplicity_function'] = M._multiplicity
        super(BasisArithmeticMatroid, self).__init__(M=M, *args, **kwargs)

    def _repr_(self):
        return "Basis arithmetic matroid of rank %d on %d elements" % (self.full_rank(), len(self.groundset()))



class MinorArithmeticMatroid(ArithmeticMatroidMixin, MinorMatroid):
    """
    Minor of an arithmetic matroid.
    """
    def __init__(self, *args, **kwargs):
        super(ArithmeticMatroidMixin, self).__init__(*args, **kwargs)

    def _repr_(self):
        s = "M"
        if self._contractions:
            s += r" / " + setprint_s(self._contractions)
        if self._deletions:
            s += r" \ " + setprint_s(self._deletions)
        s += ", where M is " + repr(self._matroid)
        return s

    def __eq__(self, other):
        return (self._contractions == other._contractions) and (self._deletions == other._deletions) and (self._matroid == other._matroid)

    def __reduce__(self):
        return super(ArithmeticMatroidMixin, self).__reduce__()

    def _multiplicity(self, X):
        return self._matroid._multiplicity(self._contractions.union(X))

    def minor(self, contractions=None, deletions=None):
        if contractions is None:
            contractions = []
        if deletions is None:
            deletions = []
        return MinorArithmeticMatroid(self._matroid, self._contractions.union(contractions), self._deletions.union(deletions))



class DualArithmeticMatroid(ArithmeticMatroidMixin, DualMatroid):
    """
    Dual of an arithmetic matroid.
    """
    def __init__(self, matroid):
        if not isinstance(matroid, ArithmeticMatroidMixin):
            raise TypeError("no arithmetic matroid provided to take dual of.")
        self._matroid = matroid

    def _repr_(self):
        return "Dual of '" + repr(self._matroid) + "'"

    def __eq__(self, other):
        return self._matroid == other._matroid

    def __reduce__(self):
        return super(ArithmeticMatroidMixin, self).__reduce__()


    def _multiplicity(self, X):
        return self._matroid._multiplicity(self.groundset().difference(X))

    def dual(self):
        return self._matroid

    def minor(self, contractions=None, deletions=None):
        # Assumption: if self._matroid cannot make a dual, neither can its minor.
        return DualArithmeticMatroid(self._matroid.minor(contractions=deletions, deletions=contractions))



class ToricArithmeticMatroid(ArithmeticMatroidMixin, Matroid):
    """
    Arithmetic matroid defined by a given representation.
    The representation is defined up to equivalence.
    """

    def __init__(self, arrangement_matrix, torus_matrix=None, ordered_groundset=None):
        self._A = arrangement_matrix
        self._Q = torus_matrix if torus_matrix is not None else matrix(ZZ, arrangement_matrix.nrows(), 0)

        assert self._A.nrows() == self._Q.nrows()
        self._normalize()

        if ordered_groundset is None:
            ordered_groundset = range(arrangement_matrix.ncols())
        else:
            assert len(set(ordered_groundset)) == len(ordered_groundset)
            assert len(ordered_groundset) == arrangement_matrix.ncols()

        self._groundset = frozenset(ordered_groundset) # non-ordered groundset
        self._E = ordered_groundset                    # ordered groundset
        self._groundset_to_index = {e: i for (i,e) in enumerate(ordered_groundset)} # dictionary groundset -> columns

    def _normalize(self):
        """
        Put Q in Smith normal form (this also changes A).
        """
        D, U, V = self._Q.smith_form()  # D = U*Q*V
        self._A = U * self._A
        self._Q = D[:, :D.rank()] # remove zero columns from Q

        while self._Q.ncols() > 0 and self._Q[0,0] == 1:
            # delete first row from Q and A
            self._Q = self._Q[1:,1:]
            self._A = self._A[1:,:]


    def groundset(self):
        return self._groundset

    def ordered_groundset(self):
        return self._E

    def arrangement_matrix(self):
        return self._A

    def torus_matrix(self):
        return self._Q


    def _rank(self, X):
        T = block_matrix(ZZ, [[self._A[:, [self._groundset_to_index[e] for e in X]], self._Q]])
        return T.rank() - self._Q.ncols()

    def _multiplicity(self, X):
        T = block_matrix(ZZ, [[self._A[:, [self._groundset_to_index[e] for e in X]], self._Q]])
        return reduce(operator.mul, [d for d in T.elementary_divisors() if d != 0], 1)


    def _repr_(self):
        return "Toric arithmetic matroid of rank %d on %d elements" % (self.full_rank(), len(self.groundset()))


    def __hash__(self):
        return hash((self._E, self._A, self._Q))

    def __eq__(self, other):
        if not isinstance(other, ToricArithmeticMatroid):
            return False

        return self._E == other._E and self._A == other._A and self._Q == other._Q


    def __copy__(self):
        return ToricArithmeticMatroid(arrangement_matrix=self._A, torus_matrix=self._Q, ordered_groundset=self._E)

    def __deepcopy__(self, *args, **kwargs):
        return ToricArithmeticMatroid(arrangement_matrix=copy.deepcopy(self._A), torus_matrix=copy.deepcopy(self._Q), ordered_groundset=copy.deepcopy(self._E))

    def __reduce__(self):
        # TODO
        raise NotImplementedError

    def is_valid(self):
        return True

    def minor(self, contractions=None, deletions=None):
        contractions = list(contractions) if contractions else []
        deletions = list(deletions) if deletions else []

        new_groundset = [e for e in self._E if e not in contractions+deletions]
        A2 = copy.copy(self._A[:, [self._groundset_to_index[e] for e in new_groundset]])
        Q2 = block_matrix(ZZ, [[self._A[:, [self._groundset_to_index[e] for e in contractions]], self._Q]])
        return ToricArithmeticMatroid(arrangement_matrix=A2, torus_matrix=Q2, ordered_groundset=new_groundset)


    def dual(self):
        T = block_matrix(ZZ, [[self._A[:, [self._groundset_to_index[e] for e in self._E]], self._Q]]).transpose()
        I = identity_matrix(ZZ, T.nrows())

        # create temporary names for new groundset elements
        temp_elements = [-i for i in range(self._Q.ncols())]
        while len(frozenset(temp_elements).intersection(self.groundset())) > 0:
            temp_elements = [e-1 for e in temp_elements]

        M = ToricArithmeticMatroid(arrangement_matrix=I, torus_matrix=T, ordered_groundset=list(self._E)+temp_elements)
        return M.minor(deletions=temp_elements)


    def representation(self, ordered_groundset=None):
        """
        Compute any essential representation.
        Return None if the matroid is not representable.
        """
        # first try to return self._A
        if self._Q.ncols() == 0:
            if ordered_groundset is not None:
                # use the groundset in the given order
                E = ordered_groundset
                assert frozenset(E) == self.groundset()
                assert len(E) == len(self.groundset())
            else:
                # use self._E (which is ordered)
                E = self._E

            # return self._A, with shuffled columns
            return self._A[:, [self._groundset_to_index[e] for e in E]]
            # TODO should also return Q, when m({}) > 1?

        return super(ToricArithmeticMatroid, self).representation(ordered_groundset=ordered_groundset)


    def is_orientable(self):
        """
        Determine if the matroid is an orientable arithmetic matroid, as defined in [Pag18].
        """
        if self._Q.ncols() == 0:
            return True

        else:
            return super(ToricArithmeticMatroid, self).is_orientable()


    def is_equivalent(self, other, morphism=None):
        """
        Check if the two ToricArithmeticMatroids are equivalent,
        i.e. the defining representations are equivalent (see [PP19, Section 2]).
        If morphism is None, assume that the groundsets coincide.
        """
        if not isinstance(other, ToricArithmeticMatroid):
            raise TypeError("can only test for equivalence between toric arithmetic matroids.")

        if self._Q.ncols() > 0 or other._Q.ncols() > 0:
            # TODO
            raise NotImplementedError

        if morphism is None:
            assert self.groundset() == other.groundset()
            morphism = {e: e for e in self.groundset()}

        E = self._E

        # take matrices in Hermite normal form, removing zero rows
        # (copy is needed to make matrices mutable)
        M = copy.copy(self._A.echelon_form(include_zero_rows=False))
        N = copy.copy(other._A[:, [other._groundset_to_index[morphism[e]] for e in self._E]].echelon_form(include_zero_rows=False))

        # choose a basis
        B = self.basis()
        if not other.is_basis(frozenset(morphism[e] for e in B)):
            return False

        # find bipartite graph
        edges = []
        for x in B:
            for y in E:
                C = B.difference([x]).union([y])
                if y not in B and self.is_basis(C):
                    if not other.is_basis(frozenset(morphism[e] for e in C)):
                        return False

                    edges.append((x,y))

        spanning_forest = nx.Graph()

        # find spanning forest
        uf = DisjointSet(E)
        for (x,y) in edges:
            if uf.find(x) != uf.find(y):
                spanning_forest.add_edge(x,y)
                uf.union(x,y)

        B_indices = list(sorted(self._groundset_to_index[e] for e in B))

        M1 = M[:, B_indices].inverse() * M
        N1 = N[:, B_indices].inverse() * N

        def change_signs(A, A1):
            for (i,j) in nx.edge_dfs(spanning_forest):
                (x,y) = (i,j) if i in B else (j,i)
                if A1[x,y] < 0:
                    if j in B:
                        # change sign of row j and column j
                        A1[j,:] *= -1
                        A1[:,j] *= -1

                        A[j,:] *= -1
                        A[:,j] *= -1

                        assert A1 == A[:, B_indices].inverse() * A

                    else:
                        # change sign of column j
                        A1[:,j] *= -1

                        A[:,j] *= -1

        change_signs(M, M1)
        change_signs(N, N1)

        return M.echelon_form() == N.echelon_form()


    def decomposition(self):
        """
        Find the decomposition of the matroid as a direct sum of indecomposable matroids.
        Return a partition of the groundset.
        Uses the algorithm of [PP19, Section 7].
        """
        B = self.basis()
        new_groundset = list(B) + list(self.groundset().difference(B))

        # construct matrix with permuted columns
        columns = [vector(self._A[:,self._groundset_to_index[e]]) for e in new_groundset]
        A = matrix(ZZ, columns).transpose().echelon_form()

        uf = DisjointSet(self.groundset())

        for i in range(A.nrows()):
            for j in range(i+1, A.ncols()):
                if A[i,j] != 0:
                    uf.union(new_groundset[i], new_groundset[j])

        return SetPartition(uf)


    def is_decomposable(self):
        """
        Check if the matroid is decomposable.
        """
        return len(self.decomposition()) > 1


    def is_indecomposable(self):
        """
        Check if the matroid is indecomposable.
        """
        return not self.is_decomposable()


    def poset_of_layers(self):
        """
        Compute the poset of layers of the associated toric arrangement, using Lenz's algorithm [Len17a].
        """
        # TODO: implement for Q != 0
        if self._Q.ncols() > 0:
            raise NotImplementedError


        A = self._A.transpose()
        E = range(A.nrows())

        data = {}

        # compute Smith normal forms of all submatrices
        for S in powerset(E):
            D, U, V = A[S,:].smith_form()   # D == U*A[S,:]*V
            diagonal = [D[i,i] if i < D.ncols() else 0 for i in range(len(S))]
            data[tuple(S)] = (diagonal, U)

        # generate al possible elements of the poset of layers
        elements = {tuple(S): list(vector(ZZ, x) for x in itertools.product(*(range(max(data[tuple(S)][0][i], 1)) for i in range(len(S))))) for S in powerset(E)}

        for l in elements.values():
            for v in l:
                v.set_immutable()

        possible_layers = list((S, x) for (S, l) in elements.items() for x in l)
        uf = DisjointSet(possible_layers)

        cover_relations = []

        for (S, l) in elements.items():
            diagonal_S, U_S = data[S]
            rk_S = A[S,:].rank()

            for s in S:
                i = S.index(s)  # index where the element s appears in S
                T = tuple(t for t in S if t != s)

                diagonal_T, U_T = data[T]
                rk_T = A[T,:].rank()

                for x in l:
                    h = (S, x)

                    y = U_S**(-1) * x
                    z = U_T * vector(ZZ, y[:i].list() + y[i+1:].list())
                    w = vector(ZZ, (a % diagonal_T[j] if diagonal_T[j] > 0 else 0 for j, a in enumerate(z)))
                    w.set_immutable()

                    ph = (T, w)

                    if rk_S == rk_T:
                        uf.union(h, ph)

                    else:
                        cover_relations.append((ph, h))

        # find representatives for layers (we keep the representative (S,x) with maximal S)
        root_to_representative_dict = {}
        for root, subset in uf.root_to_elements_dict().items():
            S, x = max(subset, key=lambda Sx: len(Sx[0]))
            S_labeled = tuple(self._E[i] for i in S)
            root_to_representative_dict[root] = (S_labeled, x)

        # get layers and cover relations
        layers = root_to_representative_dict.values()
        cover_relations = set(
            (root_to_representative_dict[uf.find(a)], root_to_representative_dict[uf.find(b)])
            for (a,b) in cover_relations)

        return Poset(data=(layers, cover_relations), cover_relations=True)


    def arithmetic_independence_poset(self):
        """
        Compute the poset of (arithmetic) independent sets of the associated central toric arrangement.
        This is defined in [Len17b, Definition 5], [Mar18, Definitions 2.1 and 2.2], [DD18, Section 7].
        Notice that it is not the same as the independence poset of the underlying matroid.
        """
        # TODO: implement for Q != 0
        if self._Q.ncols() > 0:
            raise NotImplementedError

        A = self._A.transpose()
        data = {}

        # compute Smith normal forms of all submatrices
        for S_labeled in self.independent_sets():
            S = tuple(sorted(self._groundset_to_index[e] for e in S_labeled))
            D, U, V = A[S,:].smith_form()   # D == U*A[S,:]*V
            diagonal = [D[i,i] if i < D.ncols() else 0 for i in range(len(S))]
            data[tuple(S)] = (diagonal, U)

        # generate all elements of the poset
        elements = {S: list(vector(ZZ, x) for x in itertools.product(*(range(max(data[tuple(S)][0][i], 1)) for i in range(len(S))))) for S in data.keys()}

        for l in elements.values():
            for v in l:
                v.set_immutable()

        all_elements = list((tuple(self._E[i] for i in S), x) for (S, l) in elements.items() for x in l)
        cover_relations = []

        for (S, l) in elements.items():
            diagonal_S, U_S = data[S]
            S_labeled = tuple(self._E[i] for i in S)

            for s in S:
                i = S.index(s)  # index where the element s appears in S
                T = tuple(t for t in S if t != s)
                T_labeled = tuple(self._E[i] for i in T)

                diagonal_T, U_T = data[T]

                for x in l:
                    y = U_S**(-1) * x
                    z = U_T * vector(ZZ, y[:i].list() + y[i+1:].list())
                    w = vector(ZZ, (a % diagonal_T[j] if diagonal_T[j] > 0 else 0 for j, a in enumerate(z)))
                    w.set_immutable()

                    cover_relations.append(((T_labeled, w), (S_labeled, x)))

        return Poset(data=(all_elements, cover_relations), cover_relations=True)




def _hermite_normal_forms(r, det):
    """
    Generate all r x r integer matrices in (left) Hermite normal form
    with the given determinant.
    """
    if r == 0:
        if det == 1:
            yield matrix(ZZ, 0, [])

    else:
        for d in divisors(det):
            for A in _hermite_normal_forms(r-1, det//d):
                for column in itertools.product(range(d), repeat=r-1):
                    yield matrix(ZZ, r, r, lambda i, j:
                        A[i,j] if i < r-1 and j < r-1
                        else d if i == r-1 and j == r-1
                        else column[i] if j == r-1
                        else 0
                    )


# copied from sage.matroids.utilities
def setprint_s(X, toplevel=False):
    """
    Create the string for use by ``setprint()``.

    INPUT:

    - ``X`` -- any Python object
    - ``toplevel`` -- (default: ``False``) indicates whether this is a
      recursion or not.

    OUTPUT:

    A string representation of the object, with nice notation for sets and
    frozensets.

    EXAMPLES::

        sage: from sage.matroids.utilities import setprint_s
        sage: L = [{1, 2, 3}, {1, 2, 4}, {2, 3, 4}, {4, 1, 3}]
        sage: setprint_s(L)
        '[{1, 2, 3}, {1, 2, 4}, {1, 3, 4}, {2, 3, 4}]'

    The ``toplevel`` argument only affects strings, to mimic ``print``'s
    behavior::

        sage: X = 'abcd'
        sage: setprint_s(X)
        "'abcd'"
        sage: setprint_s(X, toplevel=True)
        'abcd'
    """
    if isinstance(X, frozenset) or isinstance(X, set):
        return '{' + ', '.join(sorted(setprint_s(x) for x in X)) + '}'
    elif isinstance(X, dict):
        return '{' + ', '.join(sorted(setprint_s(key) + ': ' + setprint_s(val)
                                      for key, val in X.items())) + '}'
    elif isinstance(X, str):
        if toplevel:
            return X
        else:
            return "'" + X + "'"
    elif hasattr(X, '__iter__') and not isinstance(X, SageObject):
        return '[' + ', '.join(sorted(setprint_s(x) for x in X)) + ']'
    else:
        return repr(X)
