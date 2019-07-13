# -*- coding: utf-8 -*-
"""
Signed Hermite normal form.

Copyright (C) 2019 Giovanni Paolini
Copyright (C) 2019 Roberto Pagaria

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

import networkx as nx
import itertools


def poset_of_layers(A):
    """
    Poset of layers of the central toric arrangement defined by the integer matrix A.
    Uses Lenz's algorithm.
    """
    A = A.transpose()
    E = range(A.nrows())

    data = {}

    # compute Smith normal forms of all submatrices
    for S in powerset(E):
        D, U, V = A[S,:].smith_form()   # D == U*A[S,:]*V
        diagonal = [D[i,i] if i < D.ncols() else 0 for i in xrange(len(S))]
        data[tuple(S)] = (diagonal, U)

    # generate al possible elements of the poset of layers
    elements = {tuple(S): list(vector(ZZ, x) for x in itertools.product(*(range(max(data[tuple(S)][0][i], 1)) for i in xrange(len(S))))) for S in powerset(E)}

    for l in elements.itervalues():
        for v in l:
            v.set_immutable()

    possible_layers = list((S, x) for (S, l) in elements.iteritems() for x in l)
    uf = DisjointSet(possible_layers)

    cover_relations = []

    for (S, l) in elements.iteritems():
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

    # find actual layers and cover relations
    layers = [a for a in possible_layers if uf.find(a) == a]
    cover_relations = set((uf.find(a), uf.find(b)) for (a,b) in cover_relations)

    return Poset(data=(layers, cover_relations), cover_relations=True)
