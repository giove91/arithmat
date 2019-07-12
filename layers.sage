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
    r = A.nrows()
    n = A.ncols()
    E = range(n)

    data = {tuple(S): A[:,S].smith_form() for S in powerset(E)}
    # D, U, V = A[:,S].smith_form()
    #   assert D == U*A[:,S]*V
    
    elements = {tuple(S): list(vector(ZZ, x) for x in itertools.product(*(range(max(data[tuple(S)][0][i,i], 1) if i < data[tuple(S)][0].nrows() else 1) for i in xrange(len(S))))) for S in powerset(E)}

    for l in elements.itervalues():
        for v in l:
            v.set_immutable()

    #print list(elements.iteritems())
    uf = DisjointSet((S, x) for (S, l) in elements.iteritems() for x in l)
    print uf

    for (S, l) in elements.iteritems():
        D, U, V = data[S]
        for s in S:
            T = tuple(t for t in S if t != s)
            for x in l:
                h = (S, x)
                ph = (T, )  # TODO




A = matrix(ZZ, [[1,0,1], [0,1,3]])
print poset_of_layers(A)