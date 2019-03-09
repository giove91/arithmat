# -*- coding: utf-8 -*-
"""
Arithmetic matroids classes.

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

import itertools
import operator
import copy

def normal_form(A):
    """
    S-normal form of an integer matrix A.
    This is a normal form up to left multiplication by an invertible matrix and change of sign of the columns.
    A matrix in S-normal form is also in left Hermite normal form.
    """
    r = A.nrows()
    n = A.ncols()
    G_basis = []    # basis of G
    m = 0   # rank of A[:,:j]

    for j in xrange(n):
        A.echelonize()
        q = A[m,j] if m < r else 0  # pivot
        print "Column:", j
        print "Pivot:", q

        phi = []
        for S in G_basis:
            H, U = (A[:,:j] * S[:j,:j]).echelon_form(transformation=True)
            assert H == A[:,:j]
            phi.append(U)

        G_basis.append(matrix(ZZ, [[0 if k != i else -1 if k == j else 1 for k in xrange(n)] for i in xrange(n)]))
        phi.append(matrix(ZZ, [[0 if k != i else -1 if k < m else 1 for k in xrange(r)] for i in xrange(r)]))

        for i in reversed(range(m)):
            # find possible values of A[i,j]
            x = Integer(A[i,j])
            if q > 0:
                x %= q

            orbit = [x]
            columns = [A[:,j]]
            construction = [identity_matrix(n)] # which element of G gives a certain element

            new_elements = True
            while new_elements:
                new_elements = False
                for h, U in enumerate(phi):
                    for k, v in enumerate(columns):
                        w = U*v
                        y = Integer(w[i,0])
                        if q > 0:
                            y %= q
                        if y not in orbit:
                            orbit.append(y)
                            columns.append(w)
                            construction.append(G_basis[h] * construction[k])
                            new_elements = True

            # select the minimal possible value
            u = min(orbit)
            k = orbit.index(u)
            S = construction[k]

            # change A
            A = (A*S).echelon_form()
            assert A[i,j] == u

            # update the stabilizer G
            # TODO

        if q != 0:
            m += 1

    return A


A = matrix(ZZ, [[5, 8], [0, 3]])
print normal_form(A)
