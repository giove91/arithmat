# -*- coding: utf-8 -*-
"""
Signed Hermite normal form.

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

from sage.matrix.special import diagonal_matrix, identity_matrix

def signed_hermite_normal_form(A):
    """
    Signed Hermite normal form of an integer matrix A, see [PP19, Section 6].
    This is a normal form up to left-multiplication by invertible matrices and change of sign of the columns.
    A matrix in signed Hermite normal form is also in left Hermite normal form.
    """
    r = A.nrows()
    n = A.ncols()
    G_basis = []    # Z_2-basis of G
    m = 0   # rank of A[:,:j]

    for j in range(n):
        A = A.echelon_form()
        q = A[m,j] if m < r else 0  # pivot

        phi = []
        for S in G_basis:
            H, U = (A[:,:j] * S[:j,:j]).echelon_form(transformation=True)
            assert H == A[:,:j]
            phi.append(U)

        G_basis.append(diagonal_matrix([-1 if i == j else 1 for i in range(n)]))
        phi.append(diagonal_matrix([-1 if i < m else 1 for i in range(r)]))
        assert len(G_basis) == len(phi)

        for i in reversed(range(m)):
            # find possible values of A[i,j]
            x = A[i,j]
            if q > 0:
                x %= q

            orbit = [x]
            columns = [A[:,j]]
            construction = [identity_matrix(n)] # which element of G gives a certain element of the orbit

            new_elements = True
            while new_elements:
                new_elements = False
                for h, U in enumerate(phi):
                    for k, v in enumerate(columns):
                        w = U*v
                        y = w[i,0]
                        if q > 0:
                            y %= q
                        if y not in orbit:
                            orbit.append(y)
                            columns.append(w)
                            construction.append(G_basis[h] * construction[k])
                            new_elements = True

            assert len(orbit) in [1,2,4]

            # find action of G on the orbit
            action = []
            for h, U in enumerate(phi):
                if q > 0:
                    action.append({x: (U*columns[k])[i,0] % q for k, x in enumerate(orbit)})
                else:
                    action.append({x: (U*columns[k])[i,0] for k, x in enumerate(orbit)})

            # select the minimal possible value
            u = min(orbit)
            k = orbit.index(u)
            S = construction[k]

            # change A
            A = (A*S).echelon_form()
            assert A[i,j] == u

            # update the stabilizer G
            G_new_basis = []    # basis for the new stabilizer
            new_phi = []    # value of phi on the new basis elements
            complement_basis = {} # dictionary of the form {x: h}, where G_basis[h] sends u to x
            for h, S in enumerate(G_basis):
                if action[h][u] == u:
                    # this basis element fixes u
                    G_new_basis.append(S)
                    new_phi.append(phi[h])

                elif action[h][u] in complement_basis:
                    # we already encountered a basis element which sends u to action[h][u]
                    G_new_basis.append(S * G_basis[complement_basis[action[h][u]]])
                    new_phi.append(phi[h] * phi[complement_basis[action[h][u]]])

                elif len(complement_basis) == 2:
                    # the product of the two basis elements of the complement sends u to action[h][u]
                    x, y = list(complement_basis)
                    G_new_basis.append(S * G_basis[complement_basis[x]] * G_basis[complement_basis[y]])
                    new_phi.append(phi[h] * phi[complement_basis[x]] * phi[complement_basis[y]])

                else:
                    # add S to the basis of the complement
                    complement_basis[action[h][u]] = h

            assert len(G_new_basis) + len(complement_basis) == len(G_basis)
            G_basis = G_new_basis
            phi = new_phi
            assert len(G_basis) == len(phi)

        if q != 0:
            m += 1

    return A
