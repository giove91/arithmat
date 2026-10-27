"""
Helpers shared by several test modules.

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

import operator
from functools import reduce

from arithmat.arithmetic_matroid import ArithmeticMatroid


def representation_to_matroid(A):
    """
    Given an integer matrix A, return the associated ArithmeticMatroid.
    """
    E = range(A.ncols())

    def rk(X):
        return A[:, list(X)].rank()

    def m(X):
        return reduce(operator.mul, [d for d in A[:, list(X)].elementary_divisors() if d != 0], 1)

    return ArithmeticMatroid(E, rk, m)
