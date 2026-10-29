"""
pytest configuration for the Arithmat tests.

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

import hashlib
import os
import secrets

import pytest

# Initialize Sage once before any test module loads.
import sage.all  # noqa: F401
from sage.misc.randstate import set_random_seed

from helpers import MATROID_NAMES, make_matroid

# Random seed unless ARITHMAT_SEED is set for reproducibility.
# Store it in the environment so pytest-xdist workers use the same seed.
SEED = int(os.environ.setdefault("ARITHMAT_SEED", str(secrets.randbelow(2**32))))


def pytest_report_header(config):
    return f"Reproduce with ARITHMAT_SEED={SEED}"


@pytest.fixture(autouse=True)
def seed_sage_random(request):
    # Reseed deterministically before every test
    digest = hashlib.sha256(f"{SEED}:{request.node.nodeid}".encode()).digest()
    set_random_seed(int.from_bytes(digest[:8], "big"))


@pytest.fixture(params=MATROID_NAMES)
def matroid(request):
    return make_matroid(request.param)
