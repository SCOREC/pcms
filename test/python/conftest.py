"""
Shared pytest fixtures for pcms Python API tests.
"""

import pytest

import pcms


@pytest.fixture(scope="session")
def omega_h_lib():
    """Create Omega_h library instance.

    Session-scoped so the (potentially expensive) Omega_h initialization
    happens exactly once per test run.
    """
    if not pcms.features["omega_h"]:
        pytest.skip("pcms was built without the Omega_h backend")
    omega_h = pytest.importorskip("PyOmega_h")
    return omega_h.OmegaHLibrary()


@pytest.fixture(scope="session")
def world(omega_h_lib):
    """Create Omega_h world communicator.

    Session-scoped; the world is obtained from the session-scoped library
    and shared by all tests that need mesh-building capabilities.
    """
    return omega_h_lib.world()
