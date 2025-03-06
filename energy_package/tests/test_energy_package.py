"""
Unit and regression test for the energy_package package.
"""

# Import package, test suite, and other packages as needed
import sys

import pytest

import energy_package


def test_energy_package_imported():
    """Sample test, will always pass so long as import statement worked."""
    assert "energy_package" in sys.modules
