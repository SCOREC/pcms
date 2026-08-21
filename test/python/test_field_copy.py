"""Test copying Omega_h field data."""

import numpy as np
import pytest
import pcms
import PyOmega_h as omega_h


class TestFieldCopy:
    """Tests for field data copying between Omega_h-backed fields."""

    @pytest.mark.parametrize("dim, order, num_components", [
        (2, 1, 1),
        (2, 2, 1),
    ])
    def test_copy(self, world, dim, order, num_components):
        """Test copying omega_h field data."""
        nx = 100
        ny = 100 if dim > 1 else 0
        nz = 100 if dim > 2 else 0

        # Build mesh
        mesh = omega_h.build_box(
            world,
            omega_h.Family.SIMPLEX,
            1.0, 1.0, 1.0,
            nx, ny, nz,
            False
        )

        # Create factory and layout
        factory = pcms.LagrangeFunctionSpace.from_mesh(
            mesh, order, num_components, pcms.CoordinateSystem.Cartesian
        )

        # Create original field and set data
        original = factory.create_field()
        ndata = original.get_num_dof_holders() * original.get_num_components()

        # Create sequential array of IDs
        ids = np.arange(ndata, dtype=np.float64)
        original.set_dof_holder_data(ids)

        # Create copied field and copy data
        copied = factory.create_field()
        copier = pcms.Copy(factory, factory)
        copier.apply(original, copied)

        # Get copied data
        copied_array = copied.get_dof_holder_data()

        # Verify the copy
        assert len(copied_array) == ndata, \
            f"Expected {ndata} elements, got {len(copied_array)}"

        # Check that all values match
        differences = np.abs(ids - copied_array)
        num_matches = np.sum(differences < 1e-12)
        assert num_matches == ndata, \
            f"Only {num_matches}/{ndata} elements matched"
