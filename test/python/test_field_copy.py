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
            mesh, order, num_components)

        # Create original field and set data
        original = factory.create_function()
        layout = original.get_layout()
        ndata = layout.get_num_owned_dof_holder() * layout.get_num_components()

        # Create sequential array of IDs
        ids = np.arange(ndata, dtype=np.float64).reshape(
            -1, original.get_layout().get_num_components())
        original.set_dof_holder_data(ids)

        # Create copied field and copy data
        copied = factory.create_function()
        copier = pcms.Copy(factory, factory)
        copier.apply(original, copied)

        # Get copied data
        copied_array = copied.get_dof_holder_data().get_values()

        # Verify the copy
        assert copied_array.shape == ids.shape, \
            f"Expected shape {ids.shape}, got {copied_array.shape}"

        # Check that all values match
        differences = np.abs(ids - copied_array)
        num_matches = np.sum(differences < 1e-12)
        assert num_matches == ndata, \
            f"Only {num_matches}/{ndata} elements matched"
