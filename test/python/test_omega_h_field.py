"""
Field-centric tests for Omega_h-backed function spaces.
"""

import numpy as np
import pytest
import pcms
import PyOmega_h as omega_h

class TestOmegaHField:
    """Tests for Omega_h-backed Field operations."""

    @staticmethod
    def _build_mesh(world, dim, nx=10):
        """Build a box mesh of the requested dimension."""
        ny = nx if dim > 1 else 0
        nz = nx if dim > 2 else 0
        return omega_h.build_box(
            world, omega_h.Family.SIMPLEX, 1.0, 1.0, 1.0, nx, ny, nz, False,
        )

    @pytest.mark.parametrize("dim, order, num_components", [
        (2, 1, 1), (2, 2, 1),
    ])
    def test_field_methods(self, world, dim, order, num_components):
        """Create an Omega_h-backed Field and exercise the public Field API."""
        mesh = self._build_mesh(world, dim)
        factory = pcms.LagrangeFunctionSpace.from_mesh(
            mesh, order, num_components, pcms.CoordinateSystem.Cartesian
        )
        field = factory.create_field()

        assert field.get_num_components() == num_components
        assert field.get_num_dof_holders() > 0

        coords = field.get_dof_holder_coordinates()
        assert coords.shape[0] == field.get_num_dof_holders()
        assert coords.shape[1] == dim

        ndata = field.get_num_dof_holders() * num_components
        test_data = np.arange(ndata, dtype=np.float64)
        expected_data = test_data.copy()
        field.set_dof_holder_data(test_data)
        test_data.fill(-999.0)
        np.testing.assert_allclose(field.get_dof_holder_data(), expected_data)

    @pytest.mark.parametrize("dim, order, num_components", [
        (2, 1, 1), (2, 2, 1),
    ])
    def test_field_transfer(self, world, dim, order, num_components):
        """Transfer data between identical Omega_h-backed function spaces."""
        if num_components != 1:
            pytest.skip("Multi-component transfer not yet supported")

        mesh = self._build_mesh(world, dim)
        source_space = pcms.LagrangeFunctionSpace.from_mesh(
            mesh, order, num_components, pcms.CoordinateSystem.Cartesian
        )
        target_space = pcms.LagrangeFunctionSpace.from_mesh(
            mesh, order, num_components, pcms.CoordinateSystem.Cartesian
        )
        source = source_space.create_field()
        target = target_space.create_field()

        coords = source.get_dof_holder_coordinates()
        source_data = np.zeros(source.get_num_dof_holders(), dtype=np.float64)
        for i in range(source.get_num_dof_holders()):
            source_data[i] = np.sum(coords[i, :dim])
        source.set_dof_holder_data(source_data)

        interp = pcms.Interpolator(source_space, target_space)
        interp.apply(source, target)
        np.testing.assert_allclose(
            target.get_dof_holder_data(), source_data, atol=1e-14
        )

    @pytest.mark.parametrize("dim, order, num_components", [
        (2, 1, 1), (2, 2, 1),
    ])
    def test_field_evaluation(self, world, dim, order, num_components):
        """Evaluate an Omega_h-backed field at explicit query points."""
        mesh = self._build_mesh(world, dim)
        factory = pcms.LagrangeFunctionSpace.from_mesh(
            mesh, order, num_components, pcms.CoordinateSystem.Cartesian
        )
        field = factory.create_field()

        mesh_coords = mesh.coords()
        num_verts = mesh.nverts()
        num_owned = field.get_num_dof_holders()
        ndata = num_owned * num_components

        def test_func(x, y, z, component):
            return np.sin(5.0 * x * y) + float(component)

        test_data = np.zeros(ndata, dtype=np.float64)
        for i in range(min(num_verts, num_owned)):
            x = mesh_coords[i * dim + 0] if dim >= 1 else 0.0
            y = mesh_coords[i * dim + 1] if dim >= 2 else 0.0
            z = mesh_coords[i * dim + 2] if dim >= 3 else 0.0
            for c in range(num_components):
                idx = i * num_components + c
                test_data[idx] = test_func(x, y, z, c)

        if order == 2 and num_owned > num_verts:
            for i in range(num_verts, num_owned):
                for c in range(num_components):
                    idx = i * num_components + c
                    test_data[idx] = float(c) + 0.5

        field.set_dof_holder_data(test_data)

        if dim == 2:
            eval_coords = np.array([
                [0.5, 0.5],
                [0.25, 0.25],
                [0.75, 0.75],
                [0.1, 0.9],
                [0.9, 0.1],
            ], dtype=np.float64)
        else:
            eval_coords = np.array([
                [0.5, 0.5, 0.5],
                [0.25, 0.25, 0.25],
                [0.75, 0.75, 0.75],
                [0.1, 0.9, 0.1],
                [0.9, 0.1, 0.9],
            ], dtype=np.float64)

        request = pcms.EvaluationRequest.from_coordinates(eval_coords)
        evaluator = factory.create_point_evaluator(request)
        eval_values = np.zeros((eval_coords.shape[0], num_components),
                               dtype=np.float64)
        evaluator.evaluate(field, eval_values)

        for i in range(eval_coords.shape[0]):
            for c in range(num_components):
                val = eval_values[i, c]
                assert not np.isnan(val)
                assert not np.isinf(val)
