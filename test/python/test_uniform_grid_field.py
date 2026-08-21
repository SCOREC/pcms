"""
Field-centric tests for uniform-grid-backed function spaces.
"""

import numpy as np
import pytest
import pcms
import PyOmega_h as omega_h

class TestUniformGridField:
    """Tests for UniformGrid-backed Field operations (no Omega_h needed)."""

    def test_field_creation(self):
        """Create a 2D field from a uniform-grid function space."""
        grid = pcms.UniformGrid2D()
        grid.bot_left = [0.0, 0.0]
        grid.edge_length = [10.0, 10.0]
        grid.divisions = [4, 4]

        factory = pcms.LagrangeFunctionSpace.from_uniform_grid(
            grid, 1, pcms.CoordinateSystem.Cartesian
        )
        field = factory.create_field()
        expected = (grid.divisions[0] + 1) * (grid.divisions[1] + 1)
        assert field.get_num_components() == 1
        assert field.get_num_dof_holders() == expected

    def test_data_operations(self):
        """Set and get flat DOF data through Field."""
        grid = pcms.UniformGrid2D()
        grid.bot_left = [0.0, 0.0]
        grid.edge_length = [10.0, 10.0]
        grid.divisions = [2, 2]

        factory = pcms.LagrangeFunctionSpace.from_uniform_grid(
            grid, 1, pcms.CoordinateSystem.Cartesian
        )
        field = factory.create_field()
        data = np.arange(field.get_num_dof_holders(), dtype=np.float64)
        expected_data = data.copy()
        field.set_dof_holder_data(data)
        data.fill(-999.0)
        np.testing.assert_allclose(field.get_dof_holder_data(), expected_data)

    def test_coordinates_2d(self):
        """Expose 2D DOF-holder coordinates through Field."""
        grid = pcms.UniformGrid2D()
        grid.bot_left = [0.0, 0.0]
        grid.edge_length = [10.0, 10.0]
        grid.divisions = [2, 3]

        factory = pcms.LagrangeFunctionSpace.from_uniform_grid(
            grid, 1, pcms.CoordinateSystem.Cartesian
        )
        field = factory.create_field()
        coords = field.get_dof_holder_coordinates()
        expected = (grid.divisions[0] + 1) * (grid.divisions[1] + 1)

        assert isinstance(coords, np.ndarray)
        assert coords.shape == (expected, 2)
        np.testing.assert_allclose(coords[0], [0.0, 0.0])
        np.testing.assert_allclose(coords[-1], [10.0, 10.0])

    def test_closest_cell(self):
        """UniformGrid helpers remain available for grid setup."""
        grid = pcms.UniformGrid2D()
        grid.bot_left = [0.0, 0.0]
        grid.edge_length = [10.0, 10.0]
        grid.divisions = [4, 4]

        cell_id = grid.closest_cell_id(np.array([5.0, 5.0]))
        cell_id_outside = grid.closest_cell_id(np.array([-1.0, -1.0]))
        assert isinstance(cell_id, (int, np.integer))
        assert isinstance(cell_id_outside, (int, np.integer))

    def test_coordinates_3d(self):
        """Expose 3D DOF-holder coordinates and data through Field."""
        grid = pcms.UniformGrid3D()
        grid.bot_left = [0.0, 0.0, 0.0]
        grid.edge_length = [10.0, 10.0, 10.0]
        grid.divisions = [2, 2, 2]

        factory = pcms.LagrangeFunctionSpace.from_uniform_grid(
            grid, 1, pcms.CoordinateSystem.Cartesian
        )
        field = factory.create_field()
        coords = field.get_dof_holder_coordinates()
        expected = (
            (grid.divisions[0] + 1)
            * (grid.divisions[1] + 1)
            * (grid.divisions[2] + 1)
        )
        assert coords.shape == (expected, 3)

        data = np.arange(expected, dtype=np.float64)
        expected_data = data.copy()
        field.set_dof_holder_data(data)
        data.fill(-999.0)
        np.testing.assert_allclose(field.get_dof_holder_data(), expected_data)

    def test_field_evaluation(self):
        """Evaluate a uniform-grid field at explicit query points."""
        grid = pcms.UniformGrid2D()
        grid.bot_left = [0.0, 0.0]
        grid.edge_length = [1.0, 1.0]
        grid.divisions = [10, 10]

        factory = pcms.LagrangeFunctionSpace.from_uniform_grid(
            grid, 1, pcms.CoordinateSystem.Cartesian
        )
        field = factory.create_field()
        coords = field.get_dof_holder_coordinates()
        num_dofs = field.get_num_dof_holders()
        data = np.zeros(num_dofs, dtype=np.float64)
        for i in range(num_dofs):
            data[i] = coords[i, 0] + 2.0 * coords[i, 1]
        field.set_dof_holder_data(data)

        query_pts = np.array(
            [[0.1, 0.2], [0.5, 0.5], [0.9, 0.8]], dtype=np.float64
        )
        request = pcms.EvaluationRequest.from_coordinates(query_pts)
        evaluator = factory.create_point_evaluator(request)
        results = np.zeros((len(query_pts), 1), dtype=np.float64)
        evaluator.evaluate(field, results)

        for i, pt in enumerate(query_pts):
            expected = pt[0] + 2.0 * pt[1]
            assert np.abs(results[i, 0] - expected) < 1e-6


class TestUniformGridOmegaHWorkflow:
    """Tests combining UniformGrid with Omega_h meshes (needs world)."""

    def test_uniform_grid_to_omega_h_workflow(self, world):
        """
        Test complete UniformGrid workflow with Omega_h mesh and field
        interpolation (Omega_h → UniformGrid).
        """
        # Create a simple 2D box mesh: 1.0 x 1.0 domain with 4x4 elements.
        mesh = omega_h.build_box(
            world, omega_h.Family.SIMPLEX, 1.0, 1.0, 0.0, 4, 4, 0, False,
        )

        grid = pcms.create_uniform_grid_from_mesh(mesh, [4, 4])
        mask_field = pcms.create_uniform_grid_binary_field(mesh, [4, 4])
        mask_data = mask_field.get_dof_holder_data()

        # Build a mesh-backed source field: f(x,y) = x + 2y.
        omega_h_factory = pcms.LagrangeFunctionSpace.from_mesh(mesh, 1)
        omega_h_field = omega_h_factory.create_field()
        coords = omega_h_field.get_dof_holder_coordinates()
        omega_h_data = np.zeros(omega_h_field.get_num_dof_holders(),
                                dtype=np.float64)
        for i in range(omega_h_field.get_num_dof_holders()):
            omega_h_data[i] = coords[i, 0] + 2.0 * coords[i, 1]
        omega_h_field.set_dof_holder_data(omega_h_data)

        ug_factory = pcms.LagrangeFunctionSpace.from_uniform_grid(
            grid, 1, pcms.CoordinateSystem.Cartesian
        )
        ug_field = ug_factory.create_field()

        # Transfer from unstructured mesh to uniform grid.
        interp = pcms.Interpolator(omega_h_factory, ug_factory)
        interp.apply(omega_h_field, ug_field)

        ug_field_data = ug_field.get_dof_holder_data()
        ug_coords = ug_field.get_dof_holder_coordinates()

        assert len(mask_data) == 25
        assert len(ug_field_data) == 25

        for vertex_id in range(len(ug_field_data)):
            x, y = ug_coords[vertex_id, 0], ug_coords[vertex_id, 1]
            expected = x + 2.0 * y
            actual = ug_field_data[vertex_id]
            assert abs(actual - expected) < 1e-10

    def test_omega_h_to_omega_h_transfer(self, world):
        """Transfer a field from one Omega_h mesh to another Omega_h mesh."""
        src_mesh = omega_h.build_box(
            world, omega_h.Family.SIMPLEX, 1.0, 1.0, 0.0, 4, 4, 0, False,
        )
        tgt_mesh = omega_h.build_box(
            world, omega_h.Family.SIMPLEX, 1.0, 1.0, 0.0, 8, 8, 0, False,
        )

        src_factory = pcms.LagrangeFunctionSpace.from_mesh(src_mesh, 1)
        tgt_factory = pcms.LagrangeFunctionSpace.from_mesh(tgt_mesh, 1)
        src_field = src_factory.create_field()
        tgt_field = tgt_factory.create_field()

        src_coords = src_field.get_dof_holder_coordinates()
        src_data = np.zeros(src_field.get_num_dof_holders())
        for i in range(len(src_data)):
            src_data[i] = src_coords[i, 0] + 2.0 * src_coords[i, 1]
        src_field.set_dof_holder_data(src_data)

        interp = pcms.Interpolator(src_factory, tgt_factory)
        interp.apply(src_field, tgt_field)

        tgt_coords = tgt_field.get_dof_holder_coordinates()
        tgt_data = tgt_field.get_dof_holder_data()
        for i in range(tgt_field.get_num_dof_holders()):
            expected = tgt_coords[i, 0] + 2.0 * tgt_coords[i, 1]
            assert abs(expected - tgt_data[i]) < 1e-10

