"""
Tests for MLS-backed polynomial reconstruction spaces.
"""

import numpy as np
try:
    import PyOmega_h as omega_h
except ImportError:  # pcms built without Omega_h
    omega_h = None
import pytest
import pcms

pytestmark = pytest.mark.skipif(
    not pcms.features["polynomial_reconstruction"],
    reason="pcms was built without POLYNOMIAL_RECONSTRUCTION")


def poly_value(x, y, degree):
    """Evaluate a simple polynomial at (x, y)."""
    if degree == 0:
        return 3.0
    if degree == 1:
        return x + y
    if degree == 2:
        return x**2 + y**2
    if degree == 3:
        return x**3 + y**3
    raise ValueError(f"Unsupported polynomial degree: {degree}")


class TestMLSInterpolation:
    """Tests for MLS polynomial reproduction."""

    def test_polynomial_reproduction(self):
        """
        Verify that MLS reproduces polynomials up to the configured degree.

        Source points are a 4x4 grid on [0,1]^2.  Target points are six
        arbitrary interior locations.  For each interpolation degree d, MLS
        must exactly reproduce any polynomial of degree <= d (up to the
        configured tolerance).
        """
        tolerance = 5e-4

        grid_vals = np.linspace(0.0, 1.0, 4)
        source_xy = np.array(
            [(x, y) for x in grid_vals for y in grid_vals],
            dtype=np.float64,
        )
        target_xy = np.array(
            [
                [0.11, 0.19],
                [0.23, 0.77],
                [0.51, 0.49],
                [0.87, 0.26],
                [0.69, 0.91],
                [0.35, 0.62],
            ],
            dtype=np.float64,
        )
        num_targets = target_xy.shape[0]

        for interp_degree in range(1, 4):
            opts = pcms.MLSOptions()
            opts.degree = interp_degree
            opts.basis = pcms.RadialBasisFunction.NO_OP
            neighbourhood = pcms.NeighbourhoodOptions()
            neighbourhood.radius = 1.5
            neighbourhood.adapt_radius = False

            factory = pcms.create_polynomial_reconstruction_space(
                pcms.point_cloud.create_geometry(
                    pcms.CoordinateView(source_xy), neighbourhood),
                opts,
            )
            field = factory.create_function()
            evaluator = factory.create_point_evaluator(
                factory.get_geometry().bind(pcms.CoordinateView(target_xy)))
            results = np.zeros((num_targets, 1), dtype=np.float64)

            for func_degree in range(interp_degree, -1, -1):
                source_values = np.array(
                    [poly_value(x, y, func_degree) for (x, y) in source_xy],
                    dtype=np.float64,
                )
                exact_target_values = np.array(
                    [poly_value(x, y, func_degree) for (x, y) in target_xy],
                    dtype=np.float64,
                )

                field.set_dof_holder_data(source_values)
                evaluator.evaluate(field, results)
                approx_target_values = results[:, 0]

                max_abs_err = np.max(
                    np.abs(exact_target_values - approx_target_values)
                )
                assert max_abs_err < tolerance, (
                    f"MLS failed for interp_degree={interp_degree}, "
                    f"func_degree={func_degree}: max_abs_err={max_abs_err}"
                )


class TestMLSOnMeshEntities:
    """MLS reconstruction from the entities of an Omega_h mesh."""

    def test_element_centroids_to_vertices(self, world):
        """
        Reconstruct a linear function from the element centroids of a box
        mesh at its vertices; the supports come from the mesh adjacency.
        """
        mesh = omega_h.build_box(
            world, omega_h.Family.SIMPLEX, 1.0, 1.0, 0.0, 4, 4, 0, False)
        face_dim = 2
        mesh_geometry = pcms.omega_h.create_geometry(mesh)
        geometry = pcms.omega_h.create_point_cloud_geometry(mesh_geometry, face_dim)
        opts = pcms.MLSOptions()
        opts.degree = 1
        opts.basis = pcms.RadialBasisFunction.NO_OP
        space = pcms.create_polynomial_reconstruction_space(geometry, opts)

        centroids = space.get_layout().get_dof_holder_coordinates().get_values()
        assert centroids.shape == (32, 2)
        field = space.create_function()
        field.set_dof_holder_data(
            (1.0 + 2.0 * centroids[:, 0] - 3.0 * centroids[:, 1]).reshape(-1, 1))

        # a space on the same geometry object: binding the cloud at its DOF
        # holders forwards the hints of that geometry to the mesh search
        vertex_space = pcms.create_lagrange_space(mesh_geometry, 1, 1)
        vertices = np.array(
            vertex_space.get_layout().get_dof_holder_coordinates().get_values())
        evaluator = space.create_point_evaluator(geometry.bind(vertex_space))
        results = np.zeros((vertices.shape[0], 1), dtype=np.float64)
        evaluator.evaluate(field, results)

        np.testing.assert_allclose(
            results[:, 0], 1.0 + 2.0 * vertices[:, 0] - 3.0 * vertices[:, 1],
            atol=1e-4)

    def test_entity_dim_out_of_range(self, world):
        mesh = omega_h.build_box(
            world, omega_h.Family.SIMPLEX, 1.0, 1.0, 0.0, 2, 2, 0, False)
        with pytest.raises(ValueError):
            pcms.omega_h.create_point_cloud_geometry(
                pcms.omega_h.create_geometry(mesh),
                3)

    def test_geometry_must_be_an_omega_h_geometry(self):
        grid = pcms.UniformGrid2D()
        grid.bot_left = [0.0, 0.0]
        grid.edge_length = [1.0, 1.0]
        grid.divisions = [2, 2]
        with pytest.raises(ValueError, match="Omega_h geometry"):
            pcms.omega_h.create_point_cloud_geometry(
                pcms.uniform_grid.create_geometry(grid),
                0)
