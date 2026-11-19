"""
Objects created from others hold what they need: the inputs can go out of
scope (and be collected) before the result is used.
"""

import gc

import numpy as np
import pytest
import pcms
import PyOmega_h as omega_h


def _linear(x, y):
    return 1.0 + 2.0 * x + 3.0 * y


def _mesh(world, n=8):
    return omega_h.build_box(
        world, omega_h.Family.SIMPLEX, 1.0, 1.0, 0.0, n, n, 0, False)


def _p1(mesh, **options):
    return pcms.create_lagrange_space(pcms.omega_h.create_geometry(mesh), 1, 1,
                                      **options)


def _grid_space():
    grid = pcms.UniformGrid2D()
    grid.bot_left = [0.0, 0.0]
    grid.edge_length = [1.0, 1.0]
    grid.divisions = [4, 4]
    return pcms.create_lagrange_space(pcms.uniform_grid.create_geometry(grid),
                                      1, 1)


def _linear_function(space):
    f = space.create_function()
    xy = f.get_layout().get_dof_holder_coordinates().get_values()
    f.set_dof_holder_data(_linear(xy[:, 0], xy[:, 1]))
    return f


def _assert_linear(field, atol=1e-12):
    xy = field.get_layout().get_dof_holder_coordinates().get_values()
    np.testing.assert_allclose(field.get_dof_holder_data().get_values()[:, 0],
                               _linear(xy[:, 0], xy[:, 1]), atol=atol)


def test_point_evaluator_outlives_its_source_and_bound(world):
    mesh = _mesh(world)
    space = _p1(mesh)
    f = _linear_function(space)
    points = [[0.25, 0.5], [0.75, 0.1]]
    bound = space.get_geometry().bind(pcms.CoordinateView(points))
    evaluator = space.create_point_evaluator(bound)
    del mesh, space, bound
    gc.collect()
    out = np.zeros(2)
    evaluator.evaluate(f, out)
    np.testing.assert_allclose(out, [_linear(*p) for p in points], atol=1e-12)


def test_point_cloud_geometry_outlives_the_mesh_geometry(world):
    mesh = _mesh(world)
    geometry = pcms.omega_h.create_geometry(mesh)
    cloud = pcms.omega_h.create_point_cloud_geometry(geometry, 2)
    del mesh, geometry
    gc.collect()
    opts = pcms.MLSOptions()
    opts.degree = 1
    opts.basis = pcms.RadialBasisFunction.NO_OP
    f = _linear_function(pcms.create_polynomial_reconstruction_space(cloud,
                                                                     opts))
    points = [[0.3, 0.6], [0.55, 0.45]]
    out = np.zeros(2)
    pcms.evaluate(f, pcms.CoordinateView(points), out)
    np.testing.assert_allclose(out, [_linear(*p) for p in points], atol=1e-4)


def test_interpolator_outlives_its_spaces(world):
    mesh = _mesh(world)
    source = _p1(mesh)
    target = _grid_space()
    f = _linear_function(source)
    out = target.create_function()
    op = pcms.create_interpolator(source, target)
    del mesh, source, target
    gc.collect()
    op.apply(f, out)
    _assert_linear(out)


def test_interpolator_with_bound_outlives_its_inputs(world):
    mesh = _mesh(world)
    source = _p1(mesh)
    target = _grid_space()
    f = _linear_function(source)
    out = target.create_function()
    bound = source.get_geometry().bind(target)
    op = pcms.create_interpolator(source, target, bound)
    del mesh, source, target, bound
    gc.collect()
    op.apply(f, out)
    _assert_linear(out)


@pytest.mark.skipif(not hasattr(pcms, "create_conservative_projection"),
                    reason="needs PETSc and MeshFields")
def test_conservative_projection_outlives_its_spaces(world):
    source_mesh = _mesh(world, 6)
    target_mesh = _mesh(world, 5)
    # the projections need Omega_h layouts, not the MeshFields default
    source = _p1(source_mesh, backend=pcms.LagrangeBackend.OmegaH)
    target = _p1(target_mesh, backend=pcms.LagrangeBackend.OmegaH)
    f = _linear_function(source)
    out = target.create_function()
    op = pcms.create_conservative_projection(source, target)
    del source_mesh, target_mesh, source, target
    gc.collect()
    op.apply(f, out)
    _assert_linear(out, atol=1e-10)


@pytest.mark.skipif(not hasattr(pcms, "create_control_variate_projection"),
                    reason="needs PETSc and MeshFields")
def test_control_variate_projection_outlives_its_inputs(world):
    source_mesh = _mesh(world, 6)
    target_mesh = _mesh(world, 5)
    # the projections need Omega_h layouts, not the MeshFields default
    source = _p1(source_mesh, backend=pcms.LagrangeBackend.OmegaH)
    target = _p1(target_mesh, backend=pcms.LagrangeBackend.OmegaH)
    f = _linear_function(source)
    out = target.create_function()
    op = pcms.create_control_variate_projection(source, target, 4)
    del source_mesh, target_mesh, source, target
    gc.collect()
    op.apply(f, out)
    _assert_linear(out, atol=1e-10)
