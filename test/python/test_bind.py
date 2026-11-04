"""
Tests for binding sources to query points from Python.
"""

import gc
import subprocess
import sys
import textwrap
import weakref

import numpy as np
import pytest
import pcms
import PyOmega_h as omega_h


def _linear(x, y):
    return 1.0 + 2.0 * x + 3.0 * y


@pytest.fixture
def mesh(world):
    return omega_h.build_box(
        world, omega_h.Family.SIMPLEX, 1.0, 1.0, 0.0, 8, 8, 0, False
    )


@pytest.fixture
def geometry(mesh):
    return pcms.backend.omega_h.create_geometry(mesh)


@pytest.fixture
def space(geometry):
    return pcms.create_lagrange_space(
        geometry, 1, 1, backend=pcms.LagrangeBackend.OmegaH)


@pytest.fixture
def field(space):
    f = space.create_function()
    xy = f.get_layout().get_dof_holder_coordinates().get_values()
    f.set_dof_holder_data(_linear(xy[:, 0], xy[:, 1]))
    return f


def test_evaluators_share_a_bound_source(space, field):
    # the second point is outside the unit square
    bound = space.get_geometry().bind(
        pcms.CoordinateView([[0.25, 0.5], [2.0, 2.0], [0.75, 0.1]]),
        pcms.OutOfBoundsPolicy(pcms.OutOfBoundsMode.FILL, -1.0))
    assert bound.num_points() == 3
    np.testing.assert_array_equal(bound.found_points(), [0, 2])
    np.testing.assert_array_equal(bound.missing_points(), [1])
    assert not bound.found_points().flags.writeable

    first = space.create_point_evaluator(bound)
    second = space.create_point_evaluator(bound)
    for evaluator in (first, second):
        out = np.zeros(3)
        evaluator.evaluate(field, out)
        np.testing.assert_allclose(
            out, [_linear(0.25, 0.5), -1.0, _linear(0.75, 0.1)], atol=1e-12)


def test_points_outside_raise_in_error_mode(space):
    with pytest.raises(RuntimeError, match="outside the source's domain"):
        space.get_geometry().bind(pcms.CoordinateView([[2.0, 2.0]]))


def test_bound_source_belongs_to_its_geometry(mesh, space):
    other = pcms.create_lagrange_space(
        pcms.backend.omega_h.create_geometry(mesh), 1, 1,
        backend=pcms.LagrangeBackend.OmegaH)
    bound = space.get_geometry().bind(pcms.CoordinateView([[0.5, 0.5]]))
    with pytest.raises(RuntimeError, match="another geometry"):
        other.create_point_evaluator(bound)


def test_spaces_on_one_geometry_share_a_bind(mesh):
    geometry = pcms.backend.omega_h.create_geometry(mesh)
    backend = pcms.LagrangeBackend.OmegaH
    scalars = pcms.create_lagrange_space(
        geometry, 1, 1, backend=backend)
    vectors = pcms.create_lagrange_space(
        geometry, 1, 3, backend=backend)
    bound = geometry.bind(pcms.CoordinateView([[0.25, 0.5], [0.75, 0.1]]))
    scalar = scalars.create_function()
    xy = scalar.get_layout().get_dof_holder_coordinates().get_values()
    scalar.set_dof_holder_data(_linear(xy[:, 0], xy[:, 1]))
    # three scalar components per vertex
    vector = vectors.create_function()
    vector.set_dof_holder_data(
        np.stack([_linear(xy[:, 0], xy[:, 1])] * 3, axis=1))
    expected = [_linear(0.25, 0.5), _linear(0.75, 0.1)]

    out = np.zeros(2)
    scalars.create_point_evaluator(bound).evaluate(scalar, out)
    np.testing.assert_allclose(out, expected, atol=1e-12)
    out3 = np.zeros((2, 3))
    vectors.create_point_evaluator(bound).evaluate(vector, out3)
    np.testing.assert_allclose(out3, np.stack([expected] * 3, axis=1),
                               atol=1e-12)


def test_derivatives_need_a_source_that_has_them(space):
    bound = space.get_geometry().bind(pcms.CoordinateView([[0.5, 0.5]]))
    assert space.supports(pcms.DifferentialOperator.Value)
    assert not space.supports(pcms.DifferentialOperator.Curl)
    with pytest.raises(RuntimeError, match="Gradient is not supported"):
        space.create_point_evaluator(bound, pcms.DifferentialOperator.Gradient)
    with pytest.raises(RuntimeError,
                       match="Apply: Gradient is not supported by the "
                             "wrapped source"):
        pcms.gradient(space)


def test_interpolator_reuses_a_bound_source(geometry, space, field):
    target = pcms.create_lagrange_space(
        geometry, 1, 1, backend=pcms.LagrangeBackend.OmegaH)
    bound = geometry.bind(target)
    out = target.create_function()
    pcms.Interpolator(space, target, bound).apply(field, out)
    np.testing.assert_allclose(out.get_dof_holder_data().get_values(),
                               field.get_dof_holder_data().get_values(),
                               atol=1e-12)


def test_get_geometry_outlives_the_original_wrappers(world):
    mesh = omega_h.build_box(
        world, omega_h.Family.SIMPLEX, 1.0, 1.0, 0.0, 4, 4, 0, False)
    geometry = pcms.backend.omega_h.create_geometry(mesh)
    space = pcms.create_lagrange_space(
        geometry, 1, 1, backend=pcms.LagrangeBackend.OmegaH)
    returned = space.get_geometry()
    assert returned == geometry
    space_ref = weakref.ref(space)
    geometry_ref = weakref.ref(geometry)
    mesh_ref = weakref.ref(mesh)

    del space, geometry, mesh
    gc.collect()
    assert space_ref() is None
    assert geometry_ref() is None
    assert mesh_ref() is None
    bound = returned.bind(pcms.CoordinateView([[0.5, 0.5]]))
    assert bound.num_points() == 1


@pytest.mark.parametrize("get_layout", [
    lambda space, field: space.get_layout(),
    lambda space, field: field.get_layout(),
], ids=["space", "field"])
def test_returned_layout_outlives_the_original_wrappers(world, get_layout):
    mesh = omega_h.build_box(
        world, omega_h.Family.SIMPLEX, 1.0, 1.0, 0.0, 4, 4, 0, False)
    geometry = pcms.backend.omega_h.create_geometry(mesh)
    space = pcms.create_lagrange_space(
        geometry, 1, 1, backend=pcms.LagrangeBackend.OmegaH)
    field = space.create_function()
    layout = get_layout(space, field)
    space_ref = weakref.ref(space)
    field_ref = weakref.ref(field)
    mesh_ref = weakref.ref(mesh)
    del field, space, geometry, mesh
    gc.collect()
    assert field_ref() is None
    assert space_ref() is None
    assert mesh_ref() is None
    coords = layout.get_dof_holder_coordinates().get_values()
    assert coords.shape == (25, 2)
    cloud = pcms.create_point_cloud_geometry(layout)
    reconstruction = pcms.create_polynomial_reconstruction_space(cloud)
    bound = cloud.bind(pcms.CoordinateView([[0.5, 0.5]]))
    np.testing.assert_array_equal(bound.found_points(), [0])
    reconstruction.create_point_evaluator(bound)


def test_function_get_space_owns_mesh_and_runtime():
    # A fresh process gives this test sole ownership of the Omega_h runtime.
    script = textwrap.dedent("""
        import gc
        import weakref
        import pcms
        import PyOmega_h as omega_h

        library = omega_h.OmegaHLibrary()
        world = library.world()
        mesh = omega_h.build_box(
            world, omega_h.Family.SIMPLEX, 1.0, 1.0, 0.0, 4, 4, 0, False)
        geometry = pcms.backend.omega_h.create_geometry(mesh)
        space = pcms.create_lagrange_space(
            geometry, 1, 1, backend=pcms.LagrangeBackend.OmegaH)
        field = space.create_function()
        returned = field.get_space()
        originals = [weakref.ref(obj) for obj in
                     (field, space, geometry, mesh, world, library)]
        del field, space, geometry, mesh, world, library
        gc.collect()
        assert all(ref() is None for ref in originals)
        replacement = returned.create_function()
        assert replacement.get_layout().get_num_dof_holder() == 25
        bound = returned.get_geometry().bind(pcms.CoordinateView([[0.5, 0.5]]))
        assert bound.found_points().tolist() == [0]
        del bound, replacement, returned
        gc.collect()
    """)
    subprocess.run([sys.executable, "-c", script], check=True, timeout=60)


def test_query_bound_outlives_the_original_wrappers(world):
    mesh = omega_h.build_box(
        world, omega_h.Family.SIMPLEX, 1.0, 1.0, 0.0, 4, 4, 0, False)
    geometry = pcms.backend.omega_h.create_geometry(mesh)
    query = pcms.QueryPoints(geometry)
    query.locate(pcms.CoordinateView([[0.5, 0.5]]))
    bound = query.bound
    mesh_ref = weakref.ref(mesh)
    del query, geometry, mesh
    gc.collect()
    assert mesh_ref() is None
    np.testing.assert_array_equal(bound.found_points(), [0])
    new_bound = bound.get_geometry().bind(pcms.CoordinateView([[0.25, 0.25]]))
    np.testing.assert_array_equal(new_bound.found_points(), [0])


def test_query_points_follow_moving_points(geometry, field):
    query = pcms.QueryPoints(
        geometry, pcms.OutOfBoundsPolicy(pcms.OutOfBoundsMode.FILL, -1.0))
    assert query.bound is None
    out = np.zeros(2)
    for points in ([[0.25, 0.5], [2.0, 2.0]], [[0.6, 0.6], [0.3, 0.2]]):
        query.locate(pcms.CoordinateView(points))
        query.evaluate(field, out)
        expected = [_linear(x, y) if x <= 1.0 else -1.0 for x, y in points]
        np.testing.assert_allclose(out, expected, atol=1e-12)
    np.testing.assert_array_equal(query.bound.found_points(), [0, 1])


def test_evaluate_at_points_in_one_call(field):
    out = np.zeros(2)
    pcms.evaluate(field, pcms.CoordinateView([[0.25, 0.5], [0.75, 0.1]]), out)
    np.testing.assert_allclose(out, [_linear(0.25, 0.5), _linear(0.75, 0.1)],
                               atol=1e-12)
