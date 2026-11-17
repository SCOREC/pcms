"""
Tests for numpy <-> pcms array conversions.
"""

import gc

import numpy as np
import pytest
import pcms
import PyOmega_h as omega_h

ZERO_COPY = pcms.device_memory_is_host_accessible


def _linear(x, y):
    return 1.0 + 2.0 * x + 3.0 * y


@pytest.fixture
def mesh(world):
    return omega_h.build_box(
        world, omega_h.Family.SIMPLEX, 1.0, 1.0, 0.0, 8, 8, 0, False
    )


@pytest.fixture
def p1_space(mesh):
    return pcms.LagrangeFunctionSpace.from_mesh(
        mesh, 1, 1)


def _coordinates(field):
    return field.get_layout().get_dof_holder_coordinates().get_values()


def _num_dof_holders(field):
    return field.get_layout().get_num_owned_dof_holder()


def _linear_field(space):
    field = space.create_function()
    coords = _coordinates(field)
    field.set_dof_holder_data(_linear(coords[:, 0], coords[:, 1]))
    return field


class TestFieldData:
    def test_shape_matches_layout(self, mesh, world):
        space = pcms.LagrangeFunctionSpace.from_mesh(
            mesh, 1, 2)
        field = space.create_function()
        n = _num_dof_holders(field)
        data = np.arange(2 * n, dtype=np.float64).reshape(n, 2)
        field.set_dof_holder_data(data)
        result = field.get_dof_holder_data().get_values()
        assert result.shape == (n, 2)
        np.testing.assert_array_equal(result, data)

    def test_dof_holder_coordinates_match_mesh(self, mesh, p1_space):
        field = p1_space.create_function()
        coords = _coordinates(field)
        assert coords.shape == (mesh.nverts(), 2)
        np.testing.assert_array_equal(coords, mesh.coords().reshape(-1, 2))

    def test_set_from_strided_array(self, p1_space):
        field = p1_space.create_function()
        n = _num_dof_holders(field)
        wide = np.arange(2 * n, dtype=np.float64).reshape(n, 2)
        field.set_dof_holder_data(wide[:, 1:])  # stride of two doubles
        np.testing.assert_array_equal(field.get_dof_holder_data(), wide[:, 1:])

    def test_set_from_float32_array(self, p1_space):
        field = p1_space.create_function()
        n = _num_dof_holders(field)
        data = np.arange(n, dtype=np.float32).reshape(n, 1)
        field.set_dof_holder_data(data)
        np.testing.assert_array_equal(field.get_dof_holder_data(),
                                      data.astype(np.float64))

    def test_scalar_field_accepts_1d(self, p1_space):
        field = p1_space.create_function()
        n = _num_dof_holders(field)
        data = np.arange(n, dtype=np.float64)
        field.set_dof_holder_data(data)
        result = field.get_dof_holder_data().get_values()
        assert result.shape == (n, 1)
        np.testing.assert_array_equal(result[:, 0], data)

    @pytest.mark.parametrize("shape", [
        lambda n: (n, 2), lambda n: (n + 1, 1), lambda n: (n + 1,),
        lambda n: (n, 1, 1),
    ], ids=["wrong-components", "wrong-dof-holders", "1d-wrong-dof-holders",
            "3d"])
    def test_set_rejects_wrong_shape(self, p1_space, shape):
        field = p1_space.create_function()
        with pytest.raises(ValueError):
            field.set_dof_holder_data(np.zeros(shape(_num_dof_holders(field))))

    def test_vector_field_rejects_1d(self, mesh):
        space = pcms.LagrangeFunctionSpace.from_mesh(
            mesh, 1, 2)
        field = space.create_function()
        n = _num_dof_holders(field)
        for shape in [(n,), (2 * n,)]:
            with pytest.raises(ValueError):
                field.set_dof_holder_data(np.zeros(shape))

    def test_returned_arrays_are_read_only(self, p1_space):
        field = _linear_field(p1_space)
        for arr in (field.get_dof_holder_data().get_values(),
                    _coordinates(field)):
            assert not arr.flags.writeable
            with pytest.raises(ValueError):
                arr[0, 0] = 0.0

    def test_data_is_zero_copy_on_host_builds(self, p1_space):
        field = _linear_field(p1_space)
        assert np.shares_memory(field.get_dof_holder_data().get_values(),
                                field.get_dof_holder_data().get_values()) == \
            ZERO_COPY

    def test_data_reflects_set(self, p1_space):
        # Only arrays obtained after the modification are guaranteed current.
        field = _linear_field(p1_space)
        new = field.get_dof_holder_data().get_values() + 1.0
        field.set_dof_holder_data(new)
        np.testing.assert_array_equal(field.get_dof_holder_data(), new)

    def test_coordinates_are_zero_copy_on_host_builds(self, p1_space):
        field = p1_space.create_function()
        assert np.shares_memory(_coordinates(field), _coordinates(field)) == \
            ZERO_COPY

    def test_shared_arrays_keep_field_alive(self, p1_space):
        field = _linear_field(p1_space)
        data = field.get_dof_holder_data().get_values()
        coords = _coordinates(field)
        expected_data, expected_coords = data.copy(), coords.copy()
        del field
        gc.collect()
        np.testing.assert_array_equal(data, expected_data)
        np.testing.assert_array_equal(coords, expected_coords)


class TestPointEvaluation:
    QUERY = np.array([[0.5, 0.5], [0.25, 0.75], [0.9, 0.1]])

    def _evaluator(self, space, coords):
        request = pcms.EvaluationRequest.from_coordinates(
            pcms.CoordinateView(coords))
        return space.create_point_evaluator(request)

    def test_evaluate_linear_field(self, p1_space):
        field = _linear_field(p1_space)
        out = np.zeros((len(self.QUERY), 1))
        self._evaluator(p1_space, self.QUERY).evaluate(field, out)
        np.testing.assert_allclose(
            out[:, 0], _linear(self.QUERY[:, 0], self.QUERY[:, 1]),
            atol=1e-12)

    def test_evaluate_into_1d_output(self, p1_space):
        field = _linear_field(p1_space)
        out = np.zeros(len(self.QUERY))
        self._evaluator(p1_space, self.QUERY).evaluate(field, out)
        np.testing.assert_allclose(
            out, _linear(self.QUERY[:, 0], self.QUERY[:, 1]), atol=1e-12)

    def test_evaluate_rejects_wrong_components(self, p1_space):
        field = _linear_field(p1_space)
        out = np.zeros((len(self.QUERY), 2))
        with pytest.raises(ValueError):
            self._evaluator(p1_space, self.QUERY).evaluate(field, out)

    def test_request_from_transposed_coordinates(self, p1_space):
        field = _linear_field(p1_space)
        coords = np.ascontiguousarray(self.QUERY.T).T  # same values, F-order
        assert not coords.flags.c_contiguous
        out = np.zeros((len(coords), 1))
        self._evaluator(p1_space, coords).evaluate(field, out)
        np.testing.assert_allclose(
            out[:, 0], _linear(self.QUERY[:, 0], self.QUERY[:, 1]),
            atol=1e-12)

    def test_request_copies_coordinates(self, p1_space):
        field = _linear_field(p1_space)
        coords = self.QUERY.copy()
        evaluator = self._evaluator(p1_space, coords)
        coords[:] = -1.0  # must not affect the request
        del coords
        gc.collect()
        out = np.zeros((len(self.QUERY), 1))
        evaluator.evaluate(field, out)
        np.testing.assert_allclose(
            out[:, 0], _linear(self.QUERY[:, 0], self.QUERY[:, 1]),
            atol=1e-12)

    @pytest.mark.parametrize("make_output", [
        lambda n: np.zeros((n, 1), dtype=np.float32),
        lambda n: np.zeros((n, 2))[:, :1],
    ], ids=["float32", "strided"])
    def test_evaluate_rejects_unwritable_output_layout(self, p1_space,
                                                       make_output):
        field = _linear_field(p1_space)
        out = make_output(len(self.QUERY))
        with pytest.raises(ValueError):
            self._evaluator(p1_space, self.QUERY).evaluate(field, out)

    def test_evaluate_rejects_read_only_output(self, p1_space):
        field = _linear_field(p1_space)
        out = np.zeros((len(self.QUERY), 1))
        out.flags.writeable = False
        with pytest.raises(ValueError):
            self._evaluator(p1_space, self.QUERY).evaluate(field, out)


class TestCoordinateView:
    def test_copies_input(self):
        coords = np.arange(12, dtype=np.float64).reshape(6, 2)
        expected = coords.copy()
        view = pcms.CoordinateView(coords)
        assert not np.shares_memory(view.get_values(), coords)
        coords[:] = -1.0
        del coords
        gc.collect()
        np.testing.assert_array_equal(view.get_values(), expected)

    def test_strided_input(self):
        coords = np.arange(12, dtype=np.float64).reshape(6, 2)
        view = pcms.CoordinateView(np.asfortranarray(coords))
        gc.collect()
        np.testing.assert_array_equal(view.get_values(), coords)
        assert not view.get_values().flags.writeable

    def test_coordinates_keep_view_alive(self):
        coords = np.arange(12, dtype=np.float64).reshape(6, 2)
        values = pcms.CoordinateView(coords.copy()).get_values()
        gc.collect()
        np.testing.assert_array_equal(values, coords)

    def test_default_system_is_cartesian_of_array_dimension(self):
        coords = np.zeros((4, 2))
        cartesian = pcms.csys.Cartesian.create(2)
        assert pcms.CoordinateView(coords).get_coordinate_system() == cartesian
        assert (pcms.CoordinateView(None, coords).get_coordinate_system()
                == cartesian)

    def test_system_must_match_array_dimension(self):
        with pytest.raises(Exception):
            pcms.CoordinateView(pcms.csys.CylindricalRThetaZ.create(),
                                np.zeros((4, 2)))


class TestMeshUtilities:
    def test_average_from_strided_input(self, mesh):
        nfaces = mesh.nfaces()
        values = np.arange(nfaces, dtype=np.float64)
        wide = np.stack([values, -values], axis=1)
        contiguous = pcms.map_entity_field_to_vertices_average(mesh, values, 2)
        strided = pcms.map_entity_field_to_vertices_average(mesh, wide[:, 0],
                                                            2)
        np.testing.assert_array_equal(strided, contiguous)

    def test_average_of_constant_is_constant(self, mesh):
        values = np.full(mesh.nfaces(), 3.5)
        result = pcms.map_entity_field_to_vertices_average(mesh, values, 2)
        assert result.shape == (mesh.nverts(),)
        assert not result.flags.writeable
        np.testing.assert_allclose(result, 3.5)

    @pytest.mark.parametrize("mapping", [
        lambda mesh, v: pcms.map_entity_field_to_vertices_average(mesh, v, 2),
        lambda mesh, v: pcms.map_entity_field_to_vertices_mls(mesh, v, 2),
    ], ids=["average", "mls"])
    def test_multicomponent_matches_per_component(self, mesh, mapping):
        centroids = pcms.compute_entity_centroids(mesh, 2)
        values = np.stack([_linear(centroids[:, 0], centroids[:, 1]),
                           centroids[:, 0] - centroids[:, 1]], axis=1)
        result = mapping(mesh, values)
        assert result.shape == (mesh.nverts(), 2)
        for c in range(2):
            np.testing.assert_allclose(
                result[:, c], mapping(mesh, values[:, c]), atol=1e-12)

    def test_single_component_2d_keeps_rank(self, mesh):
        values = np.full((mesh.nfaces(), 1), 3.5)
        result = pcms.map_entity_field_to_vertices_average(mesh, values, 2)
        assert result.shape == (mesh.nverts(), 1)
        np.testing.assert_allclose(result, 3.5)

    @pytest.mark.parametrize("shape", [
        lambda n: (n + 1,), lambda n: (n + 1, 2), lambda n: (n, 0),
        lambda n: (n, 1, 1),
    ], ids=["1d-wrong-entities", "2d-wrong-entities", "no-components", "3d"])
    def test_rejects_wrong_shape(self, mesh, shape):
        with pytest.raises(ValueError):
            pcms.map_entity_field_to_vertices_average(
                mesh, np.zeros(shape(mesh.nfaces())), 2)

    def test_centroids_survive_garbage_collection(self, mesh):
        centroids = pcms.compute_entity_centroids(mesh, 2)
        expected = centroids.copy()
        gc.collect()
        np.testing.assert_array_equal(centroids, expected)
        assert centroids.shape == (mesh.nfaces(), 2)


class TestUniformGrid:
    def test_closest_cell_id_accepts_sequences(self, mesh):
        grid = pcms.create_uniform_grid_from_mesh(mesh, [4, 4])
        expected = grid.closest_cell_id(np.array([0.3, 0.6]))
        assert grid.closest_cell_id([0.3, 0.6]) == expected
        assert grid.closest_cell_id(np.array([0.3, 0.6], np.float32)) == \
            expected

    def test_closest_cell_id_rejects_wrong_dimension(self, mesh):
        grid = pcms.create_uniform_grid_from_mesh(mesh, [4, 4])
        with pytest.raises(ValueError):
            grid.closest_cell_id(np.array([0.3, 0.6, 0.1]))
