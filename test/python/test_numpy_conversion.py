"""
Tests for numpy <-> pcms array conversions.
"""

import gc

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
def p1_space(mesh):
    return pcms.LagrangeFunctionSpace.from_mesh(
        mesh, 1, 1, pcms.CoordinateSystem.Cartesian
    )


def _linear_field(space):
    field = space.create_field()
    coords = field.get_dof_holder_coordinates()
    field.set_dof_holder_data(_linear(coords[:, 0], coords[:, 1]))
    return field


class TestFieldData:
    def test_dof_holder_coordinates_match_mesh(self, mesh, p1_space):
        field = p1_space.create_field()
        coords = field.get_dof_holder_coordinates()
        assert coords.shape == (mesh.nverts(), 2)
        assert coords.flags.c_contiguous
        np.testing.assert_array_equal(coords, mesh.coords().reshape(-1, 2))

    def test_set_from_strided_array(self, p1_space):
        field = p1_space.create_field()
        n = field.get_num_dof_holders()
        wide = np.arange(2 * n, dtype=np.float64).reshape(n, 2)
        field.set_dof_holder_data(wide[:, 1])  # stride of two doubles
        np.testing.assert_array_equal(field.get_dof_holder_data(), wide[:, 1])

    def test_set_from_float32_array(self, p1_space):
        field = p1_space.create_field()
        n = field.get_num_dof_holders()
        data = np.arange(n, dtype=np.float32)
        field.set_dof_holder_data(data)
        np.testing.assert_array_equal(field.get_dof_holder_data(),
                                      data.astype(np.float64))

    def test_set_rejects_2d_array(self, p1_space):
        field = p1_space.create_field()
        n = field.get_num_dof_holders()
        with pytest.raises(ValueError):
            field.set_dof_holder_data(np.zeros((n, 1)))


class TestPointEvaluation:
    QUERY = np.array([[0.5, 0.5], [0.25, 0.75], [0.9, 0.1]])

    def _evaluator(self, space, coords):
        request = pcms.EvaluationRequest.from_coordinates(coords)
        return space.create_point_evaluator(request)

    def test_evaluate_linear_field(self, p1_space):
        field = _linear_field(p1_space)
        out = np.zeros((len(self.QUERY), 1))
        self._evaluator(p1_space, self.QUERY).evaluate(field, out)
        np.testing.assert_allclose(
            out[:, 0], _linear(self.QUERY[:, 0], self.QUERY[:, 1]),
            atol=1e-12)

    def test_request_from_transposed_coordinates(self, p1_space):
        field = _linear_field(p1_space)
        coords = np.ascontiguousarray(self.QUERY.T).T  # same values, F-order
        assert not coords.flags.c_contiguous
        out = np.zeros((len(coords), 1))
        self._evaluator(p1_space, coords).evaluate(field, out)
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
    def test_owns_a_copy_of_the_coordinates(self):
        coords = np.arange(12, dtype=np.float64).reshape(6, 2)
        expected = coords.copy()
        # A strided input forces a temporary copy.
        # the view must not refer to it, or to the 
        #caller's array, after construction.
        view = pcms.CoordinateView(pcms.CoordinateSystem.Cartesian,
                                   np.asfortranarray(coords))
        coords[:] = -1.0
        del coords
        gc.collect()
        np.testing.assert_array_equal(view.get_coordinates(), expected)


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
        np.testing.assert_allclose(result, 3.5)


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
