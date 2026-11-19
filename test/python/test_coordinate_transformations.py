"""
Tests for coordinate systems, value bases, and coordinate transformations.
"""

import gc

import numpy as np
import pytest
import pcms
import PyOmega_h as omega_h

# This just checks if we can get data without copying it. Once
# https://github.com/SCOREC/pcms/issues/361 is resolved, we should
# be able to have zero copy on both host and device side
ZERO_COPY = pcms.device_memory_is_host_accessible


def _cylindrical_vector():
    return pcms.ValueBasis(pcms.csys.create_cylindrical_r_theta_z(),
                           pcms.ComponentScaling.Physical, pcms.values.Vector)


def _cartesian_vector():
    return pcms.ValueBasis(pcms.csys.create_cartesian(3),
                           pcms.ComponentScaling.Physical, pcms.values.Vector)


def _coordinates(field):
    return field.get_layout().get_dof_holder_coordinates().get_values()


@pytest.fixture
def cylindrical_space(world):
    """P1 vector space on a box whose coordinates are (r, theta, z)."""
    mesh = omega_h.build_box(world, omega_h.Family.SIMPLEX, 2.0, 1.6, 1.1,
                             6, 6, 4, False)
    return pcms.create_lagrange_space(
        pcms.omega_h.create_geometry(
            mesh, pcms.csys.create_cylindrical_r_theta_z()),
        1, 3, backend=pcms.LagrangeBackend.OmegaH)


@pytest.fixture
def cartesian_space(world):
    """P1 vector space on the Cartesian unit cube."""
    mesh = omega_h.build_box(world, omega_h.Family.SIMPLEX, 1.0, 1.0, 1.0,
                             3, 3, 3, False)
    return pcms.create_lagrange_space(
        pcms.omega_h.create_geometry(mesh, pcms.csys.create_cartesian(3)),
        1, 3, backend=pcms.LagrangeBackend.OmegaH)


def _coordinates_of(space):
    return space.get_layout().get_dof_holder_coordinates().get_values()


def _radial_field(space):
    """B = r e_r + z e_z in physical cylindrical components, which is
    (x, y, z) in Cartesian ones."""
    b = space.create_function("b", pcms.values.Vector,
                              scaling=pcms.ComponentScaling.Physical)
    p = _coordinates(b)
    b.set_dof_holder_data(pcms.ValueView(
        b.get_data().get_value_basis(),
        np.stack([p[:, 0], np.zeros(len(p)), p[:, 2]], axis=1)))
    return b


class TestCoordinateSystem:
    def test_kinds(self):
        cases = [
            (pcms.csys.create_cartesian(2), "cartesian", 2),
            (pcms.csys.create_cylindrical_r_theta_z(), "cylindrical_rthetaz", 3),
            (pcms.csys.create_cylindrical_rz(), "cylindrical_rz", 2),
        ]
        for system, kind, dim in cases:
            assert isinstance(system, pcms.CoordinateSystem)
            assert system.kind() == kind
            assert system.dimension() == dim

    def test_properties(self):
        cyl = pcms.csys.create_cylindrical_r_theta_z()
        assert cyl.kind() == "cylindrical_rthetaz"
        assert cyl.has_orthogonal_basis()
        assert not cyl.has_unit_scale_factors()
        cart = pcms.csys.create_cartesian(3)
        assert cart.has_orthogonal_basis() and cart.has_unit_scale_factors()

    def test_equality_is_by_value(self):
        a = pcms.csys.create_cartesian(3)
        assert a == pcms.csys.create_cartesian(3)
        assert a != pcms.csys.create_cartesian(2)
        assert a != pcms.csys.create_cylindrical_r_theta_z()
        assert len({a, pcms.csys.create_cartesian(3)}) == 1
        assert a != "cartesian"

    def test_space_defaults_to_cartesian_of_mesh_dimension(self, world):
        mesh = omega_h.build_box(world, omega_h.Family.SIMPLEX, 1.0, 1.0, 0.0,
                                 4, 4, 0, False)
        space = pcms.create_lagrange_space(
            pcms.omega_h.create_geometry(mesh), 1, 1)
        assert (space.get_layout().get_coordinate_system()
                == pcms.csys.create_cartesian(2))

    def test_create_rejects_bad_dimension(self):
        with pytest.raises(RuntimeError):
            pcms.csys.create_cartesian(4)


class TestCoordinateView:
    def test_defaults_to_cartesian_of_array_dimension(self):
        view = pcms.CoordinateView(np.zeros((4, 2)))
        assert view.get_coordinate_system() == pcms.csys.create_cartesian(2)

    def test_rejects_wrong_dimension(self):
        with pytest.raises(RuntimeError, match="coordinate columns"):
            pcms.CoordinateView(pcms.csys.create_cylindrical_r_theta_z(),
                                np.zeros((4, 2)))

    def test_layout_coordinates_carry_the_space_system(self,
                                                       cylindrical_space):
        coords = cylindrical_space.get_layout().get_dof_holder_coordinates()
        assert (coords.get_coordinate_system()
                == pcms.csys.create_cylindrical_r_theta_z())
        assert not coords.get_values().flags.writeable
        assert cylindrical_space.get_coordinate_system() == \
            pcms.csys.create_cylindrical_r_theta_z()


class TestValueBasis:
    def test_variance_signatures(self):
        assert pcms.values.Scalar.rank() == 0
        assert list(pcms.values.Vector) == [pcms.Variance.Contravariant]
        assert list(pcms.values.Covector) == [pcms.Variance.Covariant]
        assert pcms.values.Tensor.rank() == 2
        assert pcms.values.Of([pcms.Variance.Covariant]) == \
            pcms.values.Covector
        assert pcms.VarianceSignature(
            [pcms.Variance.Contravariant]) == pcms.values.Vector
        with pytest.raises(ValueError):
            pcms.values.Of([pcms.Variance.Covariant] * 5)

    def test_equality_is_by_value(self):
        # Every scalar basis is the same, whatever its members say.
        assert pcms.ValueBasis() == pcms.ValueBasis(
            pcms.csys.create_cartesian(3))
        assert _cylindrical_vector() == _cylindrical_vector()
        assert _cylindrical_vector() != _cartesian_vector()

    def test_members(self):
        basis = _cylindrical_vector()
        assert basis.rank() == 1
        assert basis.system == pcms.csys.create_cylindrical_r_theta_z()
        assert basis.scaling == pcms.ComponentScaling.Physical
        assert basis.variance == pcms.values.Vector


class TestValueView:
    def test_checks_component_count(self):
        with pytest.raises(RuntimeError, match="require 3 components"):
            pcms.ValueView(_cylindrical_vector(), np.zeros((4, 2)))

    def test_rejects_non_2d(self):
        with pytest.raises(ValueError):
            pcms.ValueView(pcms.ValueBasis(), np.zeros(4))

    def test_aliases_the_array(self):
        values = np.zeros((4, 3))
        view = pcms.ValueView(_cylindrical_vector(), values)
        assert view.get_values() is values
        assert view.rank() == 1
        assert view.get_value_type() == pcms.FieldValueType.Vector
        assert (view.extent(0), view.extent(1), view.size()) == (4, 3, 12)
        np.testing.assert_array_equal(np.asarray(view), values)


class TestFunction:
    def test_value_basis_defaults_to_the_space_system(self, cartesian_space):
        b = cartesian_space.create_function("b", pcms.values.Vector)
        assert b.get_name() == "b"
        assert b.get_data().get_value_basis() == _cartesian_vector()
        assert b.get_data().get_value_type() == pcms.FieldValueType.Vector

    def test_borrowed_basis(self, cylindrical_space):
        b = cylindrical_space.create_function(
            "b", pcms.values.Vector, system=pcms.csys.create_cartesian(3))
        assert b.get_data().get_value_basis() == _cartesian_vector()

    def test_ambiguous_scaling_must_be_stated(self, cylindrical_space):
        with pytest.raises(RuntimeError, match="scaling must be stated"):
            cylindrical_space.create_function("b", pcms.values.Vector)

    def test_get_space_returns_the_space(self, cartesian_space):
        f = cartesian_space.create_function()
        assert f.get_space() == cartesian_space
        assert isinstance(f, pcms.Field)

    def test_dof_holder_data_is_tagged(self, cylindrical_space):
        b = _radial_field(cylindrical_space)
        data = b.get_dof_holder_data()
        assert data.get_basis() == _cylindrical_vector()
        assert not data.get_values().flags.writeable
        np.testing.assert_array_equal(data.get_values()[:, 0],
                                      _coordinates(b)[:, 0])

    def test_set_checks_the_basis(self, cylindrical_space):
        b = _radial_field(cylindrical_space)
        n = b.get_layout().get_num_dof_holder()
        with pytest.raises(RuntimeError, match="does not match"):
            b.set_dof_holder_data(
                pcms.ValueView(_cartesian_vector(), np.zeros((n, 3))))

    def test_set_checks_the_shape(self, cylindrical_space):
        b = _radial_field(cylindrical_space)
        n = b.get_layout().get_num_dof_holder()
        with pytest.raises(ValueError):
            b.set_dof_holder_data(
                pcms.ValueView(_cylindrical_vector(), np.zeros((n - 1, 3))))

    def test_plain_arrays_only_for_scalars(self, cylindrical_space):
        b = _radial_field(cylindrical_space)
        n = b.get_layout().get_num_dof_holder()
        with pytest.raises(TypeError, match="need a ValueView"):
            b.set_dof_holder_data(np.zeros((n, 3)))
        b.set_dof_holder_data_unchecked(np.ones((n, 3)))
        np.testing.assert_array_equal(b.get_dof_holder_data().get_values(),
                                      1.0)
        with pytest.raises(ValueError):
            b.set_dof_holder_data_unchecked(np.ones((n, 2)))


class TestCoordinateMap:
    def test_bind_maps_points_and_rotates_values(self):
        cart = pcms.CoordinateView(pcms.csys.create_cartesian(3),
                                   [[1.0, 0.0, 0.5], [0.0, 2.0, -1.0],
                                    [1.0, 1.0, 0.0]])
        m = pcms.create_cartesian_to_cylindrical()
        assert m.get_source_coordinate_system() == \
            pcms.csys.create_cartesian(3)
        assert m.get_target_coordinate_system() == \
            pcms.csys.create_cylindrical_r_theta_z()
        bound = m.bind(cart)
        assert bound.num_points() == 3
        mapped = bound.mapped_points()
        assert mapped.get_coordinate_system() == \
            pcms.csys.create_cylindrical_r_theta_z()
        rtz = mapped.get_values()
        np.testing.assert_allclose(rtz[:, 0], [1.0, 2.0, np.sqrt(2.0)])
        np.testing.assert_allclose(rtz[:, 1], [0.0, np.pi / 2, np.pi / 4])

        assert bound.output_basis(_cylindrical_vector()) == \
            _cartesian_vector()
        e_theta = np.tile([0.0, 1.0, 0.0], (3, 1))
        out = bound.transform_values(
            pcms.ValueView(_cylindrical_vector(), e_theta))
        assert out.get_basis() == _cartesian_vector()
        assert out.get_values() is e_theta  # transformed in place
        np.testing.assert_allclose(
            e_theta, np.stack([-np.sin(rtz[:, 1]), np.cos(rtz[:, 1]),
                               np.zeros(3)], axis=1), atol=1e-15)

    def test_round_trip(self):
        rtz = pcms.CoordinateView(pcms.csys.create_cylindrical_r_theta_z(),
                                  [[1.0, 0.3, 0.5], [2.0, 1.2, -1.0]])
        xyz = pcms.create_cylindrical_to_cartesian().bind(rtz).mapped_points()
        back = pcms.create_cartesian_to_cylindrical().bind(xyz).mapped_points()
        np.testing.assert_allclose(back.get_values(), rtz.get_values())

    def test_bind_checks_the_source_system(self):
        rtz = pcms.CoordinateView(pcms.csys.create_cylindrical_r_theta_z(),
                                  np.ones((2, 3)))
        with pytest.raises(RuntimeError):
            pcms.create_cartesian_to_cylindrical().bind(rtz)

    def test_transform_values_requires_a_writeable_buffer(self):
        cart = pcms.CoordinateView(pcms.csys.create_cartesian(3),
                                   np.ones((2, 3)))
        bound = pcms.create_cartesian_to_cylindrical().bind(cart)
        values = np.zeros((2, 3))
        values.flags.writeable = False
        with pytest.raises(ValueError, match="writeable"):
            bound.transform_values(
                pcms.ValueView(_cylindrical_vector(), values))

    def test_mapped_points_outlive_the_bound_map(self):
        cart = pcms.CoordinateView(pcms.csys.create_cartesian(3),
                                   [[0.0, 3.0, 1.0]])
        points = pcms.create_cartesian_to_cylindrical().bind(cart).mapped_points()
        gc.collect()
        np.testing.assert_allclose(points.get_values(),
                                   [[3.0, np.pi / 2, 1.0]])


def _bind_through_map(cylindrical_space, target, policy=None):
    """Bind `target` (a space or a CoordinateView in Cartesian coordinates)
    on the cylindrical space's geometry through CartesianToCylindrical."""
    if policy is None:
        policy = pcms.OutOfBoundsPolicy()
    return cylindrical_space.get_geometry().bind(
        target, policy, pcms.create_cartesian_to_cylindrical())


class TestBindWithMap:
    def _evaluator(self, cylindrical_space, cartesian_space):
        bound = _bind_through_map(cylindrical_space, cartesian_space)
        assert isinstance(bound, pcms.BoundPoints)
        assert bound.num_points() == len(_coordinates_of(cartesian_space))
        return cylindrical_space.create_point_evaluator(bound)

    def test_evaluate_rotates_into_the_query_basis(self, cylindrical_space,
                                                   cartesian_space):
        b = _radial_field(cylindrical_space)
        evaluator = self._evaluator(cylindrical_space, cartesian_space)
        basis = evaluator.output_basis(b.get_data().get_value_basis())
        assert basis == _cartesian_vector()
        n = evaluator.num_points()
        out = np.zeros((n, 3))
        evaluator.evaluate(b, pcms.ValueView(basis, out))
        np.testing.assert_allclose(out, _coordinates_of(cartesian_space),
                                   atol=1e-12)
        assert len(evaluator.filled_points()) == 0

    def test_evaluate_checks_its_output(self, cylindrical_space,
                                        cartesian_space):
        b = _radial_field(cylindrical_space)
        evaluator = self._evaluator(cylindrical_space, cartesian_space)
        n = evaluator.num_points()
        with pytest.raises(TypeError, match="need a ValueView"):
            evaluator.evaluate(b, np.zeros((n, 3)))
        with pytest.raises(RuntimeError, match="does not match"):
            evaluator.evaluate(
                b, pcms.ValueView(_cylindrical_vector(), np.zeros((n, 3))))
        with pytest.raises(ValueError, match="must be"):
            evaluator.evaluate(
                b, pcms.ValueView(_cartesian_vector(), np.zeros((n + 1, 3))))

    def test_query_points_must_be_in_the_geometry_system(self,
                                                         cylindrical_space):
        coords = pcms.CoordinateView(pcms.csys.create_cartesian(3),
                                     np.full((2, 3), 0.5))
        with pytest.raises(RuntimeError, match="as Bind's map argument"):
            cylindrical_space.get_geometry().bind(coords)

    def test_the_map_must_end_in_the_geometry_system(self, cylindrical_space):
        coords = pcms.CoordinateView(pcms.csys.create_cartesian(3),
                                     np.full((2, 3), 0.5))
        with pytest.raises(RuntimeError,
                           match="must take the query coordinate system"):
            cylindrical_space.get_geometry().bind(
                coords, map=pcms.create_cylindrical_to_cartesian())

    def test_filled_points(self, cylindrical_space):
        # The second point has theta = pi, outside the source's [0, 1.6].
        bound = _bind_through_map(
            cylindrical_space,
            pcms.CoordinateView(pcms.csys.create_cartesian(3),
                                [[0.5, 0.5, 0.5], [-0.5, 0.0, 0.5]]),
            pcms.OutOfBoundsPolicy(pcms.OutOfBoundsMode.FILL, -7.0))
        np.testing.assert_array_equal(bound.found_points(), [0])
        np.testing.assert_array_equal(bound.missing_points(), [1])
        evaluator = cylindrical_space.create_point_evaluator(bound)
        filled = evaluator.filled_points()
        np.testing.assert_array_equal(filled, [1])
        assert not filled.flags.writeable
        b = _radial_field(cylindrical_space)
        out = np.zeros((2, 3))
        evaluator.evaluate(b, pcms.ValueView(_cartesian_vector(), out))
        np.testing.assert_allclose(out[0], [0.5, 0.5, 0.5], atol=1e-12)
        np.testing.assert_array_equal(out[1], -7.0)

    def test_interpolator_transfers_across_systems(self, cylindrical_space,
                                                   cartesian_space):
        b_src = _radial_field(cylindrical_space)
        b_tgt = cartesian_space.create_function("b", pcms.values.Vector)
        op = pcms.create_interpolator(
            cylindrical_space, cartesian_space,
            _bind_through_map(cylindrical_space, cartesian_space))
        assert isinstance(op, pcms.TransferOperator)
        op.apply(b_src, b_tgt)
        np.testing.assert_allclose(b_tgt.get_dof_holder_data().get_values(),
                                   _coordinates(b_tgt), atol=1e-12)

    def test_interpolator_checks_the_target_basis(self, cylindrical_space,
                                                  cartesian_space):
        b_src = _radial_field(cylindrical_space)
        wrong = cartesian_space.create_function(
            "b", pcms.values.Vector,
            system=pcms.csys.create_cylindrical_r_theta_z(),
            scaling=pcms.ComponentScaling.Physical)
        op = pcms.create_interpolator(
            cylindrical_space, cartesian_space,
            _bind_through_map(cylindrical_space, cartesian_space))
        with pytest.raises(RuntimeError, match="target's declared basis"):
            op.apply(b_src, wrong)

    def test_plain_interpolator_refuses_to_bridge_systems(
            self, cylindrical_space, cartesian_space):
        with pytest.raises(RuntimeError, match="as Bind's map argument"):
            pcms.create_interpolator(cylindrical_space, cartesian_space)


class TestLifetimes:
    def test_space_keeps_its_mesh_alive(self, world):
        mesh = omega_h.build_box(world, omega_h.Family.SIMPLEX, 1.0, 1.0, 0.0,
                                 4, 4, 0, False)
        space = pcms.create_lagrange_space(
            pcms.omega_h.create_geometry(mesh), 1, 1,
            backend=pcms.LagrangeBackend.OmegaH)
        del mesh
        gc.collect()
        f = space.create_function()
        coords = _coordinates(f)
        f.set_dof_holder_data(coords[:, :1])
        evaluator = space.create_point_evaluator(
            space.get_geometry().bind(pcms.CoordinateView([[0.3, 0.6]])))
        out = np.zeros((1, 1))
        evaluator.evaluate(f, out)
        np.testing.assert_allclose(out, [[0.3]], atol=1e-12)

    def test_operators_keep_their_sources_alive(self, cylindrical_space,
                                                cartesian_space):
        b_src = _radial_field(cylindrical_space)
        b_tgt = cartesian_space.create_function("b", pcms.values.Vector)
        op = pcms.create_interpolator(
            cylindrical_space, cartesian_space,
            _bind_through_map(cylindrical_space, cartesian_space))
        # the BoundPoints owns its bound map; the map object may go
        cart_to_cyl = pcms.create_cartesian_to_cylindrical()
        bound = cylindrical_space.get_geometry().bind(
            cartesian_space, map=cart_to_cyl)
        evaluator = cylindrical_space.create_point_evaluator(bound)
        del cart_to_cyl, bound
        gc.collect()
        op.apply(b_src, b_tgt)
        out = np.zeros((evaluator.num_points(), 3))
        evaluator.evaluate(b_src, pcms.ValueView(_cartesian_vector(), out))
        np.testing.assert_allclose(out, b_tgt.get_dof_holder_data(),
                                   atol=1e-12)
