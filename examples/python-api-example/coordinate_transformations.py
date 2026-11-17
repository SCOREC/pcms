#!/usr/bin/env python3
"""Coordinate transformations in the PCMS Python API.

  (a) map points and re-express vector components with a CoordinateMap
  (b) transfer where the source points are cylindrical but both fields store
      Cartesian components: only the points are mapped
  (c) transfer where the source points and components are both cylindrical
      and the target's are Cartesian: points are mapped, values rotated

Each case prints the maximum error against the exact answer. This is the
Python version of examples/coordinate-transformation-example.
"""
import sys

import numpy as np
import pcms
import PyOmega_h as omega_h


def cylindrical_vector():
    return pcms.ValueBasis(pcms.csys.CylindricalRThetaZ.create(),
                           pcms.ComponentScaling.Physical, pcms.values.Vector)


def fill_field(field, f):
    """Set a 3-component field from f(p) at its DOF holders, where p are the
    holder coordinates (num_dof_holders × 3) in the field's space."""
    p = field.get_layout().get_dof_holder_coordinates().get_values()
    field.set_dof_holder_data(
        pcms.ValueView(field.get_data().get_value_basis(), f(p)))


def max_error_from_xyz(b):
    """Max over the field's DOF holders of |b - (x, y, z)|."""
    x = b.get_layout().get_dof_holder_coordinates().get_values()
    return np.max(np.abs(b.get_dof_holder_data().get_values() - x))


def build_spaces(world):
    """Source: a box whose coordinates are (r, theta, z) with r in [0, 2] and
    theta in [0, 1.6]. Target: the Cartesian unit cube, which lies inside
    it."""
    src_mesh = omega_h.build_box(world, omega_h.Family.SIMPLEX, 2.0, 1.6, 1.1,
                                 8, 8, 6, False)
    src = pcms.LagrangeFunctionSpace.from_mesh(
        src_mesh, 1, 3, pcms.csys.CylindricalRThetaZ.create(),
        backend=pcms.LagrangeFunctionSpace.Backend.OmegaH)
    tgt_mesh = omega_h.build_box(world, omega_h.Family.SIMPLEX, 1.0, 1.0, 1.0,
                                 4, 4, 4, False)
    tgt = pcms.LagrangeFunctionSpace.from_mesh(
        tgt_mesh, 1, 3, pcms.csys.Cartesian.create(3),
        backend=pcms.LagrangeFunctionSpace.Backend.OmegaH)
    return src, tgt


def case_a():
    print("(a) coordinate transformation")
    xyz = pcms.CoordinateView(pcms.csys.Cartesian.create(3),
                              [[1.0, 0.0, 0.5], [0.0, 2.0, -1.0],
                               [1.0, 1.0, 0.0]])
    bound = pcms.CartesianToCylindrical().bind(xyz)
    rtz = bound.mapped_points().get_values()

    # Values run opposite to points: this map re-expresses cylindrical
    # components as Cartesian ones, in place. Here, e_theta at each point.
    v = np.tile([0.0, 1.0, 0.0], (3, 1))
    v_cart = bound.transform_values(pcms.ValueView(cylindrical_vector(), v))
    theta = rtz[:, 1]
    err = max(np.max(np.abs(v[:, 0] + np.sin(theta))),
              np.max(np.abs(v[:, 1] - np.cos(theta))))
    for p, e in zip(rtz, v):
        print("    (r,th,z)=(%.3f, %.3f, %.3f)  e_theta=(%.3f, %.3f, %.3f)"
              % (*p, *e))
    print("    output basis '%s', max |e_theta - (-sin, cos, 0)| = %.2e"
          % (v_cart.get_basis().system.kind(), err))
    return err


def case_b(world):
    print("(b) cylindrical points, Cartesian components on both sides")
    src, tgt = build_spaces(world)

    # The space's points are (r, theta, z), but the field's
    # components are Cartesian. B = (x, y, z) = (r cos theta, r sin theta, z).
    b_src = src.create_function("b", pcms.values.Vector,
                                system=pcms.csys.Cartesian.create(3))
    fill_field(b_src, lambda p: np.stack(
        [p[:, 0] * np.cos(p[:, 1]), p[:, 0] * np.sin(p[:, 1]), p[:, 2]],
        axis=1))
    b_tgt = tgt.create_function("b", pcms.values.Vector)

    # The values are already in the target's basis, so the map passes them
    # through and only the points are mapped.
    src_as_cartesian = pcms.TransformedSource(src,
                                              pcms.CartesianToCylindrical())
    pcms.Interpolator(src_as_cartesian, tgt).apply(b_src, b_tgt)
    # B's Cartesian components are not linear in (r, theta, z), so
    # the P1 source field only approximates them.
    err = max_error_from_xyz(b_tgt)
    print("    via Interpolator:  max |B - (x,y,z)| = %.2e" % err)
    return err


def case_c(world):
    print("(c) cylindrical points and components -> Cartesian")
    src, tgt = build_spaces(world)

    # B = r e_r + z e_z in physical cylindrical components, which is (x, y, z)
    # in Cartesian ones. Its cylindrical components are linear in
    # (r, theta, z), so the transfer is exact.
    b_src = src.create_function("b", pcms.values.Vector,
                                scaling=pcms.ComponentScaling.Physical)
    fill_field(b_src, lambda p: np.stack(
        [p[:, 0], np.zeros(len(p)), p[:, 2]], axis=1))
    b_tgt = tgt.create_function("b", pcms.values.Vector)

    # `to_source` carries the target's
    # Cartesian points to (r, theta, z), and the source's cylindrical
    # components back to Cartesian.
    src_as_cartesian = pcms.TransformedSource(src,
                                              pcms.CartesianToCylindrical())
    pcms.Interpolator(src_as_cartesian, tgt).apply(b_src, b_tgt)
    err = max_error_from_xyz(b_tgt)
    print("    max |B - (x,y,z)| = %.2e" % err)
    return err


def main():
    lib = omega_h.OmegaHLibrary()
    world = lib.world()
    errors = [case_a(), case_b(world), case_c(world)]
    # (b) is a P1 approximation of a nonlinear field; (a) and (c) are exact.
    ok = errors[0] < 1e-12 and errors[1] < 2e-2 and errors[2] < 1e-12
    print("SUCCESS" if ok else "FAILED: error above tolerance")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
