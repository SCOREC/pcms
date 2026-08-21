"""
Tests for MLS-backed PolynomialReconstructionFunctionSpace.
"""

import numpy as np
import pcms


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
            opts.radius = 1.5
            opts.adapt_radius = False
            opts.basis = pcms.RadialBasisFunction.NO_OP

            factory = pcms.PolynomialReconstructionFunctionSpace.from_coords(
                source_xy, pcms.CoordinateSystem.Cartesian, opts
            )
            field = factory.create_field()
            request = pcms.EvaluationRequest.from_coordinates(target_xy)
            evaluator = factory.create_point_evaluator(request)
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
