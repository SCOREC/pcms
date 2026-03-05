#!/usr/bin/env python3

import numpy as np
import pcms


def poly_value(x, y, degree):
    if degree == 0:
        return 3.0
    if degree == 1:
        return x + y
    if degree == 2:
        return x**2 + y**2
    if degree == 3:
        return x**3 + y**3
    raise ValueError(f"Unsupported polynomial degree: {degree}")


def build_full_support(num_targets, num_sources):
    supports_ptr = np.arange(0, (num_targets + 1) * num_sources,
                             num_sources, dtype=np.int32)
    supports_idx = np.tile(np.arange(num_sources, dtype=np.int32), num_targets)
    radii2 = np.full(num_targets, 10.0, dtype=np.float64)
    return pcms.SupportResults(supports_ptr, supports_idx, radii2)


def test_mls_interpolation_polynomial_reproduction():
    tolerance = 5e-4

    grid_vals = np.linspace(0.0, 1.0, 4)
    source_xy = np.array([(x, y) for x in grid_vals for y in grid_vals],
                         dtype=np.float64)
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

    num_sources = source_xy.shape[0]
    num_targets = target_xy.shape[0]

    source_coordinates = source_xy.reshape(-1)
    target_coordinates = target_xy.reshape(-1)
    support = build_full_support(num_targets, num_sources)

    print("Testing MLS interpolation for polynomial reproduction...")

    for interp_degree in range(1, 4):
        for func_degree in range(interp_degree, -1, -1):
            source_values = np.array(
                [poly_value(x, y, func_degree) for (x, y) in source_xy],
                dtype=np.float64,
            )
            exact_target_values = np.array(
                [poly_value(x, y, func_degree) for (x, y) in target_xy],
                dtype=np.float64,
            )

            approx_target_values = pcms.mls_interpolation(
                source_values,
                source_coordinates,
                target_coordinates,
                support,
                2,
                interp_degree,
                pcms.RadialBasisFunction.NO_OP,
                1e-5,
                1e-6,
                5.0,
            )

            max_abs_err = np.max(np.abs(exact_target_values - approx_target_values))
            assert max_abs_err < tolerance, (
                f"MLS failed for interp_degree={interp_degree}, "
                f"func_degree={func_degree}: max_abs_err={max_abs_err}"
            )

if __name__ == "__main__":
    lib = pcms.OmegaHLibrary()
    world = lib.world()
    
    try:
        test_mls_interpolation_polynomial_reproduction()
        print("MLS interpolation test passed")
    except Exception as e:
        print(f"MLS interpolation test failed: {e}")
        import traceback
        traceback.print_exc()
        exit(1)
    finally:
        # Explicitly delete objects to avoid an MPI finalizing issue
        del world
        del lib
