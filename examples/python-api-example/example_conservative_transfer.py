#!/usr/bin/env python3
"""
Conservative field transfer between two Omega_h meshes.

This example:
  1. Builds a source Omega_h 2D simplex mesh, writes it to disk, and reads it
     back (demonstrating the mesh read path).
  2. Evaluates an analytic function onto an order-1 (linear) Lagrange field on
     the source mesh.
  3. Builds, writes, and reads back a *different* target Omega_h mesh that
     carries its own order-1 Lagrange field.
  4. Transfers the source field onto the target field with both the Monte
     Carlo (control-variate) and the mesh-intersection conservative projection
     methods, and reports the L2 / max errors and the conservation error
     (preservation of the integral) for each.
  5. Demonstrates the mesh-intersection projection across element orders as
     well: P1 -> P0 (projection to cell averages) and P0 -> P1, confirming the
     integral is conserved in both directions.

The mesh-intersection conservative projection supports scalar Lagrange spaces
of order 0 (piecewise constant, one DOF per element) or order 1 (piecewise
linear, one DOF per vertex) on 2D simplex meshes; source and target orders may
differ. The Monte Carlo operator requires order-1 spaces.
"""

import argparse
import os
import gc
import tempfile

import numpy as np
import pcms


# Smooth analytic field sampled onto the source mesh.
def source_function(x, y):
    return np.sin(np.pi * x) * np.sin(np.pi * y) + 0.5 * x + 0.25 * y


def build_and_roundtrip_mesh(world, lib, nx, ny, path):
    """Build a 2D box mesh, write it to disk, and read it back."""
    mesh = pcms.build_box(
        world,
        pcms.Family.SIMPLEX,
        1.0, 1.0, 0.0,   # unit square (z = 0 -> 2D)
        nx, ny, 0,
        False,
    )
    pcms.write_mesh_binary(path, mesh)
    del mesh
    gc.collect()
    return pcms.read_mesh_binary(path, lib)


def integral_of_field(mesh, field):
    """Integral of an order-1 nodal field via element-averaged quadrature."""
    values = field.get_dof_holder_data()
    elem_verts = mesh.ask_elem_verts().reshape(-1, 3)
    areas = mesh.ask_sizes()  # element areas for a 2D simplex mesh
    elem_avg = values[elem_verts].mean(axis=1)
    return float(np.sum(areas * elem_avg))


def integral_of_p0_field(mesh, field):
    """Integral of an order-0 field: sum of per-element value * element area."""
    values = field.get_dof_holder_data()  # one value per element
    areas = mesh.ask_sizes()
    return float(np.sum(areas * values))


def evaluate_function_onto_field(space, func):
    """Create a field on `space` and fill it from `func` at the DOF holders."""
    field = space.create_field()
    coords = field.get_dof_holder_coordinates()
    field.set_dof_holder_data(func(coords[:, 0], coords[:, 1]))
    return field


def report(name, target_field, target_mesh, reference_coords,
           reference_values, source_integral, integrate=integral_of_field):
    values = target_field.get_dof_holder_data()
    err = values - reference_values
    l2 = np.sqrt(np.mean(err ** 2))
    linf = np.max(np.abs(err))
    target_integral = integrate(target_mesh, target_field)
    cons = abs(target_integral - source_integral)
    print(f"  {name}")
    print(f"    L2 error vs analytic        : {l2:.3e}")
    print(f"    max error vs analytic       : {linf:.3e}")
    print(f"    target integral             : {target_integral:.6f}")
    print(f"    conservation error |dI|     : {cons:.3e}")


def _triangulation(mesh, vertex_coords):
    """matplotlib Triangulation from an Omega_h 2D simplex mesh."""
    import matplotlib.tri as mtri
    tris = mesh.ask_elem_verts().reshape(-1, 3)
    return mtri.Triangulation(vertex_coords[:, 0], vertex_coords[:, 1], tris)


def _draw_field(ax, triang, field, order, title, vmin, vmax):
    values = field.get_dof_holder_data()
    if order == 1:
        # Nodal (per-vertex) values -> smooth (Gouraud) shading.
        art = ax.tripcolor(triang, values, shading="gouraud",
                           vmin=vmin, vmax=vmax)
    else:
        # Cell (per-element) values -> flat shading.
        art = ax.tripcolor(triang, facecolors=values, shading="flat",
                          vmin=vmin, vmax=vmax)
    ax.set_title(title, fontsize=9)
    ax.set_aspect("equal")
    ax.set_xticks([])
    ax.set_yticks([])
    return art


def make_plots(out_path, source_mesh, target_mesh, source_vertex_coords,
               target_vertex_coords, panels):
    """Render before/after field panels to `out_path` (PNG).

    `panels` is a list of (mesh_key, field, order, title) where mesh_key is
    "source" or "target". Returns True on success, False if matplotlib is
    unavailable.
    """
    try:
        import matplotlib
        matplotlib.use("Agg")  # file output, no display needed
        import matplotlib.pyplot as plt
    except ImportError:
        print("\n[plot] matplotlib not available; skipping plots. "
              "Install matplotlib to enable --plot.")
        return False

    triangs = {
        "source": _triangulation(source_mesh, source_vertex_coords),
        "target": _triangulation(target_mesh, target_vertex_coords),
    }

    # Shared color scale across panels so fields are directly comparable.
    all_vals = np.concatenate([f.get_dof_holder_data() for _, f, _, _ in panels])
    vmin, vmax = float(all_vals.min()), float(all_vals.max())

    ncols = 3
    nrows = (len(panels) + ncols - 1) // ncols
    fig, axes = plt.subplots(nrows, ncols, figsize=(4 * ncols, 3.6 * nrows),
                             squeeze=False)
    art = None
    for idx, (mesh_key, field, order, title) in enumerate(panels):
        ax = axes[idx // ncols][idx % ncols]
        art = _draw_field(ax, triangs[mesh_key], field, order, title, vmin, vmax)
    for idx in range(len(panels), nrows * ncols):
        axes[idx // ncols][idx % ncols].axis("off")

    fig.colorbar(art, ax=axes, shrink=0.85, label="field value")
    fig.suptitle("Conservative projection: before / after", fontsize=12)
    fig.savefig(out_path, dpi=130, bbox_inches="tight")
    plt.close(fig)
    print(f"\n[plot] wrote {out_path}")
    return True


def main(plot_path=None):
    lib = pcms.OmegaHLibrary()
    world = lib.world()

    work_dir = tempfile.mkdtemp(prefix="pcms_conservative_transfer_")
    try:
        print("=" * 64)
        print("Conservative field transfer: MonteCarlo and mesh intersection")
        print("=" * 64)

        # 1. Source mesh (read from disk) + linear Lagrange field.
        source_path = os.path.join(work_dir, "source.osh")
        source_mesh = build_and_roundtrip_mesh(world, lib, 12, 12, source_path)
        print(f"\nSource mesh : {source_mesh.nverts()} verts, "
              f"{source_mesh.nelems()} elems (read from {source_path})")

        # The conservative-projection operators require the native Omega_h
        # Lagrange backend (not the default MeshFields layout).
        source_space = pcms.LagrangeFunctionSpace.from_mesh(
            source_mesh, 1, 1, pcms.CoordinateSystem.Cartesian(),
            backend=pcms.LagrangeFunctionSpace.Backend.OmegaH,
        )
        source_field = evaluate_function_onto_field(source_space,
                                                    source_function)
        source_integral = integral_of_field(source_mesh, source_field)
        print(f"Source field: order-1 Lagrange, "
              f"integral = {source_integral:.6f}")

        # 2. Target mesh (read from disk) with its own linear Lagrange field.
        target_path = os.path.join(work_dir, "target.osh")
        target_mesh = build_and_roundtrip_mesh(world, lib, 17, 17, target_path)
        print(f"\nTarget mesh : {target_mesh.nverts()} verts, "
              f"{target_mesh.nelems()} elems (read from {target_path})")

        target_space = pcms.LagrangeFunctionSpace.from_mesh(
            target_mesh, 1, 1, pcms.CoordinateSystem.Cartesian(),
            backend=pcms.LagrangeFunctionSpace.Backend.OmegaH,
        )

        # Analytic reference values at the target DOF holders.
        target_template = target_space.create_field()
        target_coords = target_template.get_dof_holder_coordinates()
        reference_values = source_function(target_coords[:, 0],
                                           target_coords[:, 1])

        print("\nTransferring source field to target field:")

        # 3a. Mesh-intersection conservative projection.
        mi_target = target_space.create_field()
        mesh_intersection = pcms.OmegaHConservativeProjection(
            source_space, target_space
        )
        mesh_intersection.apply(source_field, mi_target)
        report("mesh intersection (exact quadrature)", mi_target, target_mesh,
               target_coords, reference_values, source_integral)

        # 3b. Mesh-intersection projection with a row-sum lumped mass matrix.
        #     The linear solve degenerates to an exact diagonal scaling; the
        #     integral is still conserved exactly, but the max nodal error
        #     grows relative to the consistent mass matrix.
        lumped_target = target_space.create_field()
        lumped_projection = pcms.OmegaHConservativeProjection(
            source_space, target_space,
            mass_matrix_type=pcms.MassMatrixType.Lumped,
        )
        lumped_projection.apply(source_field, lumped_target)
        report("mesh intersection (lumped mass)", lumped_target, target_mesh,
               target_coords, reference_values, source_integral)

        # 3c. Monte Carlo (control-variate) conservative projection.
        mc_target = target_space.create_field()
        monte_carlo = pcms.OmegaHControlVariateProjection(
            source_space,
            target_space,
            samples_per_element=64,
            sampling=pcms.MonteCarloSampling.UniformRandom,
            seed=8675309,
        )
        monte_carlo.apply(source_field, mc_target)
        report("Monte Carlo (control variate, UniformRandom)", mc_target, target_mesh,
               target_coords, reference_values, source_integral)

        # 4. Mixed-order mesh-intersection projections. The same operator works
        #    for any combination of order-0 (piecewise constant) and order-1
        #    (piecewise linear) Lagrange spaces; the integral is conserved in
        #    every case.
        print("\nMixed-order mesh-intersection projections:")

        # 4a. P1 source -> P0 target (project the linear field to cell averages).
        target_p0_space = pcms.LagrangeFunctionSpace.from_mesh(
            target_mesh, 0, 1, pcms.CoordinateSystem.Cartesian(),
            backend=pcms.LagrangeFunctionSpace.Backend.OmegaH,
        )
        p1_to_p0 = target_p0_space.create_field()
        # Analytic reference at the P0 DOF holders (element centroids).
        p0_coords = p1_to_p0.get_dof_holder_coordinates()
        p0_reference = source_function(p0_coords[:, 0], p0_coords[:, 1])
        proj_p1_p0 = pcms.OmegaHConservativeProjection(
            source_space, target_p0_space
        )
        proj_p1_p0.apply(source_field, p1_to_p0)
        report("P1 -> P0 (cell averages)", p1_to_p0, target_mesh,
               p0_coords, p0_reference, source_integral,
               integrate=integral_of_p0_field)

        # 4b. P0 source -> P1 target. Build a P0 source field by sampling the
        #     analytic function at the source element centroids.
        source_p0_space = pcms.LagrangeFunctionSpace.from_mesh(
            source_mesh, 0, 1, pcms.CoordinateSystem.Cartesian(),
            backend=pcms.LagrangeFunctionSpace.Backend.OmegaH,
        )
        source_p0_field = evaluate_function_onto_field(source_p0_space,
                                                       source_function)
        source_p0_integral = integral_of_p0_field(source_mesh, source_p0_field)
        p0_to_p1 = target_space.create_field()
        proj_p0_p1 = pcms.OmegaHConservativeProjection(
            source_p0_space, target_space
        )
        proj_p0_p1.apply(source_p0_field, p0_to_p1)
        report("P0 -> P1", p0_to_p1, target_mesh,
               target_coords, reference_values, source_p0_integral)

        print("\n✓ Field transfer completed for P1<->P1, P1->P0, and P0->P1.")

        # 5. Optional before/after visualization.
        if plot_path is not None:
            source_vertex_coords = source_field.get_dof_holder_coordinates()
            make_plots(
                plot_path, source_mesh, target_mesh,
                source_vertex_coords, target_coords,
                panels=[
                    ("source", source_field, 1, "source (P1)"),
                    ("target", mi_target, 1, "P1 -> P1 (mesh intersection)"),
                    ("target", p1_to_p0, 0, "P1 -> P0 (cell averages)"),
                    ("source", source_p0_field, 0, "source (P0)"),
                    ("target", p0_to_p1, 1, "P0 -> P1"),
                ],
            )

    finally:
        import shutil
        shutil.rmtree(work_dir, ignore_errors=True)
        #del world
        #del lib


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--plot", nargs="?", const="conservative_transfer.png", default=None,
        metavar="PATH",
        help="render before/after field panels to PATH (default "
             "conservative_transfer.png) using matplotlib, if installed",
    )
    args = parser.parse_args()
    main(plot_path=args.plot)
