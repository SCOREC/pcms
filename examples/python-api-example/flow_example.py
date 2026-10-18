#!/usr/bin/env python3
import pcms
import PyOmega_h as omega_h
import numpy as np
# import matplotlib.pyplot as plt


def print_mesh_info(mesh):
    print(f"\nDimension: {mesh.dim()}")
    print(f"Vertices: {mesh.nverts()}")
    print(f"Faces: {mesh.nfaces()}")
    print(f"Elements: {mesh.nelems()}")
    
    ent_names = ["VERTICES", "EDGES", "FACES"]
    
    for ent_dim in range(3):
        num_tags = mesh.ntags(ent_dim)
        print(f"\n{ent_names[ent_dim]}: {num_tags} tags")
        for i in range(num_tags):
            tag = mesh.get_tag(ent_dim, i)
            print(f"  [{i}] '{tag.name()}': {tag.ncomps()} comp(s)")


def demonstrate_face_field_transfer(mesh):
    print(f"\n--- Face Field Transfer ---")
    
    face_dim = 2
    face_tag = mesh.get_tag(face_dim, 14) # tag index 14
    face_field_values = mesh.get_tag(face_dim, face_tag.name())
    
    if face_tag.ncomps() > 1:
        face_field_values = face_field_values.reshape(-1, face_tag.ncomps())[:, 0]
    
    print(f"Face field: min={face_field_values.min():.6f}, max={face_field_values.max():.6f}, mean={face_field_values.mean():.6f}")
    
    # Use the C++ implementation from pcms Python API
    # vertex_field_values = pcms.map_entity_field_to_vertices_mls(
    #     mesh,
    #     face_field_values,
    #     face_dim,
    #     interp_degree=2,
    # )
    # Alternative method (simple averaging):
    vertex_field_values = pcms.map_entity_field_to_vertices_average(
        mesh, face_field_values, face_dim
    )
    print(
        f"Vertex field: min={vertex_field_values.min():.6f}, "
        f"max={vertex_field_values.max():.6f}, "
        f"mean={vertex_field_values.mean():.6f}"
    )
    
    omega_h_factory = pcms.LagrangeFunctionSpace.from_mesh(
        mesh, 1, 1, pcms.CoordinateSystem.Cartesian
    )
    omega_h_field = omega_h_factory.create_field()
    omega_h_field.set_dof_holder_data(vertex_field_values)
    
    divisions = [1000, 500]
    grid = pcms.create_uniform_grid_from_mesh(mesh, divisions)
    ug_factory = pcms.LagrangeFunctionSpace.from_uniform_grid(
        grid, 1, pcms.CoordinateSystem.Cartesian
    )
    ug_field = ug_factory.create_field()
    
    print(f"Interpolating field from OmegaH mesh to uniform grid...")
    interp = pcms.Interpolator(
        omega_h_factory,
        ug_factory,
        pcms.OutOfBoundsPolicy(pcms.OutOfBoundsMode.FILL, 0.0),
    )
    interp.apply(omega_h_field, ug_field)
    
    transferred_data = ug_field.get_dof_holder_data()[:, 0]
    print(f"Transferred field: min={np.min(transferred_data):.6f}, max={np.max(transferred_data):.6f}, mean={np.mean(transferred_data):.6f}")
    
    # Visualize the transferred data
    # nx, ny = divisions[0] + 1, divisions[1] + 1
    # grid_data = transferred_data.reshape(ny, nx)
    
    # plt.figure(figsize=(10, 8))
    # plt.imshow(grid_data, origin='lower', cmap='viridis', interpolation='nearest', vmin=0, vmax=500)
    # plt.colorbar(label='Field Value')
    # plt.title('Transferred Field Data on Uniform Grid')
    # plt.xlabel('Grid X Index')
    # plt.ylabel('Grid Y Index')
    # plt.savefig('transferred_field.png', dpi=150, bbox_inches='tight')
    # print(f"Saved visualization to transferred_field.png")
    # plt.close()

def main(exo_path):
    if not hasattr(omega_h, "exodus_open"):
        print("PyOmega_h was built without Exodus support "
              "(rebuild Omega_h with Omega_h_USE_SEACASExodus=ON)")
        return 1

    lib = omega_h.OmegaHLibrary()
    world = lib.world()

    exo_handle = omega_h.exodus_open(exo_path, verbose=False)
    num_steps = omega_h.exodus_get_num_time_steps(exo_handle)
    print(f"Time steps: {num_steps}")
    
    mesh = omega_h.OmegaHMesh(lib)
    mesh.set_comm(world)
    
    omega_h.read_mesh_exodus(exo_handle, mesh, verbose=True)
    omega_h.read_exodus_nodal_fields(exo_handle, mesh, 50, verbose=True) # timestep 50
    omega_h.read_exodus_element_fields(exo_handle, mesh, 50, verbose=True)
    omega_h.exodus_close(exo_handle)
    
    print_mesh_info(mesh)
    demonstrate_face_field_transfer(mesh)
    
    print(f"\nSUCCESS")
    del mesh  # Delete mesh before world to free Kokkos allocations properly
    del world
    del lib  # Delete lib last to finalize Kokkos after all objects are freed
    return 0


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "exo_path", nargs="?", default="flow.exo",
        help="Exodus file to read (e.g. pcms_testcases/exodus/flow/flow.exo)",
    )
    args = parser.parse_args()
    exit(main(args.exo_path))
