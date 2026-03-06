import pcms
import numpy as np

def test_copy(world, dim, order, num_components):
    """Test copying omega_h_field2 data"""
    nx = 100
    ny = 100 if dim > 1 else 0
    nz = 100 if dim > 2 else 0
    print(f"\nStarting test: dim={dim}, order={order}, num_components={num_components}")

    # Build mesh
    mesh = pcms.build_box(
        world, 
        pcms.Family.SIMPLEX, 
        1.0, 1.0, 1.0, 
        nx, ny, nz, 
        False
    )
    print(f"  Built {dim}D mesh with {nx}x{ny}x{nz} elements")
    print(f"  Mesh type: {type(mesh)}")
    print(f"  Mesh object: {mesh}")

    # Create layout
    print("  About to create layout...")
    layout = pcms.create_lagrange_layout(
        mesh, 
        order, 
        num_components,
        pcms.CoordinateSystem.Cartesian
    )
    print(f"  Layout created successfully")
    print(f"Testing dim={dim}, order={order}, num_components={num_components}...")

    # Get number of data points
    ndata = layout.get_num_owned_dof_holder() * num_components
    print(f"  Number of data points: {ndata}")

    # Create sequential array of IDs
    ids = np.arange(ndata, dtype=np.float64)
    print(f"  Created array of IDs from 0 to {ndata-1}")

    # Create original field and set data
    original = layout.create_field()
    print("  Created original field")
    original.set_dof_holder_data(ids)
    print("  Set data in original field")

    # Create copied field and copy data
    copied = layout.create_field()
    pcms.copy_field(original, copied)
    print("  Copied data to new field")

    # Get copied data
    copied_array = copied.get_dof_holder_data()

    # Verify the copy
    assert len(copied_array) == ndata, f"Expected {ndata} elements, got {len(copied_array)}"

    # Check that all values match
    differences = np.abs(ids - copied_array)
    num_matches = np.sum(differences < 1e-12)

    assert num_matches == ndata, f"Only {num_matches}/{ndata} elements matched"

    print(f"✓ Test passed: dim={dim}, order={order}, num_components={num_components}")

def main():
    """Run all test cases"""
    print("Testing copy omega_h_field2 data...")

    # Initialize Omega_h library
    lib = pcms.OmegaHLibrary()
    world = lib.world()
    print("Initialized Omega_h library and world")

    # Run test cases
    test_copy(world, 2, 1, 1)
    print("first passed")
    test_copy(world, 2, 2, 1)

    print("\nAll tests passed!")

    del world
    del lib

if __name__ == "__main__":
    main()