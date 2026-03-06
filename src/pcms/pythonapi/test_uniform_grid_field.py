"""
Example test file for uniform grid field Python bindings.

This demonstrates how to use the UniformGridField and UniformGridFieldLayout
classes in Python, similar to how OmegaHField is used.
"""
import numpy as np
import pcms


def test_uniform_grid_field_creation():
    """Test creating a 2D uniform grid field."""
    # Create a 2D uniform grid
    grid = pcms.UniformGrid2D()
    grid.bot_left = [0.0, 0.0]
    grid.edge_length = [10.0, 10.0]
    grid.divisions = [4, 4]  # 4x4 cells
    
    print(f"Grid cells: {grid.get_num_cells()}")  # Should be 16
    
    # Create field layout with 1 component (scalar field)
    layout = pcms.UniformGridFieldLayout2D(
        grid, 1, pcms.CoordinateSystem.Cartesian
    )
    
    print(f"Num components: {layout.get_num_components()}")
    print(f"Num owned DOF holders: {layout.get_num_owned_dof_holder()}")
    print(f"Num vertices: {layout.get_num_vertices()}")  # 5x5 = 25 vertices
    
    # Create field
    field = layout.create_field()
    print(f"Field created: {field is not None}")
    
    return grid, layout, field


def test_uniform_grid_field_data_operations():
    """Test setting and getting field data."""
    # Create a simple grid
    grid = pcms.UniformGrid2D()
    grid.bot_left = [0.0, 0.0]
    grid.edge_length = [10.0, 10.0]
    grid.divisions = [2, 2]  # 2x2 cells, 3x3 vertices
    
    layout = pcms.UniformGridFieldLayout2D(
        grid, 1, pcms.CoordinateSystem.Cartesian
    )
    ug_field = layout.create_field()
    
    # Get the number of DOF holders (vertices)
    num_vertices = layout.get_num_vertices()
    print(f"Number of vertices: {num_vertices}")  # Should be 9
    
    # Create data for vertices (initialize with values)
    data = np.arange(num_vertices, dtype=np.float64)
    print(f"Setting data: {data}")
    
    # Set field data
    ug_field.set_dof_holder_data(data)
    
    # Get field data back
    retrieved_data = ug_field.get_dof_holder_data()
    print(f"Retrieved data: {retrieved_data}")
    
    # Verify they match
    assert np.allclose(data, retrieved_data)
    print("Data successfully set and retrieved!")


def test_uniform_grid_field_mdspan_2d():
    """Test 2D mdspan output as numpy array."""
    grid = pcms.UniformGrid2D()
    grid.bot_left = [0.0, 0.0]
    grid.edge_length = [10.0, 10.0]
    grid.divisions = [2, 3]  # 3x4 vertices

    layout = pcms.UniformGridFieldLayout2D(
        grid, 1, pcms.CoordinateSystem.Cartesian
    )
    ug_field = layout.create_field()

    num_vertices = layout.get_num_vertices()
    data = np.arange(num_vertices, dtype=np.float64)
    ug_field.set_dof_holder_data(data)

    mdspan = ug_field.to_numpy()
    expected = data.reshape((grid.divisions[0] + 1, grid.divisions[1] + 1))
    assert type(mdspan) == np.ndarray

    assert mdspan.shape == expected.shape
    assert np.allclose(mdspan, expected)
    print("2D mdspan output verified!")


def test_uniform_grid_field_evaluation():
    """Test evaluating field at specific coordinates."""
    # Create grid
    grid = pcms.UniformGrid2D()
    grid.bot_left = [0.0, 0.0]
    grid.edge_length = [10.0, 10.0]
    grid.divisions = [2, 2]
    
    layout = pcms.UniformGridFieldLayout2D(
        grid, 1, pcms.CoordinateSystem.Cartesian
    )
    ug_field = layout.create_field()
    
    # Set simple linear field data: f(x,y) = x + y
    num_vertices = layout.get_num_vertices()
    coords = layout.get_dof_holder_coordinates()
    data = np.zeros(num_vertices)
    for i in range(num_vertices):
        data[i] = coords[i, 0] + coords[i, 1]
    
    ug_field.set_dof_holder_data(data)
    
    # Evaluate at some points
    eval_coords = np.array([[2.5, 2.5], [7.5, 7.5]])
    
    # Get localization hint
    hint = ug_field.get_localization_hint(
        eval_coords, pcms.CoordinateSystem.Cartesian
    )
    
    # Evaluate field
    results = np.zeros(len(eval_coords))
    ug_field.evaluate(
        hint, results, pcms.CoordinateSystem.Cartesian
    )
    
    print(f"Evaluation results: {results}")
    print(f"Expected (approximately): [5.0, 15.0]")


def test_uniform_grid_closest_cell():
    """Test finding closest cell to a point."""
    grid = pcms.UniformGrid2D()
    grid.bot_left = [0.0, 0.0]
    grid.edge_length = [10.0, 10.0]
    grid.divisions = [4, 4]
    
    # Test point in middle of grid
    point = np.array([5.0, 5.0])
    cell_id = grid.closest_cell_id(point)
    print(f"Point {point} is in cell {cell_id}")
    
    # Test point outside grid (should clamp to closest cell)
    point_outside = np.array([-1.0, -1.0])
    cell_id_outside = grid.closest_cell_id(point_outside)
    print(f"Point {point_outside} (outside) maps to cell {cell_id_outside}")


def test_3d_uniform_grid():
    """Test 3D uniform grid field."""
    grid = pcms.UniformGrid3D()
    grid.bot_left = [0.0, 0.0, 0.0]
    grid.edge_length = [10.0, 10.0, 10.0]
    grid.divisions = [2, 2, 2]  # 2x2x2 cells
    
    print(f"3D Grid cells: {grid.get_num_cells()}")  # Should be 8
    
    layout = pcms.UniformGridFieldLayout3D(
        grid, 1, pcms.CoordinateSystem.Cartesian
    )
    
    print(f"3D Grid vertices: {layout.get_num_vertices()}")  # 3x3x3 = 27
    
    ug_field = layout.create_field()
    
    # Set and get data
    num_vertices = layout.get_num_vertices()
    data = np.ones(num_vertices) * 42.0
    ug_field.set_dof_holder_data(data)
    
    retrieved = ug_field.get_dof_holder_data()
    assert np.allclose(data, retrieved)
    print("3D field data successfully set and retrieved!")


def test_uniform_grid_field_mdspan_3d():
    """Test 3D mdspan output as numpy array."""
    grid = pcms.UniformGrid3D()
    grid.bot_left = [0.0, 0.0, 0.0]
    grid.edge_length = [10.0, 10.0, 10.0]
    grid.divisions = [2, 1, 3]  # 3x2x4 vertices

    layout = pcms.UniformGridFieldLayout3D(
        grid, 1, pcms.CoordinateSystem.Cartesian
    )
    ug_field = layout.create_field()

    num_vertices = layout.get_num_vertices()
    data = np.arange(num_vertices, dtype=np.float64)
    ug_field.set_dof_holder_data(data)

    mdspan = ug_field.to_numpy()
    expected = data.reshape(
        (grid.divisions[0] + 1, grid.divisions[1] + 1, grid.divisions[2] + 1)
    )

    assert mdspan.shape == expected.shape
    assert np.allclose(mdspan, expected)
    print("3D mdspan output verified!")


def test_uniform_grid_workflow(world):
    """
    Test complete UniformGrid workflow with Omega_h mesh and field interpolation.
    
    This test:
    1. Creates an Omega_h mesh
    2. Creates a uniform grid from the mesh
    3. Creates a binary mask field
    4. Creates and initializes an Omega_h field with f(x,y) = x + 2*y
    5. Transfers the field to uniform grid via interpolation
    6. Verifies the transferred values
    """
    # Create a simple 2D box mesh: 1.0 x 1.0 domain with 4x4 elements
    mesh = pcms.build_box(
        world,
        pcms.Family.SIMPLEX,
        1.0, 1.0, 0.0,
        4, 4, 0,
        False
    )
    
    # Create uniform grid from mesh with 4x4 divisions
    grid = pcms.create_uniform_grid_from_mesh(mesh, [4, 4])
    print(f"Created uniform grid with {grid.get_num_cells()} cells")
    
    # Create binary mask field (returns tuple of (layout, field))
    mask_layout, mask_field = pcms.create_uniform_grid_binary_field(mesh, [4, 4])
    mask_data = mask_field.get_dof_holder_data()
    print(f"Created mask field with {len(mask_data)} vertices")
    
    # Create Omega_h field layout with linear elements
    omega_h_layout = pcms.create_lagrange_layout(
        mesh, 1
    )
    omega_h_field = omega_h_layout.create_field()
    
    # Initialize omega_h field with f(x,y) = x + 2*y
    coords = omega_h_layout.get_dof_holder_coordinates()
    num_nodes = omega_h_layout.get_num_owned_dof_holder()
    omega_h_data = np.zeros(num_nodes)
    
    for i in range(num_nodes):
        x = coords[i, 0]
        y = coords[i, 1]
        omega_h_data[i] = x + 2.0 * y
    
    omega_h_field.set_dof_holder_data(omega_h_data)
    print(f"Initialized Omega_h field with {num_nodes} nodes")
    
    # Create uniform grid field layout
    ug_layout = pcms.UniformGridFieldLayout2D(
        grid, 1, pcms.CoordinateSystem.Cartesian
    )
    ug_field = ug_layout.create_field()
    
    # Transfer from omega_h field to uniform grid field using interpolation
    pcms.interpolate_field(omega_h_field, ug_field)
    print("Field interpolation completed")
    
    # Get uniform grid field data and coordinates
    ug_field_data = ug_field.get_dof_holder_data()
    ug_coords = ug_layout.get_dof_holder_coordinates()
    
    # Verify ug_field values
    print("\nVerifying uniform grid field values:")
    num_errors = 0
    for j in range(grid.divisions[1] + 1):
        for i in range(grid.divisions[0] + 1):
            vertex_id = j * (grid.divisions[0] + 1) + i
            x = ug_coords[vertex_id, 0]
            y = ug_coords[vertex_id, 1]
            expected = x + 2.0 * y
            actual = ug_field_data[vertex_id]
            
            error = abs(expected - actual)
            if error > 1e-10:
                num_errors += 1
                print(f"  Vertex ({i}, {j}) at ({x:.2f}, {y:.2f}): "
                      f"expected {expected:.4f}, got {actual:.4f}, error {error:.2e}")
            
    if num_errors == 0:
        print("All uniform grid field values verified successfully!")
    else:
        print(f"Found {num_errors} verification errors")
        
    # Verify mask field values (all should be 1 for vertices inside the mesh)
    print("\nVerifying mask field values:")
    mask_errors = 0
    nx, ny = grid.divisions[0] + 1, grid.divisions[1] + 1
    print(f"\nBinary mask field ({nx}x{ny} vertices) visualization:")
    print("(Each position shows the mask value at that vertex)")
    for j in range(ny - 1, -1, -1):  # Print from top to bottom
        row = []
        for i in range(nx):
            vertex_id = j * nx + i
            mask_value = mask_data[vertex_id]
            row.append(str(int(mask_value)))
            if mask_value != 1.0:
                mask_errors += 1
        print(f"  Row {j}: [{' '.join(row)}]")
    
    if mask_errors == 0:
        print("All mask field values verified successfully!")
    else:
        print(f"Found {mask_errors} mask errors")
    
    # Serialize the uniform grid field
    print("\nSerializing uniform grid field:")
    num_vertices = ug_layout.get_num_vertices()
    buffer = np.zeros(num_vertices)
    permutation = np.arange(num_vertices, dtype=np.int32)  # Identity permutation
    
    bytes_written = ug_field.serialize(buffer, permutation)
    print(f"Serialized {bytes_written} values to buffer")
    
    print(f"\nSerialized uniform grid field data ({grid.divisions[0]+1}x{grid.divisions[1]+1} vertices) visualization:")
    print("(Each vertex shows its field value)")
    for j in range(grid.divisions[1], -1, -1):  # Print from top to bottom
        row_values = []
        for i in range(grid.divisions[0] + 1):
            vertex_id = j * (grid.divisions[0] + 1) + i
            value = buffer[vertex_id]
            row_values.append(f"{value:6.3f}")
        print(f"  Row {j}: [{' '.join(row_values)}]")
    
    # Also print with coordinates for reference
    print(f"\nDetailed vertex information:")
    for j in range(grid.divisions[1] + 1):
        for i in range(grid.divisions[0] + 1):
            vertex_id = j * (grid.divisions[0] + 1) + i
            x = ug_coords[vertex_id, 0]
            y = ug_coords[vertex_id, 1]
            value = buffer[vertex_id]
            print(f"  V[{i},{j}] (id={vertex_id:2d}) at ({x:.3f}, {y:.3f}): value = {value:.6f}")
    
    # Test deserialization
    print("\nTesting deserialization:")
    ug_field_copy = ug_layout.create_field()
    ug_field_copy.deserialize(buffer, permutation)
    
    ug_field_copy_data = ug_field_copy.get_dof_holder_data()
    if np.allclose(ug_field_data, ug_field_copy_data):
        print("✓ Deserialization successful - data matches!")
    else:
        print("✗ Deserialization failed - data mismatch!")
        print(f"  Max difference: {np.max(np.abs(ug_field_data - ug_field_copy_data))}")
        
    assert num_errors == 0, f"Field verification failed with {num_errors} errors"
    assert mask_errors == 0, f"Mask verification failed with {mask_errors} errors"


def test_omega_h_to_omega_h_transfer_workflow(world):
    """
    Test field transfer from one Omega_h mesh to another Omega_h mesh.

    This test:
    1. Creates a source Omega_h mesh
    2. Creates a target Omega_h mesh
    3. Creates Lagrange layouts and fields on both
    4. Initializes the source with f(x,y) = x + 2*y
    5. Transfers the field to the target mesh
    6. Verifies the transferred values at target DOF holders
    """
    # Source mesh (coarser)
    src_mesh = pcms.build_box(
        world,
        pcms.Family.SIMPLEX,
        1.0, 1.0, 0.0,
        4, 4, 0,
        False
    )

    # Target mesh (finer)
    tgt_mesh = pcms.build_box(
        world,
        pcms.Family.SIMPLEX,
        1.0, 1.0, 0.0,
        8, 8, 0,
        False
    )

    # Create Lagrange layouts and fields
    src_layout = pcms.create_lagrange_layout(
        src_mesh, 1
    )
    tgt_layout = pcms.create_lagrange_layout(
        tgt_mesh, 1
    )
    src_field = src_layout.create_field()
    tgt_field = tgt_layout.create_field()

    # Initialize source field with f(x,y) = x + 2*y
    src_coords = src_layout.get_dof_holder_coordinates()
    src_num_nodes = src_layout.get_num_owned_dof_holder()
    src_data = np.zeros(src_num_nodes)
    for i in range(src_num_nodes):
        x = src_coords[i, 0]
        y = src_coords[i, 1]
        src_data[i] = x + 2.0 * y
    src_field.set_dof_holder_data(src_data)

    # Transfer field to target mesh
    pcms.interpolate_field(src_field, tgt_field)

    # Verify target field values at target DOF holders
    tgt_coords = tgt_layout.get_dof_holder_coordinates()
    tgt_data = tgt_field.get_dof_holder_data()
    errors = 0
    for i in range(tgt_layout.get_num_owned_dof_holder()):
        x = tgt_coords[i, 0]
        y = tgt_coords[i, 1]
        expected = x + 2.0 * y
        actual = tgt_data[i]
        if abs(expected - actual) > 1e-10:
            errors += 1
    assert errors == 0, f"Omega_h transfer verification failed with {errors} errors"


if __name__ == "__main__":
    lib = pcms.OmegaHLibrary()
    world = lib.world()
    print("=" * 60)
    print("Testing UniformGrid Field Creation")
    print("=" * 60)
    test_uniform_grid_field_creation()
    
    print("\n" + "=" * 60)
    print("Testing UniformGrid Field Data Operations")
    print("=" * 60)
    test_uniform_grid_field_data_operations()

    print("\n" + "=" * 60)
    print("Testing UniformGrid Field mdspan (2D)")
    print("=" * 60)
    test_uniform_grid_field_mdspan_2d()
    
    print("\n" + "=" * 60)
    print("Testing UniformGrid Field Evaluation")
    print("=" * 60)
    test_uniform_grid_field_evaluation()
    
    print("\n" + "=" * 60)
    print("Testing Closest Cell ID")
    print("=" * 60)
    test_uniform_grid_closest_cell()
    
    print("\n" + "=" * 60)
    print("Testing 3D UniformGrid")
    print("=" * 60)
    test_3d_uniform_grid()

    print("\n" + "=" * 60)
    print("Testing UniformGrid Field mdspan (3D)")
    print("=" * 60)
    test_uniform_grid_field_mdspan_3d()
    
    print("\n" + "=" * 60)
    print("Testing UniformGrid Workflow")
    print("=" * 60)
    test_uniform_grid_workflow(world)

    print("\n" + "=" * 60)
    print("Testing Omega_h to Omega_h Transfer Workflow")
    print("=" * 60)
    test_omega_h_to_omega_h_transfer_workflow(world)
    
    print("\n" + "=" * 60)
    print("All tests passed!")
    print("=" * 60)

    del world
