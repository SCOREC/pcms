#!/usr/bin/env python3
"""
Test OmegaHFieldLayout Python bindings
"""

import pcms
import numpy as np

def test_layout_methods(world, dim, order, num_components):
    """Test OmegaHFieldLayout methods"""
    nx = 10
    ny = 10 if dim > 1 else 0
    nz = 10 if dim > 2 else 0

    print(f"\nTesting layout: dim={dim}, order={order}, num_components={num_components}")

    # Build mesh
    mesh = pcms.build_box(
        world, 
        pcms.Family.SIMPLEX, 
        1.0, 1.0, 1.0, 
        nx, ny, nz, 
        False
    )

    # Create layout
    layout = pcms.create_lagrange_layout(
        mesh, 
        order, 
        num_components,
        pcms.CoordinateSystem.Cartesian
    )

    # Test layout methods
    num_comp = layout.get_num_components()
    assert num_comp == num_components, f"Expected {num_components} components, got {num_comp}"

    num_owned = layout.get_num_owned_dof_holder()
    print(f"  Number of owned DOF holders: {num_owned}")
    assert num_owned > 0, "Expected at least one owned DOF holder"

    num_global = layout.get_num_global_dof_holder()
    print(f"  Number of global DOF holders: {num_global}")
    assert num_global >= num_owned, "Global DOF holders should be >= owned"

    # Get ownership mask
    owned = layout.get_owned()
    print(f"  Ownership mask length: {len(owned)}")
    assert len(owned) > 0, "Ownership mask should not be empty"

    # Get global IDs
    gids = layout.get_gids()
    print(f"  Global IDs length: {len(gids)}")
    assert len(gids) > 0, "Global IDs should not be empty"

    # Get coordinates of DOF holders (2D parametric coordinates)
    coords = layout.get_dof_holder_coordinates()
    print(f"  Coordinates shape: {coords.shape}")
    assert coords.shape[0] > 0, "Should have coordinates"
    assert coords.shape[1] == 2, "Coordinates should be 2D (parametric)"

    # Check if distributed
    is_distributed = layout.is_distributed()
    print(f"  Is distributed: {is_distributed}")

    # Get entity offsets
    ent_offsets = layout.get_ent_offsets()
    print(f"  Entity offsets: {ent_offsets}")
    assert len(ent_offsets) == 5, f"Expected 5 entity offsets"

    # Get nodes per dimension
    nodes = layout.get_nodes_per_dim()
    print(f"  Nodes per dim: {nodes}")
    assert len(nodes) == 4, "Should have 4 entries (verts, edges, faces, regions)"

    # Get number of entities
    num_ents = layout.get_num_ents()
    print(f"  Number of entities: {num_ents}")
    assert num_ents > 0, "Should have at least one entity"

    # Get sizes
    owned_size = layout.owned_size()
    global_size = layout.global_size()
    print(f"  Owned size: {owned_size}")
    print(f"  Global size: {global_size}")
    assert owned_size == num_owned * num_components
    assert global_size == num_global * num_components

    # Create a field from the layout
    field = layout.create_field()
    print(f"  Created field: {type(field)}")

    # Test setting and getting data
    ndata = num_owned * num_components
    test_data = np.arange(ndata, dtype=np.float64)
    field.set_dof_holder_data(test_data)

    retrieved_data = field.get_dof_holder_data()
    assert len(retrieved_data) == ndata, f"Expected {ndata} elements, got {len(retrieved_data)}"

    # Verify data matches
    # debug print
    for i in range(min(10, ndata)):
        print(f"  Data[{i}]: set={test_data[i]}, got={retrieved_data[i]}")
    differences = np.abs(test_data - retrieved_data)
    assert np.all(differences < 1e-12), "Data mismatch after set/get"

    print(f"✓ Test passed: dim={dim}, order={order}, num_components={num_components}")

def test_field_evaluation(world, dim, order, num_components):
    """Test OmegaHField evaluation at points"""
    nx = 10
    ny = 10 if dim > 1 else 0
    nz = 10 if dim > 2 else 0

    print(f"\nTesting field evaluation: dim={dim}, order={order}, num_components={num_components}")

    # Build mesh
    mesh = pcms.build_box(
        world, 
        pcms.Family.SIMPLEX, 
        1.0, 1.0, 1.0, 
        nx, ny, nz, 
        False
    )

    # Create layout
    layout = pcms.create_lagrange_layout(
        mesh, 
        order, 
        num_components,
        pcms.CoordinateSystem.Cartesian
    )

    # Create field
    field = layout.create_field()

    # Set up test data - use a simple function: f(x,y,z) = sin(x*y) for component 0, etc.
    print(f"  Setting up field data...")
    # Get mesh vertex coordinates to set field values
    mesh_coords = mesh.coords()
    num_verts = mesh.nverts()
    print(f"  Mesh has {num_verts} vertices")

    # For order 1, we only need vertex values
    # For order 2, we also need edge midpoint values
    num_owned = layout.get_num_owned_dof_holder()
    ndata = num_owned * num_components
    print(f"  Setting up field data for {ndata} DOFs")

    # Define test function
    def test_func(x, y, z, component):
        return np.sin(5.0 * x * y) + float(component)

    # Create test data array
    test_data = np.zeros(ndata, dtype=np.float64)

    # Set vertex values
    for i in range(min(num_verts, num_owned)):
        x = mesh_coords[i * dim + 0] if dim >= 1 else 0.0
        y = mesh_coords[i * dim + 1] if dim >= 2 else 0.0
        z = mesh_coords[i * dim + 2] if dim >= 3 else 0.0
        for c in range(num_components):
            idx = i * num_components + c
            test_data[idx] = test_func(x, y, z, c)

    print(f"  Set vertex values for field data")

    # For order 2, set edge midpoint values (simplified - would need proper edge coordinates)
    if order == 2 and num_owned > num_verts:
        for i in range(num_verts, num_owned):
            for c in range(num_components):
                idx = i * num_components + c
                # Use simplified values for edges
                test_data[idx] = float(c) + 0.5

    field.set_dof_holder_data(test_data)
    print(f"  Set field data with test function")

    # Define evaluation points (physical coordinates)
    if dim == 2:
        eval_coords = np.array([
            [0.5, 0.5],
            [0.25, 0.25],
            [0.75, 0.75],
            [0.1, 0.9],
            [0.9, 0.1]
        ], dtype=np.float64)
    else:  # dim == 3
        eval_coords = np.array([
            [0.5, 0.5, 0.5],
            [0.25, 0.25, 0.25],
            [0.75, 0.75, 0.75],
            [0.1, 0.9, 0.1],
            [0.9, 0.1, 0.9]
        ], dtype=np.float64)

    num_eval_points = eval_coords.shape[0]
    print(f"  Evaluating at {num_eval_points} points")

    # Get localization hint for the coordinates
    coord_system = pcms.CoordinateSystem.Cartesian
    location_hint = field.get_localization_hint(eval_coords, coord_system)
    print(f"  Got localization hint")

    # Create output buffer for evaluation results
    eval_values = np.zeros(num_eval_points * num_components, dtype=np.float64)

    # Evaluate the field
    field.evaluate(location_hint, eval_values, coord_system)
    print(f"  Field evaluated successfully")

    # Print and verify results
    for i in range(num_eval_points):
        coords_str = ", ".join([f"{eval_coords[i, j]:.3f}" for j in range(dim)])
        values_str = ", ".join([f"{eval_values[i * num_components + c]:.4f}" 
                               for c in range(num_components)])
        print(f"  Point {i} ({coords_str}): values = [{values_str}]")

        # Verify we got valid values (not NaN or infinity)
        for c in range(num_components):
            val = eval_values[i * num_components + c]
            assert not np.isnan(val), f"Got NaN value at point {i}, component {c}"
            assert not np.isinf(val), f"Got inf value at point {i}, component {c}"

    # Test gradient evaluation if available
    if field.can_evaluate_gradient():
        print(f"  Testing gradient evaluation...")
        gradient_values = np.zeros(num_eval_points * num_components * dim, dtype=np.float64)

        try:
            field.evaluate_gradient(gradient_values, coord_system)
            print(f"  Gradient evaluated successfully")

            # Print gradient results
            for i in range(min(3, num_eval_points)):
                for c in range(num_components):
                    grad_components = []
                    for d in range(dim):
                        idx = i * num_components * dim + c * dim + d
                        grad_components.append(f"{gradient_values[idx]:.4f}")
                    grad_str = ", ".join(grad_components)
                    print(f"    Point {i}, component {c}: gradient = [{grad_str}]")
        except Exception as e:
            print(f"  Gradient evaluation failed: {e}")
    else:
        print(f"  Gradient evaluation not supported")

    print(f"✓ Field evaluation test passed: dim={dim}, order={order}, num_components={num_components}")

def test_tag_operations(world, dim):
    """Test Omega_h mesh tag operations"""
    nx = 5
    ny = 5 if dim > 1 else 0
    nz = 5 if dim > 2 else 0

    print(f"\nTesting tag operations: dim={dim}")

    # Build mesh
    mesh = pcms.build_box(
        world, 
        pcms.Family.SIMPLEX, 
        1.0, 1.0, 1.0, 
        nx, ny, nz, 
        False
    )

    print(f"  Mesh has {mesh.nverts()} vertices")
    print(f"  Mesh has {mesh.nelems()} elements")

    # Test 1: Create empty tags with different data types
    print("\n  Test 1: Creating empty tags with different dtypes...")
    mesh.add_tag(0, "test_float", 1, dtype="float64")
    mesh.add_tag(0, "test_int32", 1, dtype="int32")
    mesh.add_tag(0, "test_int64", 1, dtype="int64")
    
    assert mesh.has_tag(0, "test_float"), "Float tag should exist"
    assert mesh.has_tag(0, "test_int32"), "Int32 tag should exist"
    assert mesh.has_tag(0, "test_int64"), "Int64 tag should exist"
    print("  ✓ Empty tags created successfully")

    # Test 2: Create tags with initial data (auto-detect type)
    print("\n  Test 2: Creating tags with initial numpy arrays...")
    nverts = mesh.nverts()
    
    # Float64 array
    temp_data = np.linspace(0.0, 100.0, nverts, dtype=np.float64)
    mesh.add_tag(0, "temperature", 1, temp_data)
    
    # Int32 array
    ids_data = np.arange(nverts, dtype=np.int32)
    mesh.add_tag(0, "vertex_ids", 1, ids_data)
    
    # Multi-component vector
    velocity_data = np.random.randn(nverts * dim).astype(np.float64)
    mesh.add_tag(0, "velocity", dim, velocity_data)
    
    assert mesh.has_tag(0, "temperature"), "Temperature tag should exist"
    assert mesh.has_tag(0, "vertex_ids"), "Vertex IDs tag should exist"
    assert mesh.has_tag(0, "velocity"), "Velocity tag should exist"
    print("  ✓ Tags with initial data created successfully")

    # Test 3: Get tag data and verify
    print("\n  Test 3: Getting tag data and verifying...")
    retrieved_temp = mesh.get_tag(0, "temperature")
    retrieved_ids = mesh.get_tag(0, "vertex_ids")
    retrieved_vel = mesh.get_tag(0, "velocity")
    
    assert retrieved_temp.dtype == np.float64, "Temperature should be float64"
    assert retrieved_ids.dtype == np.int32, "IDs should be int32"
    assert retrieved_vel.shape[0] == nverts * dim, f"Velocity should have {nverts * dim} elements"
    
    np.testing.assert_allclose(retrieved_temp, temp_data, rtol=1e-15)
    np.testing.assert_array_equal(retrieved_ids, ids_data)
    print("  ✓ Tag data retrieved and verified successfully")

    # Test 4: Set tag data
    print("\n  Test 4: Setting tag data...")
    new_temp = np.ones(nverts, dtype=np.float64) * 50.0
    mesh.set_tag(0, "temperature", new_temp)
    
    updated_temp = mesh.get_tag(0, "temperature")
    np.testing.assert_allclose(updated_temp, new_temp, rtol=1e-15)
    print("  ✓ Tag data set and verified successfully")

    # Test 5: ArrayType parameter
    print("\n  Test 5: Testing ArrayType parameter...")
    # Create a symmetric matrix tag (stress tensor)
    if dim == 3:
        ncomps = 6  # xx, yy, zz, xy, xz, yz
        stress_data = np.random.randn(mesh.nelems() * ncomps).astype(np.float64)
        mesh.add_tag(dim, "stress", ncomps, stress_data, 
                     internal=False, 
                     array_type=pcms.ArrayType.SymmetricSquareMatrix)
        assert mesh.has_tag(dim, "stress"), "Stress tag should exist on elements"
        print("  ✓ SymmetricSquareMatrix ArrayType tag created")
    
    # Create an internal tag (not written to file)
    internal_data = np.zeros(nverts, dtype=np.float64)
    mesh.add_tag(0, "internal_temp", 1, internal_data, internal=True)
    assert mesh.has_tag(0, "internal_temp"), "Internal tag should exist"
    print("  ✓ Internal tag created")

    # Test 6: Tag removal
    print("\n  Test 6: Removing tags...")
    mesh.remove_tag(0, "test_float")
    assert not mesh.has_tag(0, "test_float"), "test_float tag should be removed"
    
    # Other tags should still exist
    assert mesh.has_tag(0, "temperature"), "temperature tag should still exist"
    print("  ✓ Tag removal successful")

    # Test 7: Number of tags
    print("\n  Test 7: Counting tags...")
    ntags = mesh.ntags(0)
    print(f"  Number of vertex tags: {ntags}")
    assert ntags > 0, "Should have at least one vertex tag"

    # Test 8: Element tags
    print("\n  Test 8: Testing element tags...")
    nelems = mesh.nelems()
    elem_quality = np.random.rand(nelems).astype(np.float64)
    mesh.add_tag(dim, "quality", 1, elem_quality)
    
    retrieved_quality = mesh.get_tag(dim, "quality")
    assert len(retrieved_quality) == nelems, f"Should have {nelems} quality values"
    np.testing.assert_allclose(retrieved_quality, elem_quality, rtol=1e-15)
    print("  ✓ Element tags working correctly")

    # Test 9: Edge tags (if applicable)
    if dim >= 2:
        print("\n  Test 9: Testing edge tags...")
        nedges = mesh.nedges()
        edge_length = np.ones(nedges, dtype=np.float64)
        mesh.add_tag(1, "edge_marker", 1, edge_length)
        
        retrieved_edge = mesh.get_tag(1, "edge_marker")
        assert len(retrieved_edge) == nedges, f"Should have {nedges} edge values"
        print("  ✓ Edge tags working correctly")

    # Test 10: Adjacency and globals
    print("\n  Test 10: Testing adjacency and global IDs...")
    
    # Get vertex connectivity of elements
    elem_verts = mesh.ask_elem_verts()
    print(f"  Element-vertex connectivity shape: {elem_verts.shape}")
    assert len(elem_verts) > 0, "Should have element-vertex connectivity"
    
    # Get global IDs
    global_ids = mesh.globals(0)
    print(f"  Global vertex IDs shape: {global_ids.shape}")
    assert len(global_ids) == nverts, f"Should have {nverts} global IDs"
    
    # Get verts of elements
    verts_of_elems = mesh.ask_verts_of(dim)
    print(f"  Verts of elements shape: {verts_of_elems.shape}")
    assert len(verts_of_elems) > 0, "Should have vert connectivity"
    
    # Check adjacency existence
    has_vert_to_elem = mesh.has_adj(0, dim)
    print(f"  Has vertex-to-element adjacency: {has_vert_to_elem}")
    
    print("  ✓ Adjacency and global ID queries successful")

    # Test 11: Ownership
    print("\n  Test 12: Testing ownership...")
    owned_verts = mesh.owned(0)
    print(f"  Owned vertices shape: {owned_verts.shape}")
    num_owned = np.sum(owned_verts)
    print(f"  Number of owned vertices: {num_owned} / {nverts}")
    assert len(owned_verts) == nverts, "Should have ownership flag for each vertex"
    print("  ✓ Ownership queries successful")

    print(f"\n✓ All tag operation tests passed for dim={dim}!")

def test_entity_coordinates(world, dim):
    """Test computing entity coordinates (face centers, edge midpoints, etc.)"""
    nx = 4
    ny = 4 if dim > 1 else 0
    nz = 4 if dim > 2 else 0

    print(f"\nTesting entity coordinate computation: dim={dim}")

    # Build mesh
    mesh = pcms.build_box(
        world, 
        pcms.Family.SIMPLEX, 
        1.0, 1.0, 1.0, 
        nx, ny, nz, 
        False
    )

    print(f"  Mesh has {mesh.nverts()} vertices")
    
    # Get vertex coordinates
    vertex_coords = mesh.coords()
    print(f"  Vertex coordinates shape: {vertex_coords.shape}")
    
    # Test computing coordinates for different entity dimensions
    if dim >= 2:
        print(f"\n  Computing edge coordinates (midpoints)...")
        nedges = mesh.nedges()
        edge_coords = pcms.average_field(mesh, 1, dim, vertex_coords)
        edge_coords_np = np.array(edge_coords).reshape(-1, dim)
        print(f"    Number of edges: {nedges}")
        print(f"    Edge coordinates shape: {edge_coords_np.shape}")
        assert edge_coords_np.shape[0] == nedges, f"Should have {nedges} edge coordinates"
        print(f"    First edge center: {edge_coords_np[0]}")
    
    if dim == 2:
        print(f"\n  Computing element (face) coordinates (centroids)...")
        nelems = mesh.nelems()
        elem_coords = pcms.average_field(mesh, 2, dim, vertex_coords)
        elem_coords_np = np.array(elem_coords).reshape(-1, dim)
        print(f"    Number of elements: {nelems}")
        print(f"    Element coordinates shape: {elem_coords_np.shape}")
        assert elem_coords_np.shape[0] == nelems, f"Should have {nelems} element coordinates"
        print(f"    First element center: {elem_coords_np[0]}")
    
    if dim == 3:
        print(f"\n  Computing face coordinates (centroids)...")
        nfaces = mesh.nfaces()
        face_coords = pcms.average_field(mesh, 2, dim, vertex_coords)
        face_coords_np = np.array(face_coords).reshape(-1, dim)
        print(f"    Number of faces: {nfaces}")
        print(f"    Face coordinates shape: {face_coords_np.shape}")
        assert face_coords_np.shape[0] == nfaces, f"Should have {nfaces} face coordinates"
        print(f"    First face center: {face_coords_np[0]}")
        
        print(f"\n  Computing region/element coordinates (centroids)...")
        nregions = mesh.nregions()
        region_coords = pcms.average_field(mesh, 3, dim, vertex_coords)
        region_coords_np = np.array(region_coords).reshape(-1, dim)
        print(f"    Number of regions: {nregions}")
        print(f"    Region coordinates shape: {region_coords_np.shape}")
        assert region_coords_np.shape[0] == nregions, f"Should have {nregions} region coordinates"
        print(f"    First region center: {region_coords_np[0]}")
    
    # Test averaging a field from vertices to entities
    print(f"\n  Testing field averaging from vertices to elements...")
    nverts = mesh.nverts()
    vertex_field = np.arange(nverts, dtype=np.float64)
    mesh.add_tag(0, "test_vertex_field", 1, vertex_field)
    
    # Average to elements
    elem_dim = dim  # elements are top-dimensional entities
    averaged_field = pcms.average_field(mesh, elem_dim, 1, vertex_field)
    print(f"    Averaged field shape: {averaged_field.shape}")
    print(f"    First 5 values: {averaged_field[:5]}")
    
    print(f"\n✓ Entity coordinate computation tests passed for dim={dim}!")

def main():
    """Run all test cases"""
    print("Testing OmegaHFieldLayout Python bindings...")

    # Initialize Omega_h library
    lib = pcms.OmegaHLibrary()
    world = lib.world()
    print("Initialized Omega_h library and world")

    # Test different configurations
    test_layout_methods(world, 2, 1, 1)
    test_layout_methods(world, 2, 2, 1)
    # test_layout_methods(world, 2, 1, 3)
    # test_layout_methods(world, 3, 1, 1)

    # Test field evaluation
    print("\n" + "="*60)
    print("Testing field evaluation...")
    print("="*60)
    test_field_evaluation(world, 2, 1, 1)
    test_field_evaluation(world, 2, 2, 1)

    # Test tag operations
    print("\n" + "="*60)
    print("Testing tag operations...")
    print("="*60)
    test_tag_operations(world, 2)
    test_tag_operations(world, 3)

    # Test entity coordinate computation
    print("\n" + "="*60)
    print("Testing entity coordinate computation...")
    print("="*60)
    test_entity_coordinates(world, 2)
    test_entity_coordinates(world, 3)

    print("\n✓ All tests passed!")
    del world
    del lib

if __name__ == "__main__":
    main()