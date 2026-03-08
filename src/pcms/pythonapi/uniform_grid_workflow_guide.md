# Unstructured mesh field transfer to structured mesh field workflow guide

This guide demonstrates how to transfer unstructured mesh fields to uniform grid fields in PCMS Python API, including mesh creation, field interpolation, and data serialization.

## Overview

The workflow allows you to:
1. Create a structured uniform grid as the bounding box as an unstructured Omega_h mesh
2. Create binary mask fields indicating which cells are inside the mesh
3. Interpolate field data from Omega_h mesh fields to uniform grid fields
4. Serialize and deserialize field data for I/O operations

## Complete Workflow Example

### Step 1: Initialize Omega_h Library

Initialize the Omega_h library which manages parallel communication.

```python
import pcms
import numpy as np

# Create library and world communicator
lib = pcms.OmegaHLibrary()
world = lib.world()
```


---

### Step 2: Create an Omega_h Mesh

You can either create a mesh or read an existing mesh from a file.

#### Option A: Create Mesh Programmatically

```python
# Create a 2D box mesh: 1.0 x 1.0 domain with 4x4 elements
mesh = pcms.build_box(
    world,                          # Communicator
    pcms.Family.SIMPLEX,        # Element type
    1.0, 1.0, 0.0,                 # Domain size: x, y, z
    4, 4, 0,                       # Number of elements: nx, ny, nz
    False                          # Symmetric flag
)
```

#### Option B: Read Mesh from File

```python
# Auto-detect file format and read
mesh = pcms.read_mesh_file("my_mesh.osh", world)

# Or use format-specific readers for explicit control:

# Binary format (.osh)
mesh = pcms.read_mesh_binary("mesh.osh", lib)

# Gmsh format (.msh)
mesh = pcms.read_mesh_gmsh("mesh.msh", world)

# VTK format (.pvtu for parallel files)
mesh = pcms.read_mesh_parallel_vtk("mesh.pvtu", world)

# Exodus format (.exo) - if SEACAS support is enabled
# Using file handle API
exo_handle = pcms.exodus_open("mesh.exo", verbose=False)
num_steps = pcms.exodus_get_num_time_steps(exo_handle)
mesh = pcms.OmegaHMesh(lib)
mesh.set_comm(world)
pcms.read_mesh_exodus(exo_handle, mesh, verbose=False)
pcms.exodus_close(exo_handle)

# ADIOS2 format (.bp) - if ADIOS2 support is enabled
mesh = pcms.read_mesh_adios2("mesh.bp", lib, prefix="")

# MESHB format (.mesh) - if LIBMESHB support is enabled
mesh = pcms.OmegaHMesh(lib)
mesh.set_comm(world)
pcms.read_mesh_meshb(mesh, "mesh.mesh")
# Read solution/resolution data if available
pcms.read_meshb_sol(mesh, "mesh.sol", sol_name="resolution")
```

**PS**: File formats such as Exodus or libMeshB require that PCMS be built with the respective libraries enabled.

Checking Mesh Properties:

```python
print(f"Mesh dimension: {mesh.dim()}")
print(f"Number of vertices: {mesh.nverts()}")
print(f"Number of elements: {mesh.nelems()}")
print(f"Mesh family: {mesh.family()}")
```

---

### Step 3: Create Uniform Grid from Mesh

Generate a structured uniform grid that covers the bounding box of the mesh.

```python
# Create uniform grid with 4x4 cell divisions
grid = pcms.create_uniform_grid_from_mesh(mesh, [4, 4])
print(f"Created uniform grid with {grid.get_num_cells()} cells")
```

---

### Step 4: Create Binary Mask Field

Create a mask where each uniform grid vertex is marked as 1 (inside mesh) or 0 (outside mesh).

```python
# Create binary mask indicating which vertices are inside the mesh
# Returns tuple of (layout, field) - layout must be kept alive while using field
mask_layout, mask_field = pcms.create_uniform_grid_binary_field(mesh, [4, 4])
print(f"Created mask field with {mask_layout.get_num_vertices()} vertices")
```

**Accessing Values**:
```python
# Access mask value for vertex (i, j)
vertex_id = j * (grid.divisions[0] + 1) + i
mask_data = mask_field.get_dof_holder_data()
mask_value = mask_data[vertex_id]  # Returns 0.0 or 1.0
```

---

### Step 5: Create and Initialize Omega_h Field

#### Option A: Create Field Programmatically

```python
# Create field layout with linear (order=1) elements and 1 component (scalar)
omega_h_layout = pcms.create_lagrange_layout(
    mesh, 
    1,
    1,
    pcms.CoordinateSystem.Cartesian
)
omega_h_field = omega_h_layout.create_field()

# Initialize field with f(x,y) = x + 2*y
coords = omega_h_layout.get_dof_holder_coordinates()
num_nodes = omega_h_layout.get_num_owned_dof_holder()
omega_h_data = np.zeros(num_nodes)

for i in range(num_nodes):
    x = coords[i, 0]
    y = coords[i, 1]
    omega_h_data[i] = x + 2.0 * y

omega_h_field.set_dof_holder_data(omega_h_data)

# Set out-of-bounds behavior (FILL with 0.0 for points outside mesh)
omega_h_field.set_out_of_bounds_mode(pcms.OutOfBoundsMode.FILL, 0.0)
```

#### Option B: Use Existing Element/Face Field from Mesh Tags

Fields data are saved as tags on mesh entities (e.g., vertices, edges, faces, elements). You can retrieve these fields and convert them to vertex-based fields for interpolation.

```python
# Get field from mesh tag (e.g., element-centered field).
face_dim = 2
tag_index = 0  # Index of the tag to use
face_tag = mesh.get_tag(face_dim, tag_index)  # Get tag by index
face_field = mesh.get_tag(face_dim, face_tag.name())

# Convert element field to vertex field using averaging
vertex_field = pcms.map_entity_field_to_vertices_average(mesh, face_field, face_dim)

# Create vertex-based field layout and set data
omega_h_layout = pcms.create_lagrange_layout(mesh, 1, 1, pcms.CoordinateSystem.Cartesian)
omega_h_field = omega_h_layout.create_field()
omega_h_field.set_dof_holder_data(vertex_field)

# Set out-of-bounds behavior (FILL with 0.0 for points outside mesh)
omega_h_field.set_out_of_bounds_mode(pcms.OutOfBoundsMode.FILL, 0.0)
```

---

### Step 6: Create Uniform Grid Field Layout

Set up the data structure for storing field values on the uniform grid vertices.

```python
# Create uniform grid field layout
ug_layout = pcms.UniformGridFieldLayout2D(
    grid, 
    1,                                      # Number of components
    pcms.CoordinateSystem.Cartesian
)
ug_field = ug_layout.create_field()
```

---

### Step 7: Interpolate Field Data

Interpolate field values from the unstructured mesh to the structured uniform grid vertices.

```python
# Transfer field from Omega_h mesh to uniform grid
pcms.interpolate_field(omega_h_field, ug_field)
```

---

### Step 8: Access Field Data

Retrieve interpolated values and verify correctness.

```python
# Get field data and coordinates
ug_field_data = ug_field.get_dof_holder_data()
ug_coords = ug_layout.get_dof_holder_coordinates()

# Access data at vertex (i, j)
vertex_id = j * (grid.divisions[0] + 1) + i
value = ug_field_data[vertex_id]
x, y = ug_coords[vertex_id, 0], ug_coords[vertex_id, 1]
```

---

### Step 9: Export Field Data with `to_mdspan`

Convert field data to a structured $x \times y$ (or $x \times y \times z$) array
for downstream analyses.

```python
# Convert field data to a structured array
grid_values = ug_field.to_numpy()  # 2D: (nx+1, ny+1), 3D: (nx+1, ny+1, nz+1)
mask_values = mask_field.to_numpy()  # 2D: (nx+1, ny+1), 3D: (nx+1, ny+1, nz+1)

# Save to file (example)
np.save('field_data.npy', grid_values)
```

---

## API Reference Summary

### Grid Creation
- `create_uniform_grid_from_mesh(mesh, divisions)` - Create grid from mesh
- `create_uniform_grid_binary_field(mesh, divisions)` - Create inside/outside mask

### Field Layout
- `create_lagrange_layout(mesh, order, num_components, coord_system)` - Omega_h field layout
- `UniformGridFieldLayout2D(grid, num_components, coord_system)` - 2D grid layout
- `UniformGridFieldLayout3D(grid, num_components, coord_system)` - 3D grid layout

### Field Operations
- `layout.create_field()` - Create field from layout
- `field.set_dof_holder_data(data)` - Set field values
- `field.get_dof_holder_data()` - Get field values
- `field.to_mdspan()` - Get field values as a structured array
- `field.set_out_of_bounds_mode(mode, fill_value=0.0)` - Set behavior for points outside mesh
- `interpolate_field(source_field, target_field)` - Interpolate between fields
- `map_entity_field_to_vertices_average(mesh, field_data, entity_dim)` - Convert element/face field to vertex field by averaging

### Out-of-Bounds Modes
- `OutOfBoundsMode.ERROR` - Raise error when points are outside mesh
- `OutOfBoundsMode.FILL` - Fill with specified value when points are outside mesh
- `OutOfBoundsMode.NEAREST_BOUNDARY` - Clamp to nearest boundary cell (extrapolate)

### Mesh Tags
- `mesh.ntags(dim)` - Number of tags on entities of dimension `dim`
- `mesh.get_tag(dim, index_or_name)` - Get tag by index or name
- `tag.name()` - Tag name
- `tag.ncomps()` - Number of components in tag

### Exodus I/O (if OMEGA_H_USE_SEACASEXODUS enabled)
- `exodus_open(filepath, verbose=False)` - Open Exodus file, returns file handle
- `exodus_close(exodus_file)` - Close Exodus file handle
- `exodus_get_num_time_steps(exodus_file)` - Get number of time steps in file
- `read_mesh_exodus(exodus_file, mesh, verbose=False, classify_with=...)` - Read mesh from file handle
- `read_exodus_nodal_fields(exodus_file, mesh, time_step, prefix="", postfix="", verbose=False)` - Read nodal fields
- `read_exodus_element_fields(exodus_file, mesh, time_step, prefix="", postfix="", verbose=False)` - Read element fields
- `read_mesh_exodus_sliced(filepath, comm, verbose=False, classify_with=..., time_step=-1)` - Read mesh in parallel
- `write_mesh_exodus(filepath, mesh, verbose=False, classify_with=...)` - Write mesh to Exodus file
- `ExodusClassifyWith.NODE_SETS` - Classify using node sets
- `ExodusClassifyWith.SIDE_SETS` - Classify using side sets

### MESHB I/O (if OMEGA_H_USE_LIBMESHB enabled)
- `read_mesh_meshb(mesh, filepath)` - Read mesh from MESHB file
- `write_mesh_meshb(mesh, filepath, version)` - Write mesh to MESHB file
- `read_meshb_sol(mesh, filepath, sol_name)` - Read solution/resolution from .sol file
- `write_meshb_sol(mesh, filepath, sol_name, version)` - Write solution/resolution to .sol file

### Grid Properties
- `grid.get_num_cells()` - Total number of cells
- `grid.divisions` - Cell divisions in each dimension
- `layout.get_num_vertices()` - Total number of vertices
- `layout.get_dof_holder_coordinates()` - Vertex coordinates

---

For more examples and advanced usage, see the test files in the `pcms/pythonapi/` directory.
