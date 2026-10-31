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
import PyOmega_h as omega_h
import numpy as np

# Create library and world communicator
lib = omega_h.OmegaHLibrary()
world = lib.world()
```


---

### Step 2: Create an Omega_h Mesh

You can either create a mesh or read an existing mesh from a file.

#### Option A: Create Mesh Programmatically

```python
# Create a 2D box mesh: 1.0 x 1.0 domain with 4x4 elements
mesh = omega_h.build_box(
    world,                          # Communicator
    omega_h.Family.SIMPLEX,        # Element type
    1.0, 1.0, 0.0,                 # Domain size: x, y, z
    4, 4, 0,                       # Number of elements: nx, ny, nz
    False                          # Symmetric flag
)
```

#### Option B: Read Mesh from File

```python
# Auto-detect file format and read
mesh = omega_h.read_mesh_file("my_mesh.osh", world)

# Or use format-specific readers for explicit control:

# Binary format (.osh)
mesh = omega_h.read_mesh_binary("mesh.osh", lib)

# Gmsh format (.msh)
mesh = omega_h.read_mesh_gmsh("mesh.msh", world)

# VTK format (.pvtu for parallel files)
mesh = omega_h.read_mesh_parallel_vtk("mesh.pvtu", world)

# Exodus format (.exo) - if SEACAS support is enabled
# Using file handle API
exo_handle = omega_h.exodus_open("mesh.exo", verbose=False)
num_steps = omega_h.exodus_get_num_time_steps(exo_handle)
mesh = omega_h.Mesh(lib)
mesh.set_comm(world)
omega_h.read_mesh_exodus(exo_handle, mesh, verbose=False)
omega_h.exodus_close(exo_handle)

# ADIOS2 format (.bp) - if ADIOS2 support is enabled
mesh = omega_h.read_mesh_adios2("mesh.bp", lib, prefix="")

# MESHB format (.mesh) - if LIBMESHB support is enabled
mesh = omega_h.OmegaHMesh(lib)
mesh.set_comm(world)
omega_h.read_mesh_meshb(mesh, "mesh.mesh")
# Read solution/resolution data if available
omega_h.read_meshb_sol(mesh, "mesh.sol", sol_name="resolution")
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
mask_field = pcms.create_uniform_grid_binary_field(mesh, [4, 4])
print(f"Created mask field with {mask_field.get_layout().get_num_owned_dof_holder()} vertices")
```

**Accessing Values**:
```python
# Access mask value for vertex (i, j)
vertex_id = j * (grid.divisions[0] + 1) + i
mask_data = mask_field.get_dof_holder_data().get_values()  # shape (num_vertices, 1)
mask_value = mask_data[vertex_id, 0]  # Returns 0.0 or 1.0
```

---

### Step 5: Create and Initialize Omega_h Field

#### Option A: Create Field Programmatically

```python
# Create a concrete FunctionSpace with linear (order=1) elements and
# 1 component (scalar)
omega_h_factory = pcms.LagrangeFunctionSpace.from_mesh(
    mesh,
    1,
    1,
)
omega_h_field = omega_h_factory.create_function()

# Initialize field with f(x,y) = x + 2*y
coords = omega_h_field.get_layout().get_dof_holder_coordinates().get_values()
num_nodes = omega_h_field.get_layout().get_num_owned_dof_holder()
omega_h_data = np.zeros(num_nodes)

for i in range(num_nodes):
    x = coords[i, 0]
    y = coords[i, 1]
    omega_h_data[i] = x + 2.0 * y

omega_h_field.set_dof_holder_data(omega_h_data)

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

# Create vertex-based field and set data
omega_h_factory = pcms.LagrangeFunctionSpace.from_mesh(
    mesh, 1, 1)
omega_h_field = omega_h_factory.create_function()
omega_h_field.set_dof_holder_data(vertex_field)

```

---

### Step 6: Create Uniform Grid Field

Set up the data structure for storing field values on the uniform grid vertices.

```python
# Create a uniform-grid FunctionSpace
ug_factory = pcms.LagrangeFunctionSpace.from_uniform_grid(
    grid,
    1,                                      # Number of components
)
ug_field = ug_factory.create_function()
```

---

### Step 7: Interpolate Field Data

Interpolate field values from the unstructured mesh to the structured uniform grid vertices.

```python
# Transfer field from Omega_h mesh to uniform grid using two FunctionSpaces
interp = pcms.Interpolator(omega_h_factory, ug_factory)
interp.apply(omega_h_field, ug_field)
```

---

### Step 8: Access Field Data

Retrieve interpolated values and verify correctness.

```python
# Get field data and coordinates
ug_field_data = ug_field.get_dof_holder_data().get_values()
ug_coords = ug_field.get_layout().get_dof_holder_coordinates().get_values()

# Access data at vertex (i, j)
vertex_id = j * (grid.divisions[0] + 1) + i
value = ug_field_data[vertex_id, 0]
x, y = ug_coords[vertex_id, 0], ug_coords[vertex_id, 1]
```

---

### Step 9: Export Field Data

Retrieve the DOF arrays for downstream analyses.

```python
# Get DOF arrays, shape (num_dof_holders, num_components)
grid_values = ug_field.get_dof_holder_data().get_values()
mask_values = mask_field.get_dof_holder_data().get_values()

# Save to file (example)
np.save('field_data.npy', grid_values)
```

---

## API Reference Summary

### Grid Creation
- `create_uniform_grid_from_mesh(mesh, divisions)` - Create grid from mesh
- `create_uniform_grid_binary_field(mesh, divisions)` - Create inside/outside mask

The Python names mirror the C++ API: classes keep their C++ names, C++
namespaces are submodules (`pcms.csys`, `pcms.values`), and methods are the
snake_case form of the C++ method (`GetDOFHolderData` is
`get_dof_holder_data`).

### Coordinate Systems and Value Bases
- `csys.Cartesian.create(dim)` - Cartesian system. Coordinate system arguments default to `None`, meaning Cartesian in the dimension of the mesh, grid or coordinates they are used with
- `csys.CylindricalRThetaZ.create()`, `csys.CylindricalRZ.create()` - cylindrical systems
- `CoordinateView(system, coords)`, `CoordinateView(coords)` - coordinates `(num_points, dim)` tagged with their system (Cartesian when omitted); the array is copied
- `view.get_coordinate_system()`, `view.get_values()` - the system and a read-only coordinate array
- `ValueBasis(system, scaling, variance)` - the basis a value's components are in; `ValueBasis()` is the scalar basis
- `values.Scalar`, `values.Vector`, `values.Covector`, `values.Tensor`, `values.Of([...])` - variance signatures
- `ComponentScaling.Physical`, `ComponentScaling.Natural`, `Variance.Contravariant`, `Variance.Covariant`
- `ValueView(basis, values)` - values `(num_points, num_components)` tagged with their basis; the array is not copied

### Field Factories
- `FunctionSpace` - abstract base class for field spaces; a `PointEvaluatorFactory`
- `LagrangeFunctionSpace.from_mesh(mesh, order, num_components, mesh_coordinate_system)` - create an Omega_h-backed FunctionSpace
- `LagrangeFunctionSpace.from_uniform_grid(grid, num_components, mesh_coordinate_system, order=1)` - create a uniform-grid-backed FunctionSpace
- `PolynomialReconstructionFunctionSpace.create(coordinate_view, options)` - MLS reconstruction over a point cloud

### Field Operations
- `space.create_function(name="", variance=values.Scalar, system=None, scaling=None)` - Create a real-valued `Function`; vector and tensor values are in `system`'s basis (default: the space's)
- `field.get_layout()` - the field's `FieldLayout`
- `field.get_data().get_value_basis()` - the basis the field's values are stored in
- `field.get_layout().get_num_owned_dof_holder()` - Number of owned DOF holders (nodes/elements)
- `field.get_layout().get_num_components()` - Number of field components per DOF holder
- `field.get_layout().get_dof_holder_coordinates()` - DOF holder coordinates as a `CoordinateView`
- `field.get_dof_holder_data()` - Field values as a `ValueView` whose array (`.get_values()`, shape `(num_dof_holders, num_components)`) is read-only and valid until the field is next modified (use `.copy()` to keep values)
- `field.set_dof_holder_data(value_view)` - Set field values; the basis must match the field's. Scalar fields also accept a plain array, of shape `(num_dof_holders,)` when the field has one component
- `field.set_dof_holder_data_unchecked(data)` - Set field values without checking their basis
- `EvaluationRequest.from_coordinates(coordinate_view, policy=...)` - Build an explicit coordinate-based evaluation request
- `EvaluationRequest.from_function_space(space, policy=...)` - Build an evaluation request from another FunctionSpace's DOF-holder sites
- `factory.create_point_evaluator(request)` - Create a reusable point evaluator; the request's points must be in `factory.get_coordinate_system()`
- `evaluator.evaluate(field, output)` - Evaluate into a `ValueView` tagged `evaluator.output_basis(field.get_data().get_value_basis())` (a plain array for scalar fields, of shape `(num_points,)` when the output has one component); the output array must be writeable, C-contiguous float64
- `Interpolator(source, target_space)` - Create an interpolator from a `PointEvaluatorFactory` (cached localization)
- `operator.apply(source_field, target_field)` - Apply any `TransferOperator` (`Interpolator`, `Copy`, projections)
- `Copy(source_space, target_space)` - Create a copy operator for compatible FunctionSpaces
- `map_entity_field_to_vertices_average(mesh, field_data, entity_dim)` - Convert element/face field to vertex field by averaging

### Coordinate Transformations
- `CartesianToCylindrical()`, `CylindricalToCartesian()` - `CoordinateMap`s; points map source to target, values are re-expressed the other way
- `bound = coordinate_map.bind(coordinate_view)` - per-point work for points in the map's source system
- `bound.mapped_points()`, `bound.output_basis(basis)`, `bound.transform_values(value_view)` - mapped points, and in-place re-expression of values
- `TransformedSource(source, to_source)` - presents a `PointEvaluatorFactory` in another coordinate system; use it as an `Interpolator` source

### Out-of-Bounds Modes
- `OutOfBoundsMode.ERROR` - Raise error when points are outside mesh
- `OutOfBoundsMode.FILL` - Fill with specified value when points are outside mesh
- `OutOfBoundsMode.NEAREST_BOUNDARY` - Clamp to nearest boundary cell (extrapolate)

### Mesh Tags
- `mesh.ntags(dim)` - Number of tags on entities of dimension `dim`
- `mesh.get_tag(dim, index_or_name)` - Get tag by index or name
- `tag.name()` - Tag name
- `tag.ncomps()` - Number of components in tag

### Omega_h I/O (available via `omega_h.*` from `import PyOmega_h as omega_h`)
- `omega_h.read_mesh_file(filepath, comm)` - Auto-detect format and read mesh
- `omega_h.read_mesh_binary(filepath, lib)` - Read binary .osh mesh
- `omega_h.write_mesh_binary(filepath, mesh)` - Write binary .osh mesh
- `omega_h.read_mesh_gmsh(filepath, comm)` - Read Gmsh .msh mesh
- `omega_h.write_mesh_gmsh(filepath, mesh)` - Write Gmsh .msh mesh
- `omega_h.read_mesh_parallel_vtk(pvtupath, comm)` - Read parallel VTK .pvtu mesh
- `omega_h.write_mesh_vtu(filepath, mesh, compress=...)` - Write VTK .vtu mesh

### Exodus I/O (if OMEGA_H_USE_SEACASEXODUS enabled)
- `omega_h.exodus_open(filepath, verbose=False)` - Open Exodus file, returns file handle
- `omega_h.exodus_close(exodus_file)` - Close Exodus file handle
- `omega_h.exodus_get_num_time_steps(exodus_file)` - Get number of time steps in file
- `omega_h.read_mesh_exodus(exodus_file, mesh, verbose=False, classify_with=...)` - Read mesh from file handle
- `omega_h.read_exodus_nodal_fields(exodus_file, mesh, time_step, prefix="", postfix="", verbose=False)` - Read nodal fields
- `omega_h.read_exodus_element_fields(exodus_file, mesh, time_step, prefix="", postfix="", verbose=False)` - Read element fields
- `omega_h.read_mesh_exodus_sliced(filepath, comm, verbose=False, classify_with=..., time_step=-1)` - Read mesh in parallel
- `omega_h.write_mesh_exodus(filepath, mesh, verbose=False, classify_with=...)` - Write mesh to Exodus file
- `omega_h.ExodusClassifyWith.NODE_SETS` - Classify using node sets
- `omega_h.ExodusClassifyWith.SIDE_SETS` - Classify using side sets

### MESHB I/O (if OMEGA_H_USE_LIBMESHB enabled)
- `omega_h.read_mesh_meshb(mesh, filepath)` - Read mesh from MESHB file
- `omega_h.write_mesh_meshb(mesh, filepath, version)` - Write mesh to MESHB file
- `omega_h.read_meshb_sol(mesh, filepath, sol_name)` - Read solution/resolution from .sol file
- `omega_h.write_meshb_sol(mesh, filepath, sol_name, version)` - Write solution/resolution to .sol file

### Grid Properties
- `grid.get_num_cells()` - Total number of cells
- `grid.divisions` - Cell divisions in each dimension
- `field.get_layout().get_num_owned_dof_holder()` - Total number of field DOF holders
- `field.get_layout().get_dof_holder_coordinates().get_values()` - DOF holder coordinates

---

For more examples and advanced usage, see the test files in the `pcms/pythonapi/` directory.
