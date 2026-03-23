#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#include <Omega_h_mesh.hpp>
#include <Omega_h_file.hpp>
#include <Omega_h_library.hpp>
#include <Omega_h_build.hpp>
#include <Omega_h_comm.hpp>
#include <Omega_h_array_ops.hpp>
#include <Omega_h_defines.hpp>
#include <Omega_h_graph.hpp>
#include <Omega_h_tag.hpp>
#ifdef OMEGA_H_USE_ADIOS2
#include <Omega_h_adios2.hpp>
#endif
#include "numpy_array_transform.h"

namespace py = pybind11;

namespace pcms
{

void bind_omega_h_mesh_module(py::module& m)
{
  // Bind Omega_h::Library
  py::class_<Omega_h::Library, std::shared_ptr<Omega_h::Library>>(
    m, "OmegaHLibrary")
    .def(py::init<>(), "Default constructor")
    .def("world", &Omega_h::Library::world, py::return_value_policy::reference,
         "Get the world communicator");

  // Bind the Comm class
#ifdef OMEGA_H_USE_MPI
  using CommHandle = std::uintptr_t;
#endif

  py::class_<Omega_h::Comm, std::shared_ptr<Omega_h::Comm>>(m, "Comm")
  // Constructors
#ifdef OMEGA_H_USE_MPI
    .def(py::init([](Omega_h::Library* library, CommHandle impl_handle) {
           MPI_Comm impl;
           std::memcpy(&impl, &impl_handle, sizeof(MPI_Comm));
           return new Omega_h::Comm(library, impl);
         }),
         py::arg("library"), py::arg("impl_handle"))

    .def(py::init([](Omega_h::Library* library, CommHandle impl_handle,
                     py::array_t<const Omega_h::I32> srcs,
                     py::array_t<const Omega_h::I32> dsts) {
           MPI_Comm impl;
           std::memcpy(&impl, &impl_handle, sizeof(MPI_Comm));
           auto srcs_view = numpy_to_omega_h_read<Omega_h::I32>(srcs);
           auto dsts_view = numpy_to_omega_h_read<Omega_h::I32>(dsts);
           return new Omega_h::Comm(library, impl, srcs_view, dsts_view);
         }),
         py::arg("library"), py::arg("impl_handle"), py::arg("srcs"),
         py::arg("dsts"))

    .def(
      "get_impl_handle",
      [](Omega_h::Comm const& self) -> CommHandle {
        MPI_Comm impl = self.get_impl();
        CommHandle handle = 0;
        std::memcpy(&handle, &impl, sizeof(MPI_Comm));
        return handle;
      },
      "Get the underlying MPI communicator as an opaque integer handle")
#else
    .def(py::init<Omega_h::Library*, bool, bool>(), py::arg("library"),
         py::arg("is_graph"), py::arg("sends_to_self"))
#endif
    // Methods
    .def("library", &Omega_h::Comm::library, py::return_value_policy::reference,
         "Get the library pointer")
    .def("rank", &Omega_h::Comm::rank, "Get the rank of this process")
    .def("size", &Omega_h::Comm::size, "Get the total number of processes")
    .def("dup", &Omega_h::Comm::dup, "Duplicate the communicator")
    .def("split", &Omega_h::Comm::split, py::arg("color"), py::arg("key"),
         "Split the communicator")
    .def("graph", &Omega_h::Comm::graph, py::arg("dsts"),
         "Create a graph communicator")
    .def("graph_adjacent", &Omega_h::Comm::graph_adjacent, py::arg("srcs"),
         py::arg("dsts"), "Create an adjacent graph communicator")
    .def("graph_inverse", &Omega_h::Comm::graph_inverse,
         "Get the inverse graph communicator")
    .def("sources", &Omega_h::Comm::sources, "Get source ranks")
    .def("destinations", &Omega_h::Comm::destinations, "Get destination ranks")
    .def("reduce_or", &Omega_h::Comm::reduce_or, py::arg("x"),
         "Reduce using logical OR")
    .def("reduce_and", &Omega_h::Comm::reduce_and, py::arg("x"),
         "Reduce using logical AND")
    .def("add_int128", &Omega_h::Comm::add_int128, py::arg("x"),
         "Add Int128 values across processes")
    .def("bcast_string", &Omega_h::Comm::bcast_string, py::arg("s"),
         py::arg("root_rank") = 0, "Broadcast a string")
    .def("barrier", &Omega_h::Comm::barrier, "Synchronize all processes");
  m.def("build_box", &Omega_h::build_box, py::arg("comm"), py::arg("family"),
        py::arg("x"), py::arg("y"), py::arg("z"), py::arg("nx"), py::arg("ny"),
        py::arg("nz"), py::arg("symmetric") = false,
        "Build a box mesh with specified dimensions and discretization",
        py::return_value_policy::move);

  // Bind Omega_h_Type enum for tag data types
  py::enum_<Omega_h_Type>(m, "OmegaHType")
    .value("I8", OMEGA_H_I8)
    .value("I32", OMEGA_H_I32)
    .value("I64", OMEGA_H_I64)
    .value("F64", OMEGA_H_F64)
    .value("REAL", OMEGA_H_REAL)
    .export_values();

  // Bind Omega_h_Op enum for reduction operations
  py::enum_<Omega_h_Op>(m, "OmegaHOp")
    .value("MIN", OMEGA_H_MIN)
    .value("MAX", OMEGA_H_MAX)
    .value("SUM", OMEGA_H_SUM)
    .export_values();

  // Bind ArrayType enum for tag array types
  py::enum_<Omega_h::ArrayType>(m, "ArrayType")
    .value("VectorND", Omega_h::ArrayType::VectorND)
    .value("SymmetricSquareMatrix", Omega_h::ArrayType::SymmetricSquareMatrix)
    .export_values();

  // Bind TagBase class for tag metadata access
  py::class_<Omega_h::TagBase>(m, "TagBase")
    .def("name", &Omega_h::TagBase::name, py::return_value_policy::reference,
         "Get tag name")
    .def("ncomps", &Omega_h::TagBase::ncomps, "Get number of components")
    .def("type", &Omega_h::TagBase::type, "Get tag data type (Omega_h_Type)")
    .def("array_type", &Omega_h::TagBase::array_type,
         "Get array type (VectorND or SymmetricSquareMatrix)")
    .def(
      "class_ids",
      [](const Omega_h::TagBase& tag) {
        auto class_ids = tag.class_ids();
        return omega_h_read_to_numpy(class_ids);
      },
      "Get class IDs for rcField tags");

  py::class_<Omega_h::Graph>(m, "OmegaHGraph")
    .def_property(
      "a2ab",
      [](const Omega_h::Graph& graph) {
        return omega_h_read_to_numpy(graph.a2ab);
      },
      [](Omega_h::Graph& graph, py::array_t<Omega_h::LO> values) {
        graph.a2ab = numpy_to_omega_h_read<Omega_h::LO>(values);
      })
    .def_property(
      "ab2b",
      [](const Omega_h::Graph& graph) {
        return omega_h_read_to_numpy(graph.ab2b);
      },
      [](Omega_h::Graph& graph, py::array_t<Omega_h::LO> values) {
        graph.ab2b = numpy_to_omega_h_read<Omega_h::LO>(values);
      })
    .def("nnodes", &Omega_h::Graph::nnodes, "Get number of graph nodes");

  // Bind Omega_h::Mesh
  py::class_<Omega_h::Mesh, std::shared_ptr<Omega_h::Mesh>>(m, "OmegaHMesh")
    .def(py::init<>(), "Default constructor")
    .def(py::init<Omega_h::Library*>(), py::arg("library"),
         "Constructor with library")

    .def("set_library", &Omega_h::Mesh::set_library, py::arg("library"),
         "Set the library")

    .def("library", &Omega_h::Mesh::library, py::return_value_policy::reference,
         "Get the library")

    .def("set_comm", &Omega_h::Mesh::set_comm, py::arg("comm"),
         "Set the communicator")

    .def("comm", &Omega_h::Mesh::comm, "Get the communicator")

    .def("set_dim", &Omega_h::Mesh::set_dim, py::arg("dim"),
         "Set mesh dimension")

    .def("dim", &Omega_h::Mesh::dim, "Get mesh dimension")

    .def("set_family", &Omega_h::Mesh::set_family, py::arg("family"),
         "Set mesh family (simplex, hypercube, etc.)")

    .def("family", &Omega_h::Mesh::family, "Get mesh family")

    .def("nverts", &Omega_h::Mesh::nverts, "Get number of vertices")

    .def("nedges", &Omega_h::Mesh::nedges, "Get number of edges")

    .def("nfaces", &Omega_h::Mesh::nfaces, "Get number of faces")

    .def("nregions", &Omega_h::Mesh::nregions, "Get number of regions")

    .def("nelems", &Omega_h::Mesh::nelems, "Get number of elements")

    .def("nents", &Omega_h::Mesh::nents, py::arg("ent_dim"),
         "Get number of entities of given dimension")

    .def("nglobal_ents", &Omega_h::Mesh::nglobal_ents, py::arg("ent_dim"),
         "Get global number of entities")

    .def(
      "coords",
      [](const Omega_h::Mesh& mesh) {
        auto coords = mesh.coords();
        return omega_h_read_to_numpy(coords);
      },
      "Get mesh coordinates as numpy array")

    .def(
      "set_coords",
      [](Omega_h::Mesh& mesh, py::array_t<Omega_h::Real> coords) {
        auto coords_write = numpy_to_omega_h_write<Omega_h::Real>(coords);
        mesh.set_coords(Omega_h::Reals(coords_write));
      },
      py::arg("coords"), "Set mesh coordinates from numpy array")

    .def(
      "add_coords",
      [](Omega_h::Mesh& mesh, py::array_t<Omega_h::Real> coords) {
        auto coords_write = numpy_to_omega_h_write<Omega_h::Real>(coords);
        mesh.add_coords(Omega_h::Reals(coords_write));
      },
      py::arg("coords"), "Add mesh coordinates from numpy array")

    .def("has_tag", &Omega_h::Mesh::has_tag, py::arg("ent_dim"),
         py::arg("name"), "Check if mesh has a tag")

    .def("ntags", &Omega_h::Mesh::ntags, py::arg("ent_dim"),
         "Get number of tags for entity dimension")

    .def("nrctags", &Omega_h::Mesh::nrctags, py::arg("ent_dim"),
         "Get number of rcField tags for entity dimension")

    .def(
      "get_tag",
      [](const Omega_h::Mesh& mesh, Omega_h::Int dim, Omega_h::Int i)
        -> const Omega_h::TagBase* { return mesh.get_tag(dim, i); },
      py::arg("ent_dim"), py::arg("index"), py::return_value_policy::reference,
      "Get tag by index (for iterating through tags)")

    .def(
      "get_tagbase",
      [](const Omega_h::Mesh& mesh, Omega_h::Int dim, const std::string& name)
        -> const Omega_h::TagBase* { return mesh.get_tagbase(dim, name); },
      py::arg("ent_dim"), py::arg("name"), py::return_value_policy::reference,
      "Get tag metadata by name")

    .def("remove_tag", &Omega_h::Mesh::remove_tag, py::arg("ent_dim"),
         py::arg("name"), "Remove a tag")

    // Tag creation with dtype string
    .def(
      "add_tag",
      [](Omega_h::Mesh& mesh, Omega_h::Int dim, const std::string& name,
         Omega_h::Int ncomps, const std::string& dtype,
         Omega_h::ArrayType array_type) {
        if (dtype == "int8" || dtype == "i1") {
          mesh.add_tag<Omega_h::I8>(dim, name, ncomps, array_type);
        } else if (dtype == "int32" || dtype == "i4") {
          mesh.add_tag<Omega_h::I32>(dim, name, ncomps, array_type);
        } else if (dtype == "int64" || dtype == "i8") {
          mesh.add_tag<Omega_h::I64>(dim, name, ncomps, array_type);
        } else if (dtype == "float64" || dtype == "f8" || dtype == "double") {
          mesh.add_tag<Omega_h::Real>(dim, name, ncomps, array_type);
        } else {
          throw std::runtime_error(
            "Unsupported dtype: " + dtype +
            ". Use 'int8', 'int32', 'int64', or 'float64'");
        }
      },
      py::arg("ent_dim"), py::arg("name"), py::arg("ncomps"),
      py::arg("dtype") = "float64",
      py::arg("array_type") = Omega_h::ArrayType::VectorND,
      "Add a tag with specified dtype (int8, int32, int64, float64)")

    // Add tag with numpy array (auto-detect type)
    .def(
      "add_tag",
      [](Omega_h::Mesh& mesh, Omega_h::Int dim, const std::string& name,
         Omega_h::Int ncomps, py::array array, bool internal,
         Omega_h::ArrayType array_type) {
        auto dtype = array.dtype();
        if (dtype.is(py::dtype::of<std::int8_t>())) {
          auto array_view = numpy_to_omega_h_read<Omega_h::I8>(
            array.cast<py::array_t<std::int8_t>>());
          mesh.add_tag<Omega_h::I8>(dim, name, ncomps, array_view, internal,
                                    array_type);
        } else if (dtype.is(py::dtype::of<std::int32_t>())) {
          auto array_view = numpy_to_omega_h_read<Omega_h::I32>(
            array.cast<py::array_t<std::int32_t>>());
          mesh.add_tag<Omega_h::I32>(dim, name, ncomps, array_view, internal,
                                     array_type);
        } else if (dtype.is(py::dtype::of<std::int64_t>())) {
          auto array_view = numpy_to_omega_h_read<Omega_h::I64>(
            array.cast<py::array_t<std::int64_t>>());
          mesh.add_tag<Omega_h::I64>(dim, name, ncomps, array_view, internal,
                                     array_type);
        } else if (dtype.is(py::dtype::of<double>())) {
          auto array_view = numpy_to_omega_h_read<Omega_h::Real>(
            array.cast<py::array_t<double>>());
          mesh.add_tag<Omega_h::Real>(dim, name, ncomps, array_view, internal,
                                      array_type);
        } else {
          throw std::runtime_error(
            "Unsupported numpy dtype. Use int8, int32, int64, or float64");
        }
      },
      py::arg("ent_dim"), py::arg("name"), py::arg("ncomps"), py::arg("array"),
      py::arg("internal") = false,
      py::arg("array_type") = Omega_h::ArrayType::VectorND,
      "Add a tag with initial values from numpy array (auto-detects type)")

    // Get tag data (auto-detect type from stored tag)
    .def(
      "get_tag",
      [](Omega_h::Mesh& mesh, Omega_h::Int dim,
         const std::string& name) -> py::array {
        auto tagbase = mesh.get_tagbase(dim, name);
        auto type = tagbase->type();

        if (type == OMEGA_H_I8) {
          auto array = mesh.get_array<Omega_h::I8>(dim, name);
          return omega_h_read_to_numpy(array);
        } else if (type == OMEGA_H_I32) {
          auto array = mesh.get_array<Omega_h::I32>(dim, name);
          return omega_h_read_to_numpy(array);
        } else if (type == OMEGA_H_I64) {
          auto array = mesh.get_array<Omega_h::I64>(dim, name);
          return omega_h_read_to_numpy(array);
        } else if (type == OMEGA_H_F64) {
          auto array = mesh.get_array<Omega_h::Real>(dim, name);
          return omega_h_read_to_numpy(array);
        } else {
          throw std::runtime_error("Unsupported tag type");
        }
      },
      py::arg("ent_dim"), py::arg("name"),
      "Get tag data as numpy array by name (auto-detects type)")

    .def(
      "get_tag_array",
      [](Omega_h::Mesh& mesh, Omega_h::Int dim,
         const std::string& name) -> py::array {
        auto tagbase = mesh.get_tagbase(dim, name);
        auto type = tagbase->type();

        if (type == OMEGA_H_I8) {
          auto array = mesh.get_array<Omega_h::I8>(dim, name);
          return omega_h_read_to_numpy(array);
        } else if (type == OMEGA_H_I32) {
          auto array = mesh.get_array<Omega_h::I32>(dim, name);
          return omega_h_read_to_numpy(array);
        } else if (type == OMEGA_H_I64) {
          auto array = mesh.get_array<Omega_h::I64>(dim, name);
          return omega_h_read_to_numpy(array);
        } else if (type == OMEGA_H_F64) {
          auto array = mesh.get_array<Omega_h::Real>(dim, name);
          return omega_h_read_to_numpy(array);
        } else {
          throw std::runtime_error("Unsupported tag type");
        }
      },
      py::arg("ent_dim"), py::arg("name"),
      "Get tag data as numpy array by name (alias for clarity)")

    // Set tag data from numpy array (auto-detect type)
    .def(
      "set_tag",
      [](Omega_h::Mesh& mesh, Omega_h::Int dim, const std::string& name,
         py::array array, bool internal, Omega_h::ArrayType array_type) {
        auto dtype = array.dtype();
        if (dtype.is(py::dtype::of<std::int8_t>())) {
          auto array_view = numpy_to_omega_h_read<Omega_h::I8>(
            array.cast<py::array_t<std::int8_t>>());
          mesh.set_tag<Omega_h::I8>(dim, name, array_view, internal,
                                    array_type);
        } else if (dtype.is(py::dtype::of<std::int32_t>())) {
          auto array_view = numpy_to_omega_h_read<Omega_h::I32>(
            array.cast<py::array_t<std::int32_t>>());
          mesh.set_tag<Omega_h::I32>(dim, name, array_view, internal,
                                     array_type);
        } else if (dtype.is(py::dtype::of<std::int64_t>())) {
          auto array_view = numpy_to_omega_h_read<Omega_h::I64>(
            array.cast<py::array_t<std::int64_t>>());
          mesh.set_tag<Omega_h::I64>(dim, name, array_view, internal,
                                     array_type);
        } else if (dtype.is(py::dtype::of<double>())) {
          auto array_view = numpy_to_omega_h_read<Omega_h::Real>(
            array.cast<py::array_t<double>>());
          mesh.set_tag<Omega_h::Real>(dim, name, array_view, internal,
                                      array_type);
        } else {
          throw std::runtime_error(
            "Unsupported numpy dtype. Use int8, int32, int64, or float64");
        }
      },
      py::arg("ent_dim"), py::arg("name"), py::arg("array"),
      py::arg("internal") = false,
      py::arg("array_type") = Omega_h::ArrayType::VectorND,
      "Set tag data from numpy array (auto-detects type)")

    // Synchronize and reduce tags in parallel
    .def("sync_tag", &Omega_h::Mesh::sync_tag, py::arg("ent_dim"),
         py::arg("name"), "Synchronize tag values across processors")

    .def("reduce_tag", &Omega_h::Mesh::reduce_tag, py::arg("ent_dim"),
         py::arg("name"), py::arg("op"),
         "Reduce tag values across processors with specified operation")

    .def("has_ents", &Omega_h::Mesh::has_ents, py::arg("ent_dim"),
         "Check if mesh has entities of given dimension")

    .def(
      "ask_lengths",
      [](Omega_h::Mesh& mesh) {
        auto lengths = mesh.ask_lengths();
        return omega_h_read_to_numpy(lengths);
      },
      "Get edge lengths")

    .def(
      "ask_qualities",
      [](Omega_h::Mesh& mesh) {
        auto qualities = mesh.ask_qualities();
        return omega_h_read_to_numpy(qualities);
      },
      "Get element qualities")

    .def("min_quality", &Omega_h::Mesh::min_quality,
         "Get minimum element quality")

    .def("max_length", &Omega_h::Mesh::max_length, "Get maximum edge length")

    .def("balance", py::overload_cast<bool>(&Omega_h::Mesh::balance),
         py::arg("predictive") = false, "Balance the mesh across processors")

    .def("set_parting",
         py::overload_cast<Omega_h_Parting, bool>(&Omega_h::Mesh::set_parting),
         py::arg("parting"), py::arg("verbose") = false,
         "Set mesh partitioning")

    .def("parting", &Omega_h::Mesh::parting, "Get mesh partitioning type")

    .def("nghost_layers", &Omega_h::Mesh::nghost_layers,
         "Get number of ghost layers")

    .def(
      "owned",
      [](Omega_h::Mesh& mesh, Omega_h::Int ent_dim) {
        auto owned = mesh.owned(ent_dim);
        return omega_h_read_to_numpy(owned);
      },
      py::arg("ent_dim"), "Get ownership flags for entities")

    // Adjacency queries
    .def(
      "ask_verts_of",
      [](Omega_h::Mesh& mesh, Omega_h::Int dim) {
        auto verts = mesh.ask_verts_of(dim);
        return omega_h_read_to_numpy(verts);
      },
      py::arg("ent_dim"), "Get vertex IDs for entities of given dimension")

    .def(
      "ask_elem_verts",
      [](Omega_h::Mesh& mesh) {
        auto elem_verts = mesh.ask_elem_verts();
        return omega_h_read_to_numpy(elem_verts);
      },
      "Get vertex IDs for elements")

    .def(
      "get_vtx_patches",
      [](Omega_h::Mesh& mesh, Omega_h::LO min_patch_size, Omega_h::LO tgt_dim) {
        return mesh.get_vtx_patches(min_patch_size, tgt_dim);
      },
      py::arg("min_patch_size"), py::arg("tgt_dim"),
      py::return_value_policy::move, "Get vertex-centered patch graph")

    .def("has_adj", &Omega_h::Mesh::has_adj, py::arg("from_dim"),
         py::arg("to_dim"),
         "Check if adjacency information exists from one dimension to another")

    // Global IDs
    .def(
      "globals",
      [](Omega_h::Mesh& mesh, Omega_h::Int dim) {
        auto globals = mesh.globals(dim);
        return omega_h_read_to_numpy(globals);
      },
      py::arg("ent_dim"), "Get global IDs for entities")

    // Mesh sizing queries
    .def(
      "ask_sizes",
      [](Omega_h::Mesh& mesh) {
        auto sizes = mesh.ask_sizes();
        return omega_h_read_to_numpy(sizes);
      },
      "Get element sizes");

  // Mesh utility functions
  m.def(
    "average_field",
    [](Omega_h::Mesh* mesh, Omega_h::Int dim, Omega_h::Int ncomps,
       py::array_t<Omega_h::Real> v2x) {
      auto v2x_view = numpy_to_omega_h_read<Omega_h::Real>(v2x);
      auto result = Omega_h::average_field(mesh, dim, ncomps, v2x_view);
      return omega_h_read_to_numpy(result);
    },
    py::arg("mesh"), py::arg("ent_dim"), py::arg("ncomps"),
    py::arg("vertex_data"),
    "Average vertex field data to entity centers (e.g., get face coordinates "
    "from vertex coordinates)");

  m.def(
    "average_field",
    [](Omega_h::Mesh* mesh, Omega_h::Int dim, py::array_t<Omega_h::LO> a2e,
       Omega_h::Int ncomps, py::array_t<Omega_h::Real> v2x) {
      auto a2e_view = numpy_to_omega_h_read<Omega_h::LO>(a2e);
      auto v2x_view = numpy_to_omega_h_read<Omega_h::Real>(v2x);
      auto result =
        Omega_h::average_field(mesh, dim, a2e_view, ncomps, v2x_view);
      return omega_h_read_to_numpy(result);
    },
    py::arg("mesh"), py::arg("ent_dim"), py::arg("a2e"), py::arg("ncomps"),
    py::arg("vertex_data"),
    "Average vertex field data to subset of entities specified by a2e");

  // Bind mesh I/O functions

  // General mesh file reader (detects format)
  m.def(
    "read_mesh_file",
    [](const std::string& filepath, std::shared_ptr<Omega_h::Comm> comm) {
      return Omega_h::read_mesh_file(filepath, comm);
    },
    py::arg("filepath"), py::arg("comm"),
    "Read mesh from file (auto-detects format)", py::return_value_policy::move);

  // Binary format I/O
  m.def(
    "read_mesh_binary",
    [](const std::string& filepath, Omega_h::Library* lib) {
      Omega_h::Mesh mesh(lib);
      Omega_h::binary::read(filepath, lib->world(), &mesh);
      return mesh;
    },
    py::arg("filepath"), py::arg("library"), "Read mesh from binary file",
    py::return_value_policy::move);

  m.def(
    "read_mesh_binary",
    [](const std::string& filepath, std::shared_ptr<Omega_h::Comm> comm) {
      return Omega_h::binary::read(filepath, comm);
    },
    py::arg("filepath"), py::arg("comm"), "Read mesh from binary file",
    py::return_value_policy::move);

  m.def(
    "write_mesh_binary",
    [](const std::string& filepath, Omega_h::Mesh& mesh) {
      Omega_h::binary::write(filepath, &mesh);
    },
    py::arg("filepath"), py::arg("mesh"), "Write mesh to binary file");

  // Gmsh format I/O
  m.def(
    "read_mesh_gmsh",
    [](const std::string& filepath, std::shared_ptr<Omega_h::Comm> comm) {
      return Omega_h::gmsh::read(filepath, comm);
    },
    py::arg("filepath"), py::arg("comm"), "Read mesh from Gmsh file",
    py::return_value_policy::move);

  m.def(
    "write_mesh_gmsh",
    [](const std::string& filepath, Omega_h::Mesh& mesh) {
      Omega_h::gmsh::write(filepath, &mesh);
    },
    py::arg("filepath"), py::arg("mesh"), "Write mesh to Gmsh file");

#ifdef OMEGA_H_USE_GMSH
  m.def(
    "read_mesh_gmsh_parallel",
    [](const std::string& filepath, std::shared_ptr<Omega_h::Comm> comm) {
      return Omega_h::gmsh::read_parallel(filepath, comm);
    },
    py::arg("filepath"), py::arg("comm"), "Read parallel Gmsh mesh (MSH 4.1+)",
    py::return_value_policy::move);

  m.def(
    "write_mesh_gmsh_parallel",
    [](const std::string& filepath, Omega_h::Mesh& mesh) {
      Omega_h::gmsh::write_parallel(filepath, mesh);
    },
    py::arg("filepath"), py::arg("mesh"), "Write parallel Gmsh mesh (MSH 4.1)");
#endif

  // VTK format I/O
  m.def(
    "write_mesh_vtu",
    [](const std::string& filepath, Omega_h::Mesh& mesh, bool compress) {
      Omega_h::vtk::write_vtu(filepath, &mesh, compress);
    },
    py::arg("filepath"), py::arg("mesh"),
    py::arg("compress") = OMEGA_H_DEFAULT_COMPRESS, "Write mesh to VTU file");

  m.def(
    "write_mesh_vtu",
    [](const std::string& filepath, Omega_h::Mesh& mesh, Omega_h::Int cell_dim,
       bool compress) {
      Omega_h::vtk::write_vtu(filepath, &mesh, cell_dim, compress);
    },
    py::arg("filepath"), py::arg("mesh"), py::arg("cell_dim"),
    py::arg("compress") = OMEGA_H_DEFAULT_COMPRESS,
    "Write mesh to VTU file with specified cell dimension");

  m.def(
    "write_mesh_parallel_vtk",
    [](const std::string& filepath, Omega_h::Mesh& mesh, bool compress) {
      Omega_h::vtk::write_parallel(filepath, &mesh, compress);
    },
    py::arg("filepath"), py::arg("mesh"),
    py::arg("compress") = OMEGA_H_DEFAULT_COMPRESS,
    "Write mesh to parallel VTK files");

  m.def(
    "read_mesh_parallel_vtk",
    [](const std::string& pvtupath, std::shared_ptr<Omega_h::Comm> comm) {
      Omega_h::Mesh mesh;
      Omega_h::vtk::read_parallel(pvtupath, comm, &mesh);
      return mesh;
    },
    py::arg("pvtupath"), py::arg("comm"), "Read parallel VTK mesh",
    py::return_value_policy::move);

#ifdef OMEGA_H_USE_SEACASEXODUS
  // Exodus ClassifyWith enum
  py::enum_<Omega_h::exodus::ClassifyWith>(m, "ExodusClassifyWith")
    .value("NODE_SETS", Omega_h::exodus::NODE_SETS)
    .value("SIDE_SETS", Omega_h::exodus::SIDE_SETS)
    .export_values();

  // Exodus file operations
  m.def(
    "exodus_open",
    [](const std::string& filepath, bool verbose) {
      return Omega_h::exodus::open(filepath, verbose);
    },
    py::arg("filepath"), py::arg("verbose") = false,
    "Open an Exodus file and return file handle");

  m.def(
    "exodus_close",
    [](int exodus_file) { Omega_h::exodus::close(exodus_file); },
    py::arg("exodus_file"), "Close an Exodus file");

  m.def(
    "exodus_get_num_time_steps",
    [](int exodus_file) {
      return Omega_h::exodus::get_num_time_steps(exodus_file);
    },
    py::arg("exodus_file"), "Get the number of time steps in an Exodus file");

  m.def(
    "read_mesh_exodus",
    [](int exodus_file, Omega_h::Mesh& mesh, bool verbose, int classify_with) {
      Omega_h::exodus::read_mesh(exodus_file, &mesh, verbose, classify_with);
    },
    py::arg("exodus_file"), py::arg("mesh"), py::arg("verbose") = false,
    py::arg("classify_with") =
      Omega_h::exodus::NODE_SETS | Omega_h::exodus::SIDE_SETS,
    "Read mesh from an open Exodus file");

  m.def(
    "read_exodus_nodal_fields",
    [](int exodus_file, Omega_h::Mesh& mesh, int time_step,
       const std::string& prefix, const std::string& postfix, bool verbose) {
      Omega_h::exodus::read_nodal_fields(exodus_file, &mesh, time_step, prefix,
                                         postfix, verbose);
    },
    py::arg("exodus_file"), py::arg("mesh"), py::arg("time_step"),
    py::arg("prefix") = "", py::arg("postfix") = "", py::arg("verbose") = false,
    "Read nodal fields from an open Exodus file at a specific time step");

  m.def(
    "read_exodus_element_fields",
    [](int exodus_file, Omega_h::Mesh& mesh, int time_step,
       const std::string& prefix, const std::string& postfix, bool verbose) {
      Omega_h::exodus::read_element_fields(exodus_file, &mesh, time_step,
                                           prefix, postfix, verbose);
    },
    py::arg("exodus_file"), py::arg("mesh"), py::arg("time_step"),
    py::arg("prefix") = "", py::arg("postfix") = "", py::arg("verbose") = false,
    "Read element fields from an open Exodus file at a specific time step");

  m.def(
    "read_mesh_exodus_sliced",
    [](const std::string& filepath, std::shared_ptr<Omega_h::Comm> comm,
       bool verbose, int classify_with, int time_step) {
      return Omega_h::exodus::read_sliced(filepath, comm, verbose,
                                          classify_with, time_step);
    },
    py::arg("filepath"), py::arg("comm"), py::arg("verbose") = false,
    py::arg("classify_with") =
      Omega_h::exodus::NODE_SETS | Omega_h::exodus::SIDE_SETS,
    py::arg("time_step") = -1, "Read sliced Exodus mesh in parallel",
    py::return_value_policy::move);

  m.def(
    "write_mesh_exodus",
    [](const std::string& filepath, Omega_h::Mesh& mesh, bool verbose,
       int classify_with) {
      Omega_h::exodus::write(filepath, &mesh, verbose, classify_with);
    },
    py::arg("filepath"), py::arg("mesh"), py::arg("verbose") = false,
    py::arg("classify_with") =
      Omega_h::exodus::NODE_SETS | Omega_h::exodus::SIDE_SETS,
    "Write mesh to Exodus file");
#endif

#ifdef OMEGA_H_USE_LIBMESHB
  // MESHB format I/O
  m.def(
    "read_mesh_meshb",
    [](Omega_h::Mesh& mesh, const std::string& filepath) {
      Omega_h::meshb::read(&mesh, filepath);
    },
    py::arg("mesh"), py::arg("filepath"), "Read mesh from MESHB file");

  m.def(
    "write_mesh_meshb",
    [](Omega_h::Mesh& mesh, const std::string& filepath, int version) {
      Omega_h::meshb::write(&mesh, filepath, version);
    },
    py::arg("mesh"), py::arg("filepath"), py::arg("version") = 2,
    "Write mesh to MESHB file");

  m.def(
    "read_meshb_sol",
    [](Omega_h::Mesh& mesh, const std::string& filepath,
       const std::string& sol_name) {
      Omega_h::meshb::read_sol(&mesh, filepath, sol_name);
    },
    py::arg("mesh"), py::arg("filepath"), py::arg("sol_name"),
    "Read solution/resolution data from MESHB .sol file");

  m.def(
    "write_meshb_sol",
    [](Omega_h::Mesh& mesh, const std::string& filepath,
       const std::string& sol_name, int version) {
      Omega_h::meshb::write_sol(&mesh, filepath, sol_name, version);
    },
    py::arg("mesh"), py::arg("filepath"), py::arg("sol_name"),
    py::arg("version") = 2,
    "Write solution/resolution data to MESHB .sol file");
#endif

#ifdef OMEGA_H_USE_ADIOS2
  // ADIOS2 format I/O
  m.def(
    "read_mesh_adios2",
    [](const std::string& filepath, Omega_h::Library* lib,
       const std::string& prefix) {
      return Omega_h::adios::read(filepath, lib, prefix);
    },
    py::arg("filepath"), py::arg("library"), py::arg("prefix") = "",
    "Read mesh from ADIOS2 file", py::return_value_policy::move);

  m.def(
    "write_mesh_adios2",
    [](const std::string& filepath, Omega_h::Mesh& mesh,
       const std::string& prefix) {
      Omega_h::adios::write(filepath, &mesh, prefix);
    },
    py::arg("filepath"), py::arg("mesh"), py::arg("prefix") = "",
    "Write mesh to ADIOS2 file");
#endif

  // Other utilities
  py::enum_<Omega_h_Family>(m, "Family")
    .value("SIMPLEX", OMEGA_H_SIMPLEX)
    .value("HYPERCUBE", OMEGA_H_HYPERCUBE)
    .value("MIXED", OMEGA_H_MIXED)
    .export_values(); // This allows using values without the class prefix

  // Bind Omega_h_Parting enum
  py::enum_<Omega_h_Parting>(m, "OmegaHParting")
    .value("ELEM_BASED", OMEGA_H_ELEM_BASED)
    .value("VERT_BASED", OMEGA_H_VERT_BASED)
    .value("GHOSTED", OMEGA_H_GHOSTED)
    .export_values();
}

} // namespace pcms