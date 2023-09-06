#ifndef WDM_COUPLING_OMEGA_H_FIELD_H
#define WDM_COUPLING_OMEGA_H_FIELD_H
#include "pcms/types.h"
#include "pcms/external/span.h"
#include <Omega_h_mesh.hpp>
#include "pcms/field.h"
#include "pcms/coordinate_systems.h"
#include <Kokkos_Core.hpp>
#include <pcms/assert.h>
#include <Omega_h_for.hpp>
#include "pcms/arrays.h"
#include "pcms/array_mask.h"
#include "pcms/point_search.h"
#include <redev_variant_tools.h>
#include "pcms/transfer_field.h"
#include "pcms/memory_spaces.h"
#include "pcms/profile.h"

// FIXME add executtion spaces (don't use kokkos exe spaces directly)

namespace pcms
{

// TODO different types dependent on active OmegaHBackend
struct OmegaHMemorySpace
{
  using type = typename Kokkos::DefaultExecutionSpace::memory_space;
};

namespace detail
{
template <typename T>
struct memory_space_selector<Omega_h::Read<T>, void>
{
  using type = typename OmegaHMemorySpace::type;
};
template <typename T>
struct memory_space_selector<Omega_h::Write<T>, void>
{
  using type = typename OmegaHMemorySpace::type;
};
template <typename T>
struct memory_space_selector<Omega_h::HostRead<T>, void>
{
  using type = typename pcms::HostMemorySpace;
};
template <typename T>
struct memory_space_selector<Omega_h::HostWrite<T>, void>
{
  using type = typename pcms::HostMemorySpace;
};
template <typename T, int dim = 1>
Omega_h::Read<T> filter_array(Omega_h::Read<T> array,
                              const Omega_h::Read<LO>& mask, LO size)
{
  WDMCPL_FUNCTION_TIMER;
  static_assert(dim > 0, "array dimension must be >0");
  Omega_h::Write<T> filtered_field(size * dim);
  WDMCPL_ALWAYS_ASSERT(array.size() == mask.size() * dim);
  WDMCPL_ALWAYS_ASSERT(filtered_field.size() <= array.size());
  Omega_h::parallel_for(
    mask.size(), OMEGA_H_LAMBDA(LO i) {
      if (mask[i]) {
        const auto idx = mask[i] - 1;
        for (int j = 0; j < dim; ++j) {
          filtered_field[idx * dim + j] = array[i * dim + j];
        }
      }
    });
  return filtered_field;
}
struct GetRankOmegaH
{
  GetRankOmegaH(int i, Omega_h::HostRead<Omega_h::I8> dims,
                Omega_h::HostRead<Omega_h::ClassId> ids)
    : i_(i), ids_(ids), dims_(dims)
  {
    WDMCPL_FUNCTION_TIMER;
  }
  auto operator()(const redev::ClassPtn& ptn) const
  {
    WDMCPL_FUNCTION_TIMER;
    const auto ent = redev::ClassPtn::ModelEnt({dims_[i_], ids_[i_]});
    return ptn.GetRank(ent);
  }
  auto operator()(const redev::RCBPtn& /*unused*/)
  {
    WDMCPL_FUNCTION_TIMER;
    std::cerr << "RCB partition not handled yet\n";
    std::terminate();
    return 0;
  }
  int i_;
  Omega_h::HostRead<Omega_h::ClassId> ids_;
  Omega_h::HostRead<Omega_h::I8> dims_;
};
} // namespace detail

template <typename T,
          typename CoordinateElementType =
            Real> // CoordinateElement<Cartesian, Real>>
class OmegaHField
{
public:
  using memory_space = OmegaHMemorySpace::type;
  using value_type = T;
  using coordinate_element_type = CoordinateElementType;

  OmegaHField(std::string name, Omega_h::Mesh& mesh,
              std::string global_id_name = "", int search_nx = 10,
              int search_ny = 10)
    : name_(std::move(name)),
      mesh_(mesh),
      search_{mesh, search_nx, search_ny},
      size_(mesh.nents(0)),
      global_id_name_(std::move(global_id_name))
  {
    WDMCPL_FUNCTION_TIMER;
  }
  OmegaHField(std::string name, Omega_h::Mesh& mesh,
              Omega_h::Read<Omega_h::I8> mask, std::string global_id_name = "",
              int search_nx = 10, int search_ny = 10)
    : name_(std::move(name)),
      mesh_(mesh),
      search_{mesh, search_nx, search_ny},
      global_id_name_(std::move(global_id_name))
  {
    WDMCPL_FUNCTION_TIMER;
    if (mask.exists()) {

      using ExecutionSpace = typename memory_space::execution_space;
      auto policy = Kokkos::RangePolicy<ExecutionSpace>(0, mask.size());
      // we use a parallel scan to construct the mask mapping so that filtering
      // can happen in parallel. This method gives us the index to fill into the
      // filtered array
      WDMCPL_ALWAYS_ASSERT(mesh.nents(0) == mask.size());
      Omega_h::Write<LO> index_mask(mask.size());
      auto index_mask_view = make_array_view(index_mask);
      auto mask_view = make_const_array_view(mask);
      Kokkos::parallel_scan(
        policy, detail::ComputeMaskAV{index_mask_view, mask_view}, size_);
      Kokkos::parallel_for(policy, detail::ScaleAV{index_mask_view, mask_view});
      mask_ = index_mask;
    } else {
      size_ = mesh.nents(0);
    }
  }

  [[nodiscard]] const std::string& GetName() const noexcept { return name_; }
  [[nodiscard]] Omega_h::Mesh& GetMesh() const noexcept { return mesh_; }
  [[nodiscard]] const Omega_h::Read<LO>& GetMask() const noexcept
  {
    return mask_;
  };
  [[nodiscard]] bool HasMask() const noexcept { return mask_.exists(); };
  [[nodiscard]] LO Size() const noexcept { return size_; }
  // pass through to search function
  auto Search(Kokkos::View<Real* [2]> points) const {
    WDMCPL_FUNCTION_TIMER;
    return search_(points); }

  [[nodiscard]] Omega_h::Read<Omega_h::ClassId> GetClassIDs() const
  {
    WDMCPL_FUNCTION_TIMER;
    if (HasMask())
      return detail::filter_array(
        mesh_.get_array<Omega_h::ClassId>(0, "class_id"), GetMask(), Size());
    return mesh_.get_array<Omega_h::ClassId>(0, "class_id");
  }
  [[nodiscard]] Omega_h::Read<Omega_h::I8> GetClassDims() const
  {
    WDMCPL_FUNCTION_TIMER;
    if (HasMask())
      return detail::filter_array(mesh_.get_array<Omega_h::I8>(0, "class_dim"),
                                  GetMask(), Size());
    return mesh_.get_array<Omega_h::I8>(0, "class_dim");
  }
  [[nodiscard]] Omega_h::Read<Omega_h::GO> GetGids() const
  {
    WDMCPL_FUNCTION_TIMER;
    Omega_h::Read<Omega_h::GO> gid_array;
    if (global_id_name_.empty()) {
      gid_array = mesh_.globals(0);
    } else {
      auto tag = mesh_.get_tagbase(0, global_id_name_);
      if (Omega_h::is<GO>(tag)) {
        gid_array = mesh_.get_array<Omega_h::GO>(0, global_id_name_);
      } else if (Omega_h::is<LO>(tag)) {
        auto array = mesh_.get_array<Omega_h::LO>(0, global_id_name_);
        Omega_h::Write<Omega_h::GO> globals(array.size());
        Omega_h::parallel_for(
          array.size(), OMEGA_H_LAMBDA(int i) { globals[i] = array[i]; });
        gid_array = Omega_h::Read(globals);
      } else {
        std::cerr << "Weird tag type for global arrays.\n";
        std::abort();
      }
    }
    if (HasMask()) {
      return detail::filter_array(gid_array, GetMask(), Size());
    }
    return gid_array;
  }
private:
  std::string name_;
  Omega_h::Mesh& mesh_;
  GridPointSearch search_;
  // bitmask array that specifies a filter on the field
  Omega_h::Read<LO> mask_;
  LO size_;
  std::string global_id_name_;
};

using InternalCoordinateElement = Real;
// internal field can only be one of the types supported by Omega_h
// The coordinate element for all internal fields is the same since
// all internal fields are on the same mesh
using InternalField =
  std::variant<OmegaHField<Omega_h::I8, InternalCoordinateElement>,
               OmegaHField<Omega_h::I32, InternalCoordinateElement>,
               OmegaHField<Omega_h::I64, InternalCoordinateElement>,
               OmegaHField<Omega_h::Real, InternalCoordinateElement>>;

template <typename T, typename CoordinateElementType>
auto get_nodal_data(const OmegaHField<T, CoordinateElementType>& field)
  -> Omega_h::Read<T>
{
  WDMCPL_FUNCTION_TIMER;
  auto full_field = field.GetMesh().template get_array<T>(0, field.GetName());
  if (field.HasMask()) {
    return detail::filter_array<T>(full_field, field.GetMask(), field.Size());
  }
  return full_field;
}

// TODO since Omega_h owns coordinate data, we could potentially
// return a view of the data without lifetime issues.
template <typename T, typename CoordinateElementType>
auto get_nodal_coordinates(const OmegaHField<T, CoordinateElementType>& field)
{
  WDMCPL_FUNCTION_TIMER;
  static constexpr auto coordinate_dimension = 2;
  if constexpr (detail::HasCoordinateSystem<CoordinateElementType>::value) {
    const auto coords = field.GetMesh().coords();
    return MDArray<CoordinateElementType>{};
    // FIXME implement copy to
    throw;
  } else {
    auto coords = Omega_h::Reals{field.GetMesh().coords()};
    if (field.HasMask()) {
      return detail::filter_array<typename decltype(coords)::value_type,
                                  coordinate_dimension>(coords, field.GetMask(),
                                                        field.Size());
    }
    return coords;
  }
  // should never be here. Quash warning
  return Omega_h::Reals{};
}

/**
 * Sets the data on the entire mesh
 */
template <typename T, typename CoordinateElementType, typename U>
auto set_nodal_data(const OmegaHField<T, CoordinateElementType>& field,
                    ScalarArrayView<const U, OmegaHMemorySpace::type> data)
  -> void
{
  WDMCPL_FUNCTION_TIMER;
  static_assert(std::is_convertible_v<T, U>,
                "must be able to convert nodal data into the field types data");
  auto& mesh = field.GetMesh();
  const auto has_tag = mesh.has_tag(0, field.GetName());
  if (field.HasMask()) {
    auto& mask = field.GetMask();
    WDMCPL_ALWAYS_ASSERT(mask.size() == mesh.nents(0));
    Omega_h::Write<T> array(mask.size());
    if (has_tag) {
      auto original_data = mesh.template get_array<T>(0, field.GetName());
      WDMCPL_ALWAYS_ASSERT(original_data.size() == mask.size());
      Omega_h::parallel_for(
        mask.size(), OMEGA_H_LAMBDA(size_t i) {
          array[i] = mask[i] ? data(mask[i] - 1) : original_data[i];
        });
      mesh.set_tag(0, field.GetName(), Omega_h::Read<T>(array));
    } else {
      Omega_h::parallel_for(
        mask.size(), OMEGA_H_LAMBDA(size_t i) {
          array[i] = mask[i] ? data(mask[i] - 1) : 0;
        });
      mesh.add_tag(0, field.GetName(), 1, Omega_h::Read<T>(array));
    }
  } else {
    WDMCPL_ALWAYS_ASSERT(static_cast<LO>(data.size()) == mesh.nents(0));
    Omega_h::Write<T> array(data.size());
    Omega_h::parallel_for(
      data.size(), OMEGA_H_LAMBDA(size_t i) { array[i] = data(i); });
    if (has_tag) {
      mesh.set_tag(0, field.GetName(), Omega_h::Read<T>(array));
    } else {
      mesh.add_tag(0, field.GetName(), 1, Omega_h::Read<T>(array));
    }
  }
  WDMCPL_ALWAYS_ASSERT(mesh.has_tag(0, field.GetName()));
}

// TODO abstract out repeat parts of lagrange/nearest neighbor evaluation
template <typename T, typename CoordinateElementType>
auto evaluate(
  const OmegaHField<T, CoordinateElementType>& field, Lagrange<1> /* method */,
  ScalarArrayView<const CoordinateElementType, OmegaHMemorySpace::type>
    coordinates) -> Omega_h::Read<T>
{
  WDMCPL_FUNCTION_TIMER;
  Omega_h::Write<T> values(coordinates.size() / 2);
  auto tris2verts = field.GetMesh().ask_elem_verts();
  auto field_values = field.GetMesh().template get_array<T>(0, field.GetName());

  Kokkos::View<Real* [2]> coords("coords", coordinates.size() / 2);
  Kokkos::parallel_for(
    coordinates.size() / 2, KOKKOS_LAMBDA(LO i) {
      coords(i, 0) = coordinates(2 * i);
      coords(i, 1) = coordinates(2 * i + 1);
    });
  auto results = field.Search(coords);

  Kokkos::parallel_for(
    results.size(), KOKKOS_LAMBDA(LO i) {
      auto [elem_idx, coord] = results(i);
      // TODO deal with case for elem_idx < 0 (point outside of mesh)
      KOKKOS_ASSERT(elem_idx >= 0);
      const auto elem_tri2verts =
        Omega_h::gather_verts<3>(tris2verts, elem_idx);
      Real val = 0;
      for (int j = 0; j < 3; ++j) {
        val += field_values[elem_tri2verts[j]] * coord[j];
      }
      if constexpr (std::is_integral_v<T>) {
        val = std::round(val);
      }
      values[i] = val;
    });

  return values;
}

template <typename T, typename CoordinateElementType>
auto evaluate(
  const OmegaHField<T, CoordinateElementType>& field,
  NearestNeighbor /* method */,
  ScalarArrayView<const CoordinateElementType, OmegaHMemorySpace::type>
    coordinates) -> Omega_h::Read<T>
{
  WDMCPL_FUNCTION_TIMER;
  Omega_h::Write<T> values(coordinates.size() / 2);
  auto tris2verts = field.GetMesh().ask_elem_verts();
  auto field_values = field.GetMesh().template get_array<T>(0, field.GetName());
  // TODO reuse coordinates_data if possible
  Kokkos::View<Real* [2]> coords("coords", coordinates.size() / 2);
  Kokkos::parallel_for(
    coordinates.size() / 2, KOKKOS_LAMBDA(LO i) {
      coords(i, 0) = coordinates(2 * i);
      coords(i, 1) = coordinates(2 * i + 1);
    });
  auto results = field.Search(coords);

  Kokkos::parallel_for(
    results.size(), KOKKOS_LAMBDA(LO i) {
      auto [elem_idx, coord] = results(i);
      // TODO deal with case for elem_idx < 0 (point outside of mesh)
      KOKKOS_ASSERT(elem_idx >= 0);
      const auto elem_tri2verts =
        Omega_h::gather_verts<3>(tris2verts, elem_idx);
      // value is closest to point has the largest coordinate
      int vert = 0;
      auto max_val = coord[0];
      for (int j = 1; j <= 2; ++j) {
        auto next_val = coord[j];
        if (next_val > max_val) {
          max_val = next_val;
          vert = j;
        }
      }
      values[i] = field_values[elem_tri2verts[vert]];
    });
  return values;
}

template <typename T, typename Method, typename CoordinateElementType>
auto evaluate(
  const OmegaHField<T, CoordinateElementType>& field, Method&& m,
  ScalarArrayView<const CoordinateElementType, HostMemorySpace> coordinates)
  -> std::enable_if_t<
    !std::is_same_v<typename OmegaHMemorySpace::type, HostMemorySpace>,
    Omega_h::HostRead<T>>

{
  WDMCPL_FUNCTION_TIMER;
  auto coords_view =
    Kokkos::View<const CoordinateElementType, Kokkos::HostSpace,
                 Kokkos::MemoryTraits<Kokkos::Unmanaged>>(&coordinates[0],
                                                          coordinates.size());
  using exe_space = typename OmegaHMemorySpace::type::execution_space;
  auto coordinates_d =
    Kokkos::create_mirror_view_and_copy(exe_space(), coords_view);
  return Omega_h::HostRead<T>(evaluate(field, std::forward<Method>(m),
                                       make_const_array_view(coordinates_d)));
}

} // namespace pcms
namespace Omega_h
{
template <typename T>
auto make_array_view(const Omega_h::Read<T>& array)
  -> pcms::ScalarArrayView<const T, typename pcms::OmegaHMemorySpace::type>
{
  WDMCPL_FUNCTION_TIMER;
  pcms::ScalarArrayView<const T, typename pcms::OmegaHMemorySpace::type>
    view(array.data(), array.size());
  return view;
}
} // namespace Omega_h

namespace pcms
{

template <typename T, typename CoordinateElementType = Real>
class OmegaHFieldAdapter
{
public:
  using memory_space = OmegaHMemorySpace::type;
  using value_type = T;
  using coordinate_element_type = CoordinateElementType;
  OmegaHFieldAdapter(std::string name, Omega_h::Mesh& mesh,
                     std::string global_id_name = "", int search_nx = 10,
                     int search_ny = 10)
    : field_{std::move(name), mesh, std::move(global_id_name), search_nx,
             search_ny}
  {
    WDMCPL_FUNCTION_TIMER;
  }

  OmegaHFieldAdapter(std::string name, Omega_h::Mesh& mesh,
                     Omega_h::Read<Omega_h::I8> mask,
                     std::string global_id_name = "", int search_nx = 10,
                     int search_ny = 10)
    : field_{std::move(name),           mesh,      mask,
             std::move(global_id_name), search_nx, search_ny}
  {
    WDMCPL_FUNCTION_TIMER;
  }
  [[nodiscard]] const std::string& GetName() const noexcept
  {
    return field_.GetName();
  }
  // REQUIRED
  int Serialize(ScalarArrayView<T, pcms::HostMemorySpace> buffer,
                ScalarArrayView<const pcms::LO, pcms::HostMemorySpace>
                  permutation) const
  {
    WDMCPL_FUNCTION_TIMER;
    // host copy of filtered field data array
    const auto array_h = Omega_h::HostRead<T>(get_nodal_data(field_));
    if (buffer.size() > 0) {
      for (LO i = 0; i < array_h.size(); i++) {
        buffer[i] = array_h[permutation[i]];
      }
    }
    return array_h.size();
  }
  // REQUIRED
  void Deserialize(ScalarArrayView<const T, pcms::HostMemorySpace> buffer,
                   ScalarArrayView<const pcms::LO, pcms::HostMemorySpace>
                     permutation) const
  {
    WDMCPL_FUNCTION_TIMER;
    REDEV_ALWAYS_ASSERT(buffer.size() == permutation.size());
    Omega_h::HostWrite<T> sorted_buffer(buffer.size());
    for (size_t i = 0; i < buffer.size(); ++i) {
      sorted_buffer[permutation[i]] = buffer[i];
    }
    const auto sorted_buffer_d = Omega_h::Read<T>(sorted_buffer);
    set_nodal_data(field_, make_array_view(sorted_buffer_d));
  }

  [[nodiscard]] std::vector<GO> GetGids() const
  {
    WDMCPL_FUNCTION_TIMER;
    auto gids = field_.GetGids();
    if (gids.size() > 0) {
      auto gids_h = Omega_h::HostRead<GO>(gids);
      return {&gids_h[0], &(gids_h[gids_h.size() - 1]) + 1};
    }
    return {};
  }
  // REQUIRED
  [[nodiscard]] ReversePartitionMap GetReversePartitionMap(
    const redev::Partition& partition) const
  {
    WDMCPL_FUNCTION_TIMER;
    auto classIds_h = Omega_h::HostRead<Omega_h::ClassId>(field_.GetClassIDs());
    auto classDims_h = Omega_h::HostRead<Omega_h::I8>(field_.GetClassDims());

    // local_index number of vertices going to each destination process by
    // calling getRank - degree array
    pcms::ReversePartitionMap reverse_partition;
    pcms::LO local_index = 0;
    for (auto i = 0; i < classIds_h.size(); i++) {
      auto dr = std::visit(detail::GetRankOmegaH{i, classDims_h, classIds_h},
                           partition);
      reverse_partition[dr].emplace_back(local_index++);
    }
    return reverse_partition;
  }
  // NOT REQUIRED PART OF FieldAdapter interface
  [[nodiscard]] OmegaHField<T, CoordinateElementType>& GetField() noexcept
  {
    return field_;
  }
  // NOT REQUIRED PART OF FieldAdapter interface
  [[nodiscard]] const OmegaHField<T, CoordinateElementType>& GetField()
    const noexcept
  {
    return field_;
  }

private:
  OmegaHField<T, CoordinateElementType> field_;
};
template <typename FieldAdapter>
void ConvertFieldAdapterToOmegaH(const FieldAdapter& adapter,
                                 InternalField internal,
                                 FieldTransferMethod ftm,
                                 FieldEvaluationMethod fem)
{
  WDMCPL_FUNCTION_TIMER;
  std::visit(
    [&](auto&& internal_field) {
      transfer_field(adapter, internal_field, ftm, fem);
    },
    internal);
}

template <typename FieldAdapter>
void ConvertOmegaHToFieldAdapter(const InternalField& internal,
                                 FieldAdapter& adapter, FieldTransferMethod ftm,
                                 FieldEvaluationMethod fem)
{
  WDMCPL_FUNCTION_TIMER;
  std::visit(
    [&](auto&& internal_field) {
      transfer_field(internal_field, adapter, ftm, fem);
    },
    internal);
}
// Specializations for the Omega_h field adapter class since get/set are
// implemented on the OmegaHFieldClass which is owned by the field adapter
template <typename T, typename C>
void ConvertFieldAdapterToOmegaH(const OmegaHFieldAdapter<T, C>& adapter,
                                 InternalField internal,
                                 FieldTransferMethod ftm,
                                 FieldEvaluationMethod fem)
{
  WDMCPL_FUNCTION_TIMER;
  std::visit(
    [&](auto&& internal_field) {
      transfer_field(adapter.GetField(), internal_field, ftm, fem);
    },
    internal);
}
template <typename T, typename C>
void ConvertOmegaHToFieldAdapter(const InternalField& internal,
                                 OmegaHFieldAdapter<T, C>& adapter,
                                 FieldTransferMethod ftm,
                                 FieldEvaluationMethod fem)
{
  WDMCPL_FUNCTION_TIMER;
  std::visit(
    [&](auto&& internal_field) {
      transfer_field(internal_field, adapter.GetField(), ftm, fem);
    },
    internal);
}

} // namespace pcms

#endif // WDM_COUPLING_OMEGA_H_FIELD_H
