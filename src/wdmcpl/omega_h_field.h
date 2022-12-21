#ifndef WDM_COUPLING_OMEGA_H_FIELD_H
#define WDM_COUPLING_OMEGA_H_FIELD_H
#include "wdmcpl/types.h"
#include "wdmcpl/external/span.h"
#include <Omega_h_mesh.hpp>
#include "wdmcpl/field.h"
#include "wdmcpl/coordinate_systems.h"
#include <Kokkos_core.hpp>
#include <wdmcpl/assert.h>
#include <Omega_h_for.hpp>
#include "wdmcpl/arrays.h"
#include "wdmcpl/point_search.h"
#include <redev_variant_tools.h>
#include "wdmcpl/transfer_field.h"

// FIXME add executtion spaces (don't use kokkos exe spaces directly)

namespace wdmcpl
{

// TODO different types dependent on active OmegaHBackend
struct OmegaHMemorySpace
{
  using type = typename Kokkos::DefaultExecutionSpace::memory_space;
};


namespace detail
{
template <typename T, int dim = 1>
Omega_h::Read<T> filter_array(Omega_h::Read<T> array, Omega_h::Read<LO> mask,
                              LO size)
{
  static_assert(dim > 0, "array dimension must be >0");
  Omega_h::Write<T> filtered_field(size * dim);
  WDMCPL_ALWAYS_ASSERT(array.size() == mask.size());
  WDMCPL_ALWAYS_ASSERT(filtered_field.size() <= array.size());
  Omega_h::parallel_for(
    mask.size(), OMEGA_H_LAMBDA(LO i) {
      if (mask[i]) {
        const auto idx = mask[i] - 1;
        for (int j = 0; j < dim; ++j) {
          filtered_field[idx * dim + j] = array[i];
        }
      }
    });
  return filtered_field;
}
} // namespace detail

template <typename T,
          typename CoordinateElementType =
            Real> // CoordinateElement<Cartesian, Real>>
class OmegaHField
{
public:
  using value_type = T;
  using coordinate_element_type = CoordinateElementType;

  OmegaHField(std::string name, Omega_h::Mesh& mesh, int search_nx = 10,
              int search_ny = 10)
    : name_(std::move(name)),
      mesh_(mesh),
      search_{mesh, search_nx, search_ny},
      size_(mesh.nents(0))
  {
  }
  OmegaHField(std::string name, Omega_h::Mesh& mesh,
              Omega_h::Read<Omega_h::I8> mask, int search_nx = 10,
              int search_ny = 10)
    : name_(std::move(name)), mesh_(mesh), search_{mesh, search_nx, search_ny}
  {
    if (mask.exists()) {

      // we use a parallel scan to construct the mask mapping so that filtering
      // can happen in parallel. This method gives us the index to fill into the
      // filtered array
      WDMCPL_ALWAYS_ASSERT(mesh.nents(0) == mask.size());
      Omega_h::Write<LO> index_mask(mask.size());
      Kokkos::parallel_scan(
        mask.size(),
        KOKKOS_LAMBDA(LO i, LO & update, bool final) {
          update += (mask[i] > 0);
          if (final) {
            index_mask[i] = update;
          }
        },
        size_);
      // set index mask to 0 anywhere that the original mask is 0
      Kokkos::parallel_for(
        mask.size(), KOKKOS_LAMBDA(LO i) { index_mask[i] *= mask[i]; });

      mask_ = index_mask;
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
  template <typename... Ts>
  auto Search(Ts... args) const
  {
    return search_(std::forward<Ts>(args)...);
  }

  [[nodiscard]] Omega_h::Read<Omega_h::ClassId> GetClassIDs() const
  {
    if (HasMask())
      return detail::filter_array(
        mesh_.get_array<Omega_h::ClassId>(0, "class_id"), GetMask(), Size());
    return mesh_.get_array<Omega_h::ClassId>(0, "class_id");
  }
  [[nodiscard]] Omega_h::Read<Omega_h::I8> GetClassDims() const
  {
    if (HasMask())
      return detail::filter_array(mesh_.get_array<Omega_h::I8>(0, "class_dim"),
                                  GetMask(), Size());
    return mesh_.get_array<Omega_h::I8>(0, "class_dim");
  }

private:
  std::string name_;
  Omega_h::Mesh& mesh_;
  // FIXME GridPointSearch take mdspan
  GridPointSearch search_;
  // bitmask array that specifies a filter on the field
  Omega_h::Read<LO> mask_;
  LO size_;
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
  auto full_field = field.GetMesh().template get_array<T>(0, field.GetName());
  if (field.HasMask()) {
    return detail::filter_array<T>(full_field, field.GetMask(), field.Size());
  }
  return full_field;
}

// TODO Remove as customization point
template <typename T, typename CoordinateElementType>
auto get_nodal_gids(const OmegaHField<T, CoordinateElementType>& field)
  -> Omega_h::Read<Omega_h::GO>
{
  auto full_gids = field.GetMesh().globals(0);
  if (field.HasMask()) {
    return detail::filter_array(full_gids, field.GetMask(), field.Size());
  }
  return full_gids;
}

// TODO since Omega_h owns coordinate data, we could potentially
// return a view of the data without lifetime issues.
template <typename T, typename CoordinateElementType>
auto get_nodal_coordinates(const OmegaHField<T, CoordinateElementType>& field)
{
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
}

/**
 * Sets the data on the entire mesh
 */
template <typename T, typename CoordinateElementType, typename U>
auto set_nodal_data(const OmegaHField<T, CoordinateElementType>& field,
                    ScalarArrayView<const U, OmegaHMemorySpace::type> data)
  -> void
{
  auto& mesh = field.GetMesh();
  const auto has_tag = mesh.has_tag(0, field.GetName());
  if (field.HasMask()) {
    auto& mask = field.GetMask();
    WDMCPL_ALWAYS_ASSERT(mask.size() == mesh.nents(0));
    Omega_h::Write<U> array(mask.size());
    if (has_tag) {
      auto original_data = mesh.template get_array<U>(0, field.GetName());
      WDMCPL_ALWAYS_ASSERT(original_data.size() == mask.size());
      Omega_h::parallel_for(
        mask.size(), OMEGA_H_LAMBDA(size_t i) {
          array[i] = mask[i] ? data(mask[i] - 1) : original_data[i];
        });
      mesh.set_tag(0, field.GetName(), Omega_h::Read(array));
    } else {
      Omega_h::parallel_for(
        mask.size(), OMEGA_H_LAMBDA(size_t i) {
          array[i] = mask[i] ? data(mask[i] - 1) : 0;
        });
      mesh.add_tag(0, field.GetName(), 1, Omega_h::Read(array));
    }
  } else {
    WDMCPL_ALWAYS_ASSERT(static_cast<LO>(data.size()) == mesh.nents(0));
    Omega_h::Write<U> array(data.size());
    Omega_h::parallel_for(
      data.size(), OMEGA_H_LAMBDA(size_t i) { array[i] = data(i); });
    if (has_tag) {
      mesh.set_tag(0, field.GetName(), Omega_h::Read(array));
    } else {
      mesh.add_tag(0, field.GetName(), 1, Omega_h::Read(array));
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
  Omega_h::Write<T> values(coordinates.size() / 2);
  auto tris2verts = field.GetMesh().ask_elem_verts();
  auto field_values = field.GetMesh().template get_array<T>(0, field.GetName());
  // FIXME field search operation needs to be on GPU
  // FIXME field search takes mdspan
  // FIXME coordinate array view as nx2 span
  Omega_h::parallel_for(
    coordinates.size() / 2, OMEGA_H_LAMBDA(LO i) {
      auto [elem_idx, coord] = field.Search(
        Omega_h::Vector<2>{coordinates(2 * i), coordinates(2 * i + 1)});
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
  const OmegaHField<T, CoordinateElementType>& field, NearestNeighbor /* method */,
  ScalarArrayView<const CoordinateElementType, OmegaHMemorySpace::type>
    coordinates) -> Omega_h::Read<T>
{
  Omega_h::Write<T> values(coordinates.size() / 2);
  auto tris2verts = field.GetMesh().ask_elem_verts();
  auto field_values = field.GetMesh().template get_array<T>(0, field.GetName());
  // FIXME field search operation needs to be on GPU
  // FIXME field search takes mdspan
  // FIXME coordinate array view as nx2 span
  Omega_h::parallel_for(
    coordinates.size() / 2, OMEGA_H_LAMBDA(LO i) {
      auto [elem_idx, coord] = field.Search(
        Omega_h::Vector<2>{coordinates(2 * i), coordinates(2 * i + 1)});
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

} // namespace wdmcpl
namespace Omega_h
{
template <typename T>
auto make_array_view(const Omega_h::Read<T>& array)
  -> wdmcpl::ScalarArrayView<const T, typename wdmcpl::OmegaHMemorySpace::type>
{
  wdmcpl::ScalarArrayView<const T, typename wdmcpl::OmegaHMemorySpace::type>
    view(array.data(), array.size());
  return view;
}
} // namespace Omega_h

namespace wdmcpl
{

template <typename T, typename CoordinateElementType = Real>
class OmegaHFieldAdapter
{
public:
  using memory_space = OmegaHMemorySpace::type;
  using value_type = T;
  using coordinate_element_type = CoordinateElementType;
  OmegaHFieldAdapter(std::string name, Omega_h::Mesh& mesh, int search_nx = 10,
                  int search_ny = 10)
    : field_{std::move(name), mesh, search_nx, search_ny}
  {
  }

  OmegaHFieldAdapter(std::string name, Omega_h::Mesh& mesh,
                  Omega_h::Read<Omega_h::I8> mask, int search_nx = 10,
                  int search_ny = 10)
    : field_{std::move(name), mesh, mask, search_nx, search_ny}
  {
  }
  [[nodiscard]] const std::string& GetName() const noexcept
  {
    return field_.GetName();
  }
  // REQUIRED
  int Serialize(
    ScalarArrayView<T, memory_space> buffer,
    ScalarArrayView<const wdmcpl::LO, memory_space> permutation) const
  {
    // host copy of filtered field data array
    const auto array_h = Omega_h::HostRead(get_nodal_data(field_));
    if (buffer.size() > 0) {
      for (LO i = 0, j = 0; i < array_h.size(); i++) {
        // FIXME j==i, can't we just get rid of j, or am i missing something
        buffer[permutation[j++]] = array_h[i];
      }
    }
    return array_h.size();
  }
  // REQUIRED
  void Deserialize(
    ScalarArrayView<T, memory_space> buffer,
    ScalarArrayView<const wdmcpl::LO, memory_space> permutation) const
  {
    REDEV_ALWAYS_ASSERT(buffer.size() == permutation.size());
    Omega_h::Write<T> sorted_buffer(buffer.size());
    for (size_t i = 0; i < buffer.size(); ++i) {
      sorted_buffer[i] = buffer[permutation[i]];
    }
    set_nodal_data(field_, make_array_view(Omega_h::Read(sorted_buffer)));
  }

  // REQUIRED
  [[nodiscard]] std::vector<GO> GetGids() const
  {
    auto gids = get_nodal_gids(field_);
    if (gids.size() > 0) {
      auto gids_h = Omega_h::HostRead(gids);
      return {&gids_h[0], &(gids_h[gids_h.size() - 1]) + 1};
    }
    return {};
  }
  // REQUIRED
  [[nodiscard]] ReversePartitionMap GetReversePartitionMap(
    const redev::Partition& partition) const
  {
    auto classIds_h = Omega_h::HostRead(field_.GetClassIDs());
    auto classDims_h = Omega_h::HostRead(field_.GetClassDims());

    // local_index number of vertices going to each destination process by
    // calling getRank - degree array
    wdmcpl::ReversePartitionMap reverse_partition;
    wdmcpl::LO local_index = 0;
    for (auto i = 0; i < classIds_h.size(); i++) {
      auto dr = std::visit(
        redev::overloaded{
          [&classDims_h, &classIds_h, &i](const redev::ClassPtn& ptn) {
            const auto ent =
              redev::ClassPtn::ModelEnt({classDims_h[i], classIds_h[i]});
            return ptn.GetRank(ent);
          },
          [](const redev::RCBPtn& /*ptn*/) {
            std::cerr << "RCB partition not handled yet\n";
            std::terminate();
            return 0;
          }},
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
  std::visit(
    [&](auto&& internal_field) {
      transfer_field(internal_field, adapter.GetField(), ftm, fem);
    },
    internal);
}

} // namespace wdmcpl

#endif // WDM_COUPLING_OMEGA_H_FIELD_H
