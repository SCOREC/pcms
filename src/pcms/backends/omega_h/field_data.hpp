#ifndef PCMS_BACKENDS_OMEGA_H_FIELD_DATA_HPP
#define PCMS_BACKENDS_OMEGA_H_FIELD_DATA_HPP

#include "pcms/field/field_data.hpp"
#include "pcms/field/field_layout.hpp"
#include "pcms/backends/omega_h/omega_h_geometry.hpp"
#include "pcms/utility/arrays.hpp"
#include "pcms/utility/memory_spaces.hpp"
#include "pcms/utility/omega_h_array_utils.hpp"
#include <Omega_h_array.hpp>
#include <Omega_h_mesh.hpp>
#include <Kokkos_Core.hpp>
#include <atomic>
#include <memory>
#include <string>

namespace pcms::omega_h
{

namespace detail
{
// A tag array (one value per component per entity, entity-major) holding the
// owned values at their entities and zero at the ghosts. Local to this rank.
template <typename T, typename TagValue>
Omega_h::Read<TagValue> OwnedValuesAsTag(
  Rank2View<const T, DeviceMemorySpace> owned,
  Kokkos::View<const LO*, DeviceMemorySpace> owned_indices, LO num_entities,
  int num_components)
{
  Omega_h::Write<TagValue> tag(num_entities * num_components, TagValue(0));
  Kokkos::parallel_for(
    "pcms::OwnedValuesAsTag",
    Kokkos::RangePolicy<DeviceMemorySpace::execution_space>(
      0, static_cast<LO>(owned.extent(0))),
    KOKKOS_LAMBDA(LO o) {
      const LO i = owned_indices(o);
      for (int c = 0; c < num_components; ++c) {
        tag[i * num_components + c] = static_cast<TagValue>(owned(o, c));
      }
    });
  return tag;
}

// A name for a tag an OmegaHFieldData creates, unique within the process.
inline std::string NextOmegaHFieldDataTagName()
{
  static std::atomic<long> counter{0};
  return "pcms_field_" + std::to_string(counter++);
}
} // namespace detail

// OmegaHFieldData<T> is the FieldData<T> of an Omega_h layout. It stores the
// owned values in a [dof_holder][component] Kokkos View and also keeps the
// values of every DOF holder, owned and ghost, in an Omega_h tag, for code that
// reads ghosts (e.g. the Omega_h Lagrange evaluator).
//
// It creates a tag of its own on the mesh, named by
// detail::NextOmegaHFieldDataTagName(), and removes it on destruction, so the
// mesh must outlive it.
// SetValues replaces the tag with the owned values (ghosts zero, local to this
// rank); Synchronize (collective) fills the ghosts from their owners. Ghost
// values are correct only after Synchronize.
template <typename T>
class OmegaHFieldData
  : public FieldData<T>
  , public DOFAccess<T>
{
public:
  // Value type of the Omega_h tag.
  using tag_value_type = utility::OmegaHValue_t<T>;

  OmegaHFieldData(std::shared_ptr<const FieldLayout> layout, ValueBasis basis)
    : layout_(std::move(layout)),
      basis_(std::move(basis)),
      device_data_("omega_h_field_data_device",
                   static_cast<size_t>(layout_->GetNumDofHolder()),
                   static_cast<size_t>(layout_->GetNumComponents()))
  {
    const auto* oh =
      dynamic_cast<const OmegaHGeometry*>(layout_->GetGeometry().get());
    if (oh == nullptr) {
      throw pcms_error("OmegaHFieldData: layout is not on an Omega_h geometry");
    }
    tag_dim_ = layout_->GetDOFHolderEntityDim();
    if (tag_dim_ < 0) {
      throw pcms_error("OmegaHFieldData: layout DOF holders are not one per "
                       "mesh entity of a single dimension");
    }
    mesh_ = &oh->GetMesh();
    tag_name_ = detail::NextOmegaHFieldDataTagName();
    const int nc = layout_->GetNumComponents();
    mesh_->add_tag<tag_value_type>(
      tag_dim_, tag_name_, nc,
      Omega_h::Read<tag_value_type>(mesh_->nents(tag_dim_) * nc,
                                    tag_value_type(0)));
  }

  // The tag holds the name of this object, so it is neither copied nor moved.
  OmegaHFieldData(const OmegaHFieldData&) = delete;
  OmegaHFieldData& operator=(const OmegaHFieldData&) = delete;
  OmegaHFieldData(OmegaHFieldData&&) = delete;
  OmegaHFieldData& operator=(OmegaHFieldData&&) = delete;

  ~OmegaHFieldData() noexcept override
  {
    if (mesh_->has_tag(tag_dim_, tag_name_)) {
      mesh_->remove_tag(tag_dim_, tag_name_);
    }
  }

  FieldValueType GetValueType() const override
  {
    return pcms::detail::ValueTypeOfRank(basis_.Rank());
  }
  const ValueBasis& GetValueBasis() const override { return basis_; }

  const DOFAccess<T>* GetDOFAccess() const noexcept override { return this; }
  DOFAccess<T>* GetDOFAccess() noexcept override { return this; }

  Rank2View<const T, DeviceMemorySpace> GetValues() const override
  {
    return pcms::MakeConstRank2View(device_data_);
  }

  void SetValues(Rank2View<const T, DeviceMemorySpace> values) override
  {
    pcms::detail::CopyDeviceRank2ViewToDeviceView(device_data_, values);
    mesh_->set_tag(tag_dim_, tag_name_,
                   detail::OwnedValuesAsTag<T, tag_value_type>(
                     GetValues(), layout_->GetOwnedIndices(),
                     mesh_->nents(tag_dim_), layout_->GetNumComponents()));
  }

  // Fills the tag's ghost entries with their owners' values. Collective.
  void Synchronize() override { mesh_->sync_tag(tag_dim_, tag_name_); }

  // Sets the values from one entry per entity and component, of which only
  // the owners' entries are read, and fills the ghosts. Collective.
  void SetFromOwners(Omega_h::Read<tag_value_type> values)
  {
    const int nc = layout_->GetNumComponents();
    const auto synced = mesh_->sync_array(tag_dim_, values, nc);
    mesh_->set_tag(tag_dim_, tag_name_, synced);
    const auto owned = layout_->GetOwnedIndices();
    const auto data = device_data_;
    Kokkos::parallel_for(
      "pcms::OmegaHFieldData::SetFromOwners",
      Kokkos::MDRangePolicy<DeviceMemorySpace::execution_space,
                            Kokkos::Rank<2>>(
        {0, 0},
        {static_cast<int64_t>(data.extent(0)), static_cast<int64_t>(nc)}),
      KOKKOS_LAMBDA(const int64_t o, const int64_t c) {
        const LO i = owned(o) * nc + static_cast<LO>(c);
        data(o, c) = static_cast<T>(synced[i]);
      });
  }

  // Name of the Omega_h tag holding the values of every DOF holder, owned and
  // ghost, on entity dimension layout.GetDOFHolderEntityDim() (one value per
  // component per entity). Its owned entries
  // are always current; its ghost entries are correct only after Synchronize.
  [[nodiscard]] const std::string& GetTagName() const noexcept
  {
    return tag_name_;
  }

private:
  std::shared_ptr<const FieldLayout> layout_;
  ValueBasis basis_;
  Kokkos::View<T**, DeviceMemorySpace> device_data_;
  Omega_h::Mesh* mesh_ = nullptr;
  int tag_dim_ = -1;
  std::string tag_name_;
};

} // namespace pcms::omega_h

#endif // PCMS_BACKENDS_OMEGA_H_FIELD_DATA_HPP
