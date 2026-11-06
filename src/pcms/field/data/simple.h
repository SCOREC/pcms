#ifndef PCMS_SIMPLE_FIELD_DATA_H
#define PCMS_SIMPLE_FIELD_DATA_H

#include "../field_data.h"
#include "../field_layout.h"
#include "pcms/utility/arrays.h"
#include "pcms/utility/assert.h"
#include <Omega_h_array.hpp>
#include <Omega_h_mesh.hpp>
#include "pcms/discretization/discretization/omega_h.hpp"
#include "pcms/utility/memory_spaces.h"
#include "pcms/utility/omega_h_array_utils.h"
#include <Kokkos_Core.hpp>
#include <atomic>
#include <memory>
#include <string>
#include <type_traits>

namespace pcms
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

// A name for a tag a SimpleFieldData creates, unique within the process.
inline std::string NextSimpleFieldDataTagName()
{
  static std::atomic<long> counter{0};
  return "pcms_field_" + std::to_string(counter++);
}
} // namespace detail

// SimpleFieldData<T> is a generic concrete FieldData<T> backed by shaped
// [dof_holder][component] Kokkos Views.
//
// Ownership of the layout is shared — the layout is typically held by the
// factory that created this field data object.
template <typename T>
class SimpleFieldData
  : public FieldData<T>
  , public DOFAccess<T>
{
public:
  // Value type of the Omega_h tag, if the field has one.
  using tag_value_type = OmegaHValue_t<T>;

  // On an Omega_h layout of a distributed mesh, creates a tag of its own on
  // the mesh, named by detail::NextSimpleFieldDataTagName(), and removes it on
  // destruction.
  SimpleFieldData(std::shared_ptr<const FieldLayout> layout, ValueBasis basis)
    : layout_(std::move(layout)),
      basis_(std::move(basis)),
      device_data_("simple_field_data_device",
                   static_cast<size_t>(layout_->GetNumDofHolder()),
                   static_cast<size_t>(layout_->GetNumComponents()))
  {
    const auto* oh = dynamic_cast<const OmegaHDiscretization*>(
      layout_->GetDiscretization().get());
    const int dim = layout_->GetDOFHolderEntityDim();
    // Synchronize is collective, so every rank of a distributed mesh keeps a
    // tag, whether or not it holds ghosts itself.
    if (oh == nullptr || dim < 0 || oh->GetMesh().comm()->size() == 1) {
      return;
    }
    mesh_ = &oh->GetMesh();
    tag_dim_ = dim;
    tag_name_ = detail::NextSimpleFieldDataTagName();
    owns_tag_ = true;
    const int nc = layout_->GetNumComponents();
    mesh_->add_tag<tag_value_type>(
      tag_dim_, tag_name_, nc,
      Omega_h::Read<tag_value_type>(mesh_->nents(tag_dim_) * nc,
                                    tag_value_type(0)));
  }

  // The tag holds the name of this object, so it is neither copied nor moved.
  SimpleFieldData(const SimpleFieldData&) = delete;
  SimpleFieldData& operator=(const SimpleFieldData&) = delete;
  SimpleFieldData(SimpleFieldData&&) = delete;
  SimpleFieldData& operator=(SimpleFieldData&&) = delete;

  ~SimpleFieldData() noexcept override
  {
    if (owns_tag_ && mesh_->has_tag(tag_dim_, tag_name_)) {
      mesh_->remove_tag(tag_dim_, tag_name_);
    }
  }

  FieldValueType GetValueType() const override
  {
    return ValueTypeOfRank(basis_.Rank());
  }
  const ValueBasis& GetValueBasis() const override { return basis_; }

  const DOFAccess<T>* GetDOFAccess() const noexcept override { return this; }
  DOFAccess<T>* GetDOFAccess() noexcept override { return this; }

  Rank2View<const T, DeviceMemorySpace> GetValues() const override
  {
    return MakeConstRank2View(device_data_);
  }

  void SetValues(Rank2View<const T, DeviceMemorySpace> values) override
  {
    CopyDeviceRank2ViewToDeviceView(device_data_, values);
    if (!tag_name_.empty()) {
      mesh_->set_tag(tag_dim_, tag_name_,
                     detail::OwnedValuesAsTag<T, tag_value_type>(
                       MakeConstRank2View(device_data_),
                       layout_->GetOwnedIndices(), mesh_->nents(tag_dim_),
                       layout_->GetNumComponents()));
    }
  }

  // Fills the tag's ghost entries with their owners' values. Collective; a
  // no-op when the field has no tag.
  void Synchronize() override
  {
    if (!tag_name_.empty()) {
      mesh_->sync_tag(tag_dim_, tag_name_);
    }
  }

  // Name of the Omega_h tag holding the values of every DOF holder, owned and
  // ghost, on entity dimension layout.GetDOFHolderEntityDim() (one value per
  // component per entity); empty when the field has no tag. Its owned entries
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
  bool owns_tag_ = false;
};

} // namespace pcms

#endif // PCMS_SIMPLE_FIELD_DATA_H
