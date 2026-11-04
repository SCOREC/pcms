#ifndef PCMS_SIMPLE_FIELD_DATA_H
#define PCMS_SIMPLE_FIELD_DATA_H

#include "../field_data.h"
#include "../field_layout.h"
#include "pcms/utility/arrays.h"
#include "pcms/utility/assert.h"
#include <Omega_h_mesh.hpp>
#include "pcms/discretization/discretization/omega_h.hpp"
#include "pcms/utility/memory_spaces.h"
#include "pcms/utility/omega_h_array_utils.h"
#include <Kokkos_Core.hpp>
#include <memory>
#include <string>

namespace pcms
{

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
  // Owned values are the stored data. A layout with ghosts also gets a copy
  // with ghosts (one row per DOF holder, ghosts included: the layout's
  // ownership mask has one entry per DOF holder, ghosts included), which shares
  // the owned values and gets the ghost values from SynchronizeGhosts. Without
  // ghosts it is the owned data.
  SimpleFieldData(std::shared_ptr<const FieldLayout> layout, ValueBasis basis)
    : layout_(std::move(layout)),
      basis_(std::move(basis)),
      owned_data_("simple_field_data_device",
                  static_cast<size_t>(layout_->GetNumDofHolder()),
                  static_cast<size_t>(layout_->GetNumComponents()))
  {
    data_with_ghosts_ =
      !layout_->HasGhosts()
        ? owned_data_
        : Kokkos::View<T**, DeviceMemorySpace>(
            "simple_field_data_local",
            static_cast<size_t>(layout_->GetOwnedHost().size()),
            static_cast<size_t>(layout_->GetNumComponents()));
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
    return MakeConstRank2View(owned_data_);
  }

  void SetValues(Rank2View<const T, DeviceMemorySpace> values) override
  {
    CopyDeviceRank2ViewToDeviceView(owned_data_, values);
    if (layout_->HasGhosts()) {
      ScatterRows(MakeConstRank2View(owned_data_), layout_->GetOwnedIndices(),
                  MakeRank2View(data_with_ghosts_));
    }
  }

  // Values of every DOF holder, owned and ghost, in mesh order. For the
  // evaluators of this backend, which read the DOF holders of local elements.
  Rank2View<const T, DeviceMemorySpace> GetValuesWithGhosts() const
  {
    return MakeConstRank2View(data_with_ghosts_);
  }

  void SynchronizeGhosts() override
  {
    auto* oh = dynamic_cast<const OmegaHDiscretization*>(
      layout_->GetDiscretization().get());
    if (!oh) {
      return;
    }
    const int dim = layout_->GetDOFHolderEntityDim();
    if (dim < 0) {
      return;
    }
    const int nc = layout_->GetNumComponents();
    const LO num_rows = static_cast<LO>(data_with_ghosts_.extent(0));

    // data_with_ghosts_ may be component-major on device; flatten to dof-major
    // order for Omega_h's sync_array, then scatter the synchronized values
    // back. Owned rows keep their values; ghost rows get their owners'.
    Kokkos::View<T*, DeviceMemorySpace> flat("ghost_sync_flat",
                                             static_cast<size_t>(num_rows) *
                                               static_cast<size_t>(nc));
    CopyDeviceRank2ViewToDeviceView(flat,
                                    MakeConstRank2View(data_with_ghosts_));

    SynchronizeOmegaHBlock<T>(oh->GetMesh(), dim, nc, flat);

    CopyDeviceViewToRank2View<T>(data_with_ghosts_, flat);
  }

private:
  std::shared_ptr<const FieldLayout> layout_;
  ValueBasis basis_;
  Kokkos::View<T**, DeviceMemorySpace> owned_data_;
  // Every DOF holder with ghosts, in mesh order; the same view as owned_data_
  // when the layout has no ghosts.
  Kokkos::View<T**, DeviceMemorySpace> data_with_ghosts_;
};

} // namespace pcms

#endif // PCMS_SIMPLE_FIELD_DATA_H
