#ifndef PCMS_FIELD_SERIALIZER_H
#define PCMS_FIELD_SERIALIZER_H

#include "pcms/field/field.h"
#include "pcms/utility/arrays.h"
#include "pcms/utility/memory_spaces.h"
#include "pcms/utility/types.h"
#include <Kokkos_Core.hpp>

namespace pcms
{

template <typename T>
class FieldSerializer
{
public:
  virtual int Serialize(const FieldData<T>& field, const FieldLayout& layout,
                        Rank1View<T, HostMemorySpace> buffer,
                        Rank1View<const LO, HostMemorySpace> permutation) const
  {
    auto data = field.GetDOFHolderDataHost();
    auto owned = layout.GetOwnedHost();
    // The exchange plan is per DOF holder: owned[i] and permutation[i] are
    // indexed by holder. All num_components components of a holder share its
    // location, so they occupy one contiguous block permutation[i]*num_comp in
    // the wire buffer.
    if (buffer.size() > 0) {
      const LO num_dof = static_cast<LO>(data.extent(0));
      const LO num_comp = static_cast<LO>(data.extent(1));
      for (LO i = 0; i < num_dof; ++i) {
        if (owned[i]) {
          for (LO c = 0; c < num_comp; ++c) {
            buffer[permutation[i] * num_comp + c] = data(i, c);
          }
        }
      }
    }
    return static_cast<int>(buffer.size());
  }

  virtual void Deserialize(
    FieldData<T>& field, const FieldLayout& layout,
    Rank1View<const T, HostMemorySpace> buffer,
    Rank1View<const LO, HostMemorySpace> permutation) const
  {
    const LO num_dof = layout.GetNumOwnedDofHolder();
    const LO num_comp = layout.GetNumComponents();
    Kokkos::View<T*, HostMemorySpace> sorted("sorted", layout.OwnedSize());
    auto owned = layout.GetOwnedHost();
    for (LO i = 0; i < num_dof; ++i) {
      if (owned[i]) {
        for (LO c = 0; c < num_comp; ++c) {
          sorted[i * num_comp + c] = buffer[permutation[i] * num_comp + c];
        }
      }
    }
    field.SetDOFHolderDataHost(
      Rank2View<const T, HostMemorySpace>(sorted.data(), num_dof, num_comp));
  }

  virtual ~FieldSerializer() noexcept = default;
};

} // namespace pcms

#endif // PCMS_FIELD_SERIALIZER_H
