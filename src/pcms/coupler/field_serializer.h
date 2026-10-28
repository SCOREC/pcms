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
  virtual int Serialize(const Field<T>& field,
                        Rank1View<T, HostMemorySpace> buffer,
                        Rank1View<const LO, HostMemorySpace> permutation) const
  {
    // Only owned (rank-exclusive) DOF holders are serialized.
    const auto data = field.GetOwnedDOFHolderDataHost().GetValues();
    if (buffer.size() > 0) {
      const LO num_dof = static_cast<LO>(data.extent(0));
      const LO num_comp = static_cast<LO>(data.extent(1));
      for (LO i = 0; i < num_dof; ++i) {
        if (permutation[i] >= 0) {
          for (LO c = 0; c < num_comp; ++c) {
            buffer[permutation[i] * num_comp + c] = data(i, c);
          }
        }
      }
    }
    return static_cast<int>(buffer.size());
  }

  virtual void Deserialize(
    Field<T>& field, Rank1View<const T, HostMemorySpace> buffer,
    Rank1View<const LO, HostMemorySpace> permutation) const
  {
    const auto& layout = field.GetLayout();
    const LO num_owned = layout.GetNumOwnedDofHolder();
    const LO num_local = layout.GetNumLocalDofHolder();
    const LO num_comp = layout.GetNumComponents();
    const auto owned_to_local = layout.GetOwnedToLocalHost();
    // We initialize with local data so that data that doesn't get received
    // from deserialize doesn't get zeroed out. Long term FIXME is to make
    // Get/SetData operate over only part of the arrays.
    const auto current = field.GetDOFHolderDataHost().GetValues();
    Kokkos::View<T*, HostMemorySpace> sorted("sorted", layout.LocalSize());
    for (LO i = 0; i < num_local; ++i) {
      for (LO c = 0; c < num_comp; ++c) {
        sorted[i * num_comp + c] = current(i, c);
      }
    }
    for (LO o = 0; o < num_owned; ++o) {
      if (permutation[o] >= 0) {
        const LO local = owned_to_local.size() == 0 ? o : owned_to_local(o);
        for (LO c = 0; c < num_comp; ++c) {
          sorted[local * num_comp + c] = buffer[permutation[o] * num_comp + c];
        }
      }
    }
    field.SetDOFHolderDataUncheckedHost(
      Rank2View<const T, HostMemorySpace>(sorted.data(), num_local, num_comp));
  }

  virtual ~FieldSerializer() noexcept = default;
};

} // namespace pcms

#endif // PCMS_FIELD_SERIALIZER_H
