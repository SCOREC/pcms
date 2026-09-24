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
    // Only owned (rank-exclusive) DOF holders are serialized.
    auto data = field.GetOwnedDOFHolderDataHost();
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
    FieldData<T>& field, const FieldLayout& layout,
    Rank1View<const T, HostMemorySpace> buffer,
    Rank1View<const LO, HostMemorySpace> permutation) const
  {
    const LO num_owned = layout.GetNumOwnedDofHolder();
    const LO num_local = layout.GetNumLocalDofHolder();
    const LO num_comp = layout.GetNumComponents();
    const auto owned_to_local = layout.GetOwnedToLocalHost();
    Kokkos::View<T*, HostMemorySpace> sorted("sorted", layout.LocalSize());
    for (LO o = 0; o < num_owned; ++o) {
      if (permutation[o] >= 0) {
        const LO local = owned_to_local.size() == 0 ? o : owned_to_local(o);
        for (LO c = 0; c < num_comp; ++c) {
          sorted[local * num_comp + c] = buffer[permutation[o] * num_comp + c];
        }
      }
    }
    field.SetDOFHolderDataHost(
      Rank2View<const T, HostMemorySpace>(sorted.data(), num_local, num_comp));
  }

  virtual ~FieldSerializer() noexcept = default;
};

} // namespace pcms

#endif // PCMS_FIELD_SERIALIZER_H
