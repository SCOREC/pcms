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
    // owned/permutation/buffer index the flat node-major wire order
    // (dof * num_components + component); read the field values through the
    // 2D [dof][comp] view.
    if (buffer.size() > 0) {
      const LO num_dof = static_cast<LO>(data.extent(0));
      const LO num_comp = static_cast<LO>(data.extent(1));
      for (LO i = 0; i < num_dof; ++i) {
        for (LO c = 0; c < num_comp; ++c) {
          const LO flat = i * num_comp + c;
          if (owned[flat])
            buffer[permutation[flat]] = data(i, c);
        }
      }
    }
    return static_cast<int>(data.size());
  }

  virtual void Deserialize(
    FieldData<T>& field, const FieldLayout& layout,
    Rank1View<const T, HostMemorySpace> buffer,
    Rank1View<const LO, HostMemorySpace> permutation) const
  {
    Kokkos::View<T*, HostMemorySpace> sorted("sorted", permutation.size());
    auto owned = layout.GetOwnedHost();
    for (LO i = 0; i < static_cast<LO>(sorted.size()); ++i) {
      if (owned[i])
        sorted[i] = buffer[permutation[i]];
    }
    field.SetDOFHolderDataHost(Rank2View<const T, HostMemorySpace>(
      sorted.data(), layout.GetNumOwnedDofHolder(), layout.GetNumComponents()));
  }

  virtual ~FieldSerializer() noexcept = default;
};

} // namespace pcms

#endif // PCMS_FIELD_SERIALIZER_H
