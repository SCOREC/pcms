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
    const auto data = field.GetDOFHolderDataHost().GetValues();
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
    const LO num_dof = layout.GetNumOwnedDofHolder();
    const LO num_comp = layout.GetNumComponents();
    // We initialize with local data so that
    // data that doesn't get received from deserialized
    // doesn't get zeroed out. Long term FIXME
    // is to make Get/SetData operate over only part of the
    // arrays
    const auto current = field.GetDOFHolderDataHost().GetValues();
    Kokkos::View<T*, HostMemorySpace> sorted("sorted", layout.OwnedSize());
    for (LO i = 0; i < num_dof; ++i) {
      for (LO c = 0; c < num_comp; ++c) {
        sorted[i * num_comp + c] = current(i, c);
      }
    }
    for (LO i = 0; i < num_dof; ++i) {
      if (permutation[i] >= 0) {
        for (LO c = 0; c < num_comp; ++c) {
          sorted[i * num_comp + c] = buffer[permutation[i] * num_comp + c];
        }
      }
    }
    field.SetDOFHolderDataUncheckedHost(
      Rank2View<const T, HostMemorySpace>(sorted.data(), num_dof, num_comp));
  }

  virtual ~FieldSerializer() noexcept = default;
};

} // namespace pcms

#endif // PCMS_FIELD_SERIALIZER_H
