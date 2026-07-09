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
    LO counter = 0;
    for (LO i = 0; i < static_cast<LO>(data.size()); ++i) {
      if (!owned[i] || permutation[i] < 0) {
        continue;
      }
      ++counter;
      if (!buffer.empty()) {
        PCMS_ALWAYS_ASSERT(static_cast<size_t>(permutation[i]) < buffer.size());
        buffer[permutation[i]] = data[i];
      }
    }
    return counter;
  }

  virtual void Deserialize(
    FieldData<T>& field, const FieldLayout& layout,
    Rank1View<const T, HostMemorySpace> buffer,
    Rank1View<const LO, HostMemorySpace> permutation) const
  {
    // Seed from the field's current values so DOFs that were not received are
    // preserved rather than overwritten. This matters for masked coupling: a
    // permutation entry < 0 (kUnreceivedDof) means the sender did not include
    // that DOF's GID, so its existing value must be left untouched instead of
    // reading buffer[0].
    auto current = field.GetDOFHolderDataHost();
    Kokkos::View<T*, HostMemorySpace> sorted("sorted", permutation.size());
    for (LO i = 0; i < static_cast<LO>(sorted.size()); ++i) {
      sorted[i] = current[i];
    }
    auto owned = layout.GetOwnedHost();
    for (LO i = 0; i < static_cast<LO>(sorted.size()); ++i) {
      if (owned[i] && permutation[i] >= 0)
        sorted[i] = buffer[permutation[i]];
    }
    field.SetDOFHolderDataHost(make_const_array_view(sorted));
  }

  virtual ~FieldSerializer() noexcept = default;
};

} // namespace pcms

#endif // PCMS_FIELD_SERIALIZER_H
