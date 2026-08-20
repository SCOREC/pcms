#ifndef PCMS_XGC_FIELD_SERIALIZER_H
#define PCMS_XGC_FIELD_SERIALIZER_H

#include "pcms/field/data/xgc.h"
#include "pcms/coupler/field_serializer.h"
#include "pcms/utility/assert.h"
#include "pcms/utility/mpi_type.h"
#include <vector>

namespace pcms
{

template <typename T>
class XGCFieldSerializer : public FieldSerializer<T>
{
public:
  explicit XGCFieldSerializer(MPI_Comm plane_comm,
                              bool rank_participates = true)
    : plane_comm_(plane_comm), rank_participates_(rank_participates)
  {
  }

  int Serialize(const Field<T>& field, Rank1View<T, HostMemorySpace> buffer,
                Rank1View<const LO, HostMemorySpace> permutation) const override
  {
    if (!rank_participates_) {
      return 0;
    }

    auto const* xgc_field =
      dynamic_cast<const XGCFieldData<T>*>(&field.GetData());
    if (!xgc_field) {
      throw pcms_error("XGCFieldSerializer::Serialize: incompatible FieldData");
    }

    auto data = xgc_field->GetDOFHolderDataHost();
    // Per-holder plan: permutation[i] indexes holders; a holder's num_components
    // values form one contiguous block in the wire buffer.
    if (buffer.size() > 0) {
      const LO num_dof = static_cast<LO>(data.extent(0));
      const LO num_comp = static_cast<LO>(data.extent(1));
      for (LO i = 0; i < num_dof; ++i) {
        // A negative permutation entry marks a holder outside the exchange
        // (non-owned, or owned but outside the overlap region); it has no slot.
        if (permutation[i] >= 0) {
          for (LO c = 0; c < num_comp; ++c) {
            buffer[permutation[i] * num_comp + c] = data(i, c);
          }
        }
      }
    }
    return static_cast<int>(buffer.size());
  }

  void Deserialize(
    Field<T>& field, Rank1View<const T, HostMemorySpace> buffer,
    Rank1View<const LO, HostMemorySpace> permutation) const override
  {
    auto* xgc_field = dynamic_cast<XGCFieldData<T>*>(&field.GetData());
    if (!xgc_field) {
      throw pcms_error(
        "XGCFieldSerializer::Deserialize: incompatible FieldData");
    }

    auto current = xgc_field->GetDOFHolderDataHost();
    const auto& layout = field.GetLayout();
    const LO num_dof = static_cast<LO>(current.extent(0));
    const LO num_comp = static_cast<LO>(current.extent(1));
    std::vector<T> full_data(current.size());
    for (LO i = 0; i < num_dof; ++i) {
      for (LO c = 0; c < num_comp; ++c) {
        full_data[i * num_comp + c] = current(i, c);
      }
    }
    if (rank_participates_) {
      for (LO i = 0; i < num_dof; ++i) {
        // A negative permutation entry marks a holder outside the exchange
        // (owned but outside the overlap region); no data was received for it,
        // so it keeps its current value (pre-filled above).
        if (permutation[i] >= 0) {
          for (LO c = 0; c < num_comp; ++c) {
            full_data[i * num_comp + c] = buffer[permutation[i] * num_comp + c];
          }
        }
      }
    }

    MPI_Bcast(full_data.data(), static_cast<int>(full_data.size()),
              pcms::GetMPIType(T{}), 0, plane_comm_);

    xgc_field->SetDOFHolderDataHost(Rank2View<const T, HostMemorySpace>(
      full_data.data(), layout.GetNumOwnedDofHolder(),
      layout.GetNumComponents()));
  }

private:
  MPI_Comm plane_comm_;
  bool rank_participates_;
};

} // namespace pcms

#endif // PCMS_XGC_FIELD_SERIALIZER_H
