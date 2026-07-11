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

  int Serialize(const FieldData<T>& field, const FieldLayout& layout,
                Rank1View<T, HostMemorySpace> buffer,
                Rank1View<const LO, HostMemorySpace> permutation) const override
  {
    if (!rank_participates_) {
      return 0;
    }

    auto const* xgc_field = dynamic_cast<const XGCFieldData<T>*>(&field);
    if (!xgc_field) {
      throw pcms_error("XGCFieldSerializer::Serialize: incompatible FieldData");
    }

    auto data = xgc_field->GetDOFHolderDataHost();
    auto owned = layout.GetOwnedHost();
    if (buffer.size() > 0) {
      const LO num_dof = static_cast<LO>(data.extent(0));
      const LO num_comp = static_cast<LO>(data.extent(1));
      for (LO i = 0; i < num_dof; ++i) {
        for (LO c = 0; c < num_comp; ++c) {
          const LO flat = i * num_comp + c;
          if (owned[flat]) {
            buffer[permutation[flat]] = data(i, c);
          }
        }
      }
    }
    return static_cast<int>(buffer.size());
  }

  void Deserialize(
    FieldData<T>& field, const FieldLayout& layout,
    Rank1View<const T, HostMemorySpace> buffer,
    Rank1View<const LO, HostMemorySpace> permutation) const override
  {
    auto* xgc_field = dynamic_cast<XGCFieldData<T>*>(&field);
    if (!xgc_field) {
      throw pcms_error(
        "XGCFieldSerializer::Deserialize: incompatible FieldData");
    }

    auto current = xgc_field->GetDOFHolderDataHost();
    auto owned = layout.GetOwnedHost();
    const LO num_dof = static_cast<LO>(current.extent(0));
    const LO num_comp = static_cast<LO>(current.extent(1));
    std::vector<T> full_data(current.size());
    for (LO i = 0; i < num_dof; ++i) {
      for (LO c = 0; c < num_comp; ++c) {
        full_data[i * num_comp + c] = current(i, c);
      }
    }
    if (rank_participates_) {
      for (LO i = 0; i < static_cast<LO>(current.size()); ++i) {
        if (owned[i]) {
          full_data[i] = buffer[permutation[i]];
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
