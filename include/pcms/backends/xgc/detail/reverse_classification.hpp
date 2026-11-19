#ifndef PCMS_BACKENDS_XGC_DETAIL_REVERSE_CLASSIFICATION_HPP
#define PCMS_BACKENDS_XGC_DETAIL_REVERSE_CLASSIFICATION_HPP

#include "pcms/utility/types.hpp"

namespace pcms::xgc::detail
{

struct DimID
{
  LO dim;
  LO id;
  bool operator==(const DimID& other) const
  {
    return (dim == other.dim) && (id == other.id);
  }
};

/// Access to the private members of ReverseClassificationVertex, for PCMS
/// itself and its tests.
struct ReverseClassificationAccess
{
  template <typename RC, typename Data>
  static void Insert(RC& rc, const DimID& key, Data data)
  {
    rc.Insert(key, data);
  }

  template <typename RC>
  [[nodiscard]] static auto Serialize(const RC& rc)
  {
    return rc.Serialize();
  }

  template <typename RC, typename View>
  static void Deserialize(RC& rc, View serialized_data)
  {
    rc.Deserialize(serialized_data);
  }

  template <typename RC>
  [[nodiscard]] static auto Query(const RC& rc, const DimID& geometry)
  {
    return rc.Query(geometry);
  }

  /// the vertices of each geometric entity, in ascending order
  template <typename RC>
  [[nodiscard]] static const auto& Entries(const RC& rc)
  {
    return rc.data_;
  }

  template <typename RC>
  [[nodiscard]] static LO TotalVerts(const RC& rc)
  {
    return rc.total_verts_;
  }

  template <typename RC>
  [[nodiscard]] static bool Equal(const RC& a, const RC& b)
  {
    return a == b;
  }
};

} // namespace pcms::xgc::detail

#endif // PCMS_BACKENDS_XGC_DETAIL_REVERSE_CLASSIFICATION_HPP
