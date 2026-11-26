#ifndef PCMS_FIELD_DETAIL_TRANSFORM_CHAIN_HPP
#define PCMS_FIELD_DETAIL_TRANSFORM_CHAIN_HPP

#include "pcms/field/value_transform.hpp"
#include "pcms/field/value_view.hpp"
#include "pcms/utility/arrays.hpp"
#include "pcms/utility/memory_spaces.hpp"
#include "pcms/utility/detail/arrays.hpp"
#include <Kokkos_Core.hpp>
#include <cstddef>
#include <map>
#include <utility>
#include <vector>

namespace pcms::detail
{

/// @brief The composition T_k ... T_1 of value transforms and its transpose,
///        with scratch storage for the values whose number of components is
///        not that of the end of the chain.
class TransformChain
{
public:
  TransformChain() = default;

  /// @param input basis of the values T_1 takes
  TransformChain(ValueBasis input,
                 std::vector<const ValueTransform*> transforms)
  {
    Reset(std::move(input), std::move(transforms));
  }

  /// @brief Replaces the transforms, keeping the scratch storage.
  void Reset(ValueBasis input, std::vector<const ValueTransform*> transforms)
  {
    transforms_ = std::move(transforms);
    bases_.assign(1, std::move(input));
    for (const auto* transform : transforms_) {
      bases_.push_back(transform->OutputBasis(bases_.back()));
    }
  }

  [[nodiscard]] bool Empty() const noexcept { return transforms_.empty(); }
  [[nodiscard]] const ValueBasis& InputBasis() const noexcept
  {
    return bases_.front();
  }
  [[nodiscard]] const ValueBasis& OutputBasis() const noexcept
  {
    return bases_.back();
  }

  /// @return the storage for the input of Apply(·, out): `out` itself when
  ///         both have as many components
  [[nodiscard]] Rank2View<Real, DeviceMemorySpace> InputStorage(
    Rank2View<Real, DeviceMemorySpace> out) const
  {
    return Storage(0, out);
  }

  /// @brief out = T_k ... T_1 in.
  /// @param out [points][components of OutputBasis()]; may be `in`'s storage
  void Apply(Rank2View<const Real, DeviceMemorySpace> in,
             Rank2View<Real, DeviceMemorySpace> out) const
  {
    if (Empty() && in.data_handle() != out.data_handle()) {
      CopyDeviceRank2ViewToRank2View(out, in);
    }
    for (std::size_t k = 0; k < transforms_.size(); ++k) {
      const auto next = Storage(k + 1, out);
      transforms_[k]->Apply(
        ValueView<const Real, DeviceMemorySpace>(bases_[k], in),
        ValueView<Real, DeviceMemorySpace>(bases_[k + 1], next));
      in = next;
    }
  }

  /// @brief T_1^T ... T_k^T in.
  /// @param in [points][components of OutputBasis()]
  /// @return [points][components of InputBasis()]: `in` itself if the chain
  ///         is empty, scratch storage otherwise
  [[nodiscard]] Rank2View<const Real, DeviceMemorySpace> ApplyTranspose(
    Rank2View<const Real, DeviceMemorySpace> in) const
  {
    const std::size_t columns = in.extent(1);
    for (std::size_t k = transforms_.size(); k > 0; --k) {
      const auto next = Scratch(in.extent(0), Columns(k - 1, columns));
      transforms_[k - 1]->ApplyTranspose(
        ValueView<const Real, DeviceMemorySpace>(bases_[k], in),
        ValueView<Real, DeviceMemorySpace>(bases_[k - 1], next));
      in = next;
    }
    return in;
  }

private:
  // the number of components of values in bases_[k] when those in
  // OutputBasis() have `columns`
  [[nodiscard]] std::size_t Columns(std::size_t k, std::size_t columns) const
  {
    const auto dimension = [](const ValueBasis& basis) {
      return basis.Rank() == 0 ? 0 : basis.system->Dimension();
    };
    const int dk = dimension(bases_[k]);
    const int dout = dimension(bases_.back());
    return dk == dout ? columns : columns / dout * dk;
  }

  // the storage of values in bases_[k]: `out` when they have its components
  [[nodiscard]] Rank2View<Real, DeviceMemorySpace> Storage(
    std::size_t k, Rank2View<Real, DeviceMemorySpace> out) const
  {
    const auto columns = Columns(k, out.extent(1));
    return columns == out.extent(1) ? out : Scratch(out.extent(0), columns);
  }

  [[nodiscard]] Rank2View<Real, DeviceMemorySpace> Scratch(
    std::size_t rows, std::size_t columns) const
  {
    auto& view = scratch_[columns];
    if (view.extent(0) != rows || view.extent(1) != columns) {
      view = Kokkos::View<Real**, DeviceMemorySpace>(
        Kokkos::view_alloc(Kokkos::WithoutInitializing, "transform_chain"),
        rows, columns);
    }
    return MakeRank2View(view);
  }

  std::vector<const ValueTransform*> transforms_;
  // the input basis, then the output basis of each transform
  std::vector<ValueBasis> bases_;
  // by number of columns
  mutable std::map<std::size_t, Kokkos::View<Real**, DeviceMemorySpace>>
    scratch_;
};

} // namespace pcms::detail

#endif // PCMS_FIELD_DETAIL_TRANSFORM_CHAIN_HPP
