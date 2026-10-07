#ifndef POINT_LOCALIZATION_H
#define POINT_LOCALIZATION_H

#include <any>

#include <ArborX.hpp>
#include <ArborX_Triangle.hpp>
#include <detail/ArborX_PairValueIndex.hpp>
#include <detail/ArborX_AttachIndices.hpp>

#include <Kokkos_Core.hpp>
#include <Omega_h_mesh.hpp>
#include <Omega_h_bbox.hpp>
#include <Omega_h_shape.hpp>
#include <Omega_h_matrix.hpp>
#include <Omega_h_simplex.hpp>

#include "pcms/utility/assert.h"
#include "pcms/utility/types.h"
#include "pcms/utility/arrays.h"
#include "pcms/localization/localization.h"

namespace pcms
{

namespace detail
{

template <int dim>
class Mapping
{};

/**
 * Wrapper base class for the ArborX::BVH and Mappings classes, this is a
 * *rough* workaround for ArborX::BVH being templated on dimension
 */
struct TreeWrapper
{
  template <int dim>
  using Mappings_t =
    Kokkos::View<Mapping<dim>*, Omega_h::ExecSpace::memory_space>;
  template <int dim>
  using Tree_t =
    ArborX::BVH<Omega_h::ExecSpace::memory_space,
                ArborX::PairValueIndex<ArborX::Box<dim, double>, unsigned>>;

  template <int dim>
  TreeWrapper(const Mappings_t<dim>& mappings, const Tree_t<dim>& tree)
  {
    dim_ = dim;
    mappings_ = std::make_shared<std::any>(mappings);
    tree_ = std::make_shared<std::any>(tree);
  }
  /**
   * @brief returns a pointer to an ArborX tree
   */
  template <int dim>
  inline Tree_t<dim> get_tree()
  {
    if (dim != dim_) {
      throw pcms_error("Requested get_tree return "
                       "type does not match internal Tree_t");
    }
    return std::any_cast<Tree_t<dim>>(*tree_);
  }

  /**
   * @brief returns a pointer to a Kokkos::View of Mappings
   */
  template <int dim>
  inline Mappings_t<dim> get_mappings()
  {
    if (dim != dim_) {
      throw pcms_error("Requested get_mappings return "
                       "type does not match internal Mappings_t");
    }
    return std::any_cast<Mappings_t<dim>>(*mappings_);
  };

private:
  int dim_;
  std::shared_ptr<std::any> tree_;
  std::shared_ptr<std::any> mappings_;
};
} // namespace detail

class TreePointSearch : public PointSearch
{
public:
  using Results = PointSearch::Results;
  using Dimensionality = Results::Dimensionality;
  using ExecSpace = PointSearch::ExecSpace;
  using MemorySpace = PointSearch::MemorySpace;

  TreePointSearch(Omega_h::Mesh& mesh)
    : PointSearch(
        PointSearchTolerances("tree point search tolerances", mesh.dim())),
      mesh_(mesh)
  {
    Kokkos::deep_copy(tolerances_, 1e-12);
    tree = make_tree(mesh);
  }
  TreePointSearch(Omega_h::Mesh& mesh, const PointSearchTolerances& tolerances)
    : PointSearch(tolerances), mesh_(mesh), tree(make_tree(mesh))
  {
  }
  ~TreePointSearch() = default;
  /**
   * Given a set of points in global coordinates give the ids of the entities
   * that the points lie within and the parametric coordinates of each point
   * within a triangles or tetrahedra adjacent to the intersected entity.
   * If the point does not lie within any triangle element, the id will
   * be a negative number.
   */
  Results Apply(const CoordinateView<MemorySpace>& coords) const override;
  [[nodiscard]] LO GetOwningElementId(const Results& results, int i) override;
  [[nodiscard]] Kokkos::View<LO*> GetOwningElementIds(
    const Results& results) override;

private:
  std::unique_ptr<detail::TreeWrapper> make_tree(Omega_h::Mesh& mesh) const;
  // Reference to the input mesh
  Omega_h::Mesh& mesh_;
  std::unique_ptr<detail::TreeWrapper> tree;
};

} // namespace pcms
#endif // POINT_LOCALIZATION_H
