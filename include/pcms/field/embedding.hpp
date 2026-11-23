#ifndef PCMS_FIELD_EMBEDDING_HPP
#define PCMS_FIELD_EMBEDDING_HPP

#include "pcms/field/coordinate_map.hpp"
#include "pcms/field/detail/embedding.hpp"
#include "pcms/field/field.hpp"
#include "pcms/field/function_space.hpp"
#include "pcms/geometry/geometry.hpp"
#include <memory>
#include <optional>

namespace pcms
{

/// @brief The coordinate map of a position Function: x(xi) = position(xi),
///        from the coordinate system of the position's (parameter) geometry to
///        the coordinate system of its ValueBasis::Coordinates basis.
///
/// Bind evaluates the position and its Gradient at the parameter points (no
/// inverse): the mapped points and the Jacobian J(p, j, i) = d x_j / d xi_i.
/// Values evaluated at the mapped points flow against the map
/// (BoundCoordinateMap::ValueTransform).
/// @param position a Function with a ValueBasis::Coordinates basis, of a space
///        whose evaluator supports Gradient; the parameter geometry and the
///        coordinate system have at most 3 dimensions
/// @throws pcms_error if the position is null, its values are not
///         coordinates, its space has no Gradient or a dimension is above 3
[[nodiscard]] std::shared_ptr<const CoordinateMap> CreateCoordinateMap(
  std::shared_ptr<const Function<Real>> position);

/// @brief Options of CreateEmbeddedGeometry.
struct EmbeddingOptions
{
  /// the search: `divisions` search bins per axis (0 derives it from the
  /// cell count), `tolerance` the relative slack of the point-in-cell test
  SearchOptions search;
  /// the metric of the parameter system, if it is given rather than induced
  /// by the position: a rank-2 Function on the parameter geometry in the
  /// parameter system, g_ij (both indices covariant) or its inverse g^ij (both
  /// contravariant, inverted at the points). Its volume factor is sqrt(det
  /// g), so a metric given for a coordinate system of a toroidal wedge must
  /// include the ignored direction. Null: the metric induced by the position,
  /// if its coordinate system has one
  std::shared_ptr<const Function<Real>> metric;
  /// the volume factor per unit parameter volume, if given: a scalar Function
  /// on the parameter geometry. Null: sqrt(det g) of the metric, or without
  /// one |det J| (the coordinate measure)
  std::shared_ptr<const Function<Real>> volume;
  /// the inverse of the position or link, if given: a map from its target
  /// system to the parameter system, which Bind uses instead of a search
  /// (none is built). Null: search
  std::shared_ptr<const CoordinateMap> inverse;
};

/// @brief A geometry whose cells are those of a parameter geometry and whose
///        coordinate system is that of a position Function over it, or that
///        of a parent geometry its parameter system is linked into.
class EmbeddedGeometry final : public Geometry
{
public:
  ~EmbeddedGeometry() noexcept override;

  [[nodiscard]] const std::shared_ptr<const CoordinateSystem>&
  GetCoordinateSystem() const noexcept override
  {
    return system_;
  }

  [[nodiscard]] int TopologicalDimension() const noexcept override;

  [[nodiscard]] std::shared_ptr<const CoordinateSystem> GetParameterSystem()
    const noexcept override;

  /// @return the geometry of the cells and of the position
  [[nodiscard]] const std::shared_ptr<const Geometry>& GetParameterGeometry()
    const noexcept
  {
    return parameter_;
  }

  /// @return the position: the map from the parameter system to
  ///         GetCoordinateSystem(), or null when the geometry was created from
  ///         a link into a parent
  [[nodiscard]] const std::shared_ptr<const Function<Real>>& GetPosition()
    const noexcept
  {
    return position_;
  }

  /// @return the geometry whose parameter system the link maps into, or null
  [[nodiscard]] const std::shared_ptr<const Geometry>& GetParent()
    const noexcept
  {
    return parent_;
  }

  /// @return the metric given at creation (EmbeddingOptions::metric), or
  ///         null when it is induced by the position
  [[nodiscard]] const std::shared_ptr<const Function<Real>>& GetMetric()
    const noexcept
  {
    return metric_;
  }

  /// @return the volume factor given at creation (EmbeddingOptions::volume),
  ///         or null when it is that of the metric
  [[nodiscard]] const std::shared_ptr<const Function<Real>>& GetVolume()
    const noexcept
  {
    return volume_;
  }

protected:
  [[nodiscard]] std::shared_ptr<BoundPoints> BindImpl(
    const BindRequest& request) const override;

private:
  friend std::shared_ptr<const EmbeddedGeometry> CreateEmbeddedGeometry(
    std::shared_ptr<const Geometry> parameter_geometry,
    std::shared_ptr<const Function<Real>> position, EmbeddingOptions options);
  friend std::shared_ptr<const EmbeddedGeometry> CreateEmbeddedGeometry(
    std::shared_ptr<const Geometry> parameter_geometry,
    std::shared_ptr<const CoordinateMap> link,
    std::shared_ptr<const Geometry> parent, EmbeddingOptions options);

  [[nodiscard]] static std::shared_ptr<const EmbeddedGeometry> Create(
    std::shared_ptr<const Geometry> parameter,
    std::shared_ptr<const CoordinateMap> link,
    std::shared_ptr<const Geometry> parent,
    std::shared_ptr<const Function<Real>> position, EmbeddingOptions options);

  EmbeddedGeometry(std::shared_ptr<const Geometry> parameter,
                   std::shared_ptr<const Geometry> parent,
                   std::shared_ptr<const CoordinateSystem> system,
                   std::shared_ptr<const Function<Real>> position,
                   std::shared_ptr<const Function<Real>> metric,
                   std::shared_ptr<const Function<Real>> volume,
                   std::shared_ptr<const detail::EmbeddedSearch> search);

  std::shared_ptr<const Geometry> parameter_;
  std::shared_ptr<const Geometry> parent_;
  std::shared_ptr<const Function<Real>> position_;
  std::shared_ptr<const Function<Real>> metric_;
  std::shared_ptr<const Function<Real>> volume_;
  std::shared_ptr<const CoordinateSystem> system_;
  std::shared_ptr<const detail::EmbeddedSearch> search_;
};

/// @brief An embedded geometry over a structured parameter geometry.
///
/// Bind of physical points inverts the position: the images of the cells
/// are binned by their bounding boxes, and Gauss-Newton (Newton for a square
/// Jacobian) finds the parameter coordinates in the candidate cells of a
/// point's bin, modulo the translation identifications of the coordinate
/// system. Bind also accepts points in the parameter system
/// (GetParameterSystem), located directly (no search). The BoundPoints carry
/// the parameter coordinates, so the evaluators of spaces on the parameter
/// geometry evaluate unchanged
/// (CreateEmbeddedSpace). The geometric factors at the points are J, the
/// Gradient of the position, and the given metric (EmbeddingOptions::metric)
/// or the one the position induces, g = J^T G J (G the coordinate system's
/// metric at the image); their volume factor, the physical volume per unit
/// of parameter volume, is the given one (EmbeddingOptions::volume) or that
/// of the metric. In a coordinate system without a metric, the geometry has
/// none unless one is given, and its volume factor is |det J| unless one is
/// given. The geometry does not support
/// OutOfBoundsMode::NEAREST_BOUNDARY, and a Rebind searches anew.
/// @param parameter_geometry a uniform-grid geometry (the cells are boxes in
///        the parameter space), whose coordinate system is a parameter system
///        (csys::CreateParameterSystem) or Cartesian, and not the position's
///        coordinate system
/// @param position a Function on `parameter_geometry` with a
///        ValueBasis::Coordinates basis, whose space evaluates Gradient
/// @param options the options of the search, the metric and the volume factor
/// @throws pcms_error if the geometry is not a 2D or 3D uniform grid, the
///         position is not a Function of coordinates on it with a Gradient,
///         the parameter system has a metric other than the identity or is
///         the position's coordinate system, the inverse does not map the
///         position's coordinate system to the parameter system, or the
///         metric is not a rank-2 Function on the
///         parameter geometry in the parameter system with indices of one
///         variance, or the volume factor is not a scalar Function on the
///         parameter geometry, or neither is given for a parameter geometry
///         of a lower dimension than a coordinate system without a metric, or
///         the position's coordinate system has fewer dimensions than the
///         parameter geometry
[[nodiscard]] std::shared_ptr<const EmbeddedGeometry> CreateEmbeddedGeometry(
  std::shared_ptr<const Geometry> parameter_geometry,
  std::shared_ptr<const Function<Real>> position,
  EmbeddingOptions options = {});

/// @brief An embedded geometry whose parameter system is linked into the
///        parameter system of a parent geometry: a chain of embeddings.
///
/// The coordinate system is the parent's. Bind binds points of it in the
/// parent, then searches the images of the cells under the link for the
/// parent's parameter coordinates, as the position overload does. The factors
/// compose: J = J_parent J_link, the metric is the pullback of the parent's
/// and the volume factor the parent's times sqrt(det g / det g_parent) (both
/// unless given).
/// @param link a map from the parameter geometry's coordinate system to the
///        parent's parameter system
/// @param parent an embedded geometry
/// @throws pcms_error as the position overload, and if the link or the
///         parent is null, the link does not map the parameter geometry's
///         coordinate system to the parent's parameter system, a system has
///         more than 3 dimensions, or the parent's parameter system has fewer
///         dimensions than the parameter geometry
[[nodiscard]] std::shared_ptr<const EmbeddedGeometry> CreateEmbeddedGeometry(
  std::shared_ptr<const Geometry> parameter_geometry,
  std::shared_ptr<const CoordinateMap> link,
  std::shared_ptr<const Geometry> parent, EmbeddingOptions options = {});

/// @brief The space `space`, on the parameter geometry, as a space on the
///        embedded geometry built from it: same layout and values, evaluated
///        at BoundPoints of the embedded geometry.
///
/// Vectors and covectors whose basis is the parameter system are not yet
/// supported (their evaluation throws). Only Value is supported, and the
/// fields are plain DOF holders data (the data of a SimpleFieldData), which is
/// what uniform-grid spaces use. CreateQuadrature is `space`'s rule mapped
/// onto the embedded geometry (MapQuadrature).
/// @throws pcms_error if `embedded` is not an embedded geometry of
///         `space`'s geometry
[[nodiscard]] std::shared_ptr<const FunctionSpace> CreateEmbeddedSpace(
  std::shared_ptr<const EmbeddedGeometry> embedded,
  std::shared_ptr<const FunctionSpace> space);

} // namespace pcms

#endif // PCMS_FIELD_EMBEDDING_HPP
