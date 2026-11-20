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

/// @brief A translation symmetry of the domain of an embedded geometry: the
///        domain is invariant under adding `length` to coordinate `axis` of
///        the position's coordinate system, and its cells cover one `length`
///        (e.g. a 1/zperiod toroidal wedge: axis phi, length 2 pi / zperiod).
struct TranslationSymmetry
{
  int axis = 0;
  Real length = 0;
};

/// @brief Options of CreateEmbeddedGeometry.
struct EmbeddingOptions
{
  /// the search: `divisions` search bins per axis (0 derives it from the
  /// cell count), `tolerance` the relative slack of the point-in-cell test
  SearchOptions search;
  /// a symmetry of the domain: a query is reduced by it before the search
  /// (as by a period of the coordinate system), and the BoundPoints keep the
  /// reduced coordinates as their located points. The reduction replaces the
  /// period of the axis, so a symmetry must divide it; none by default
  std::optional<TranslationSymmetry> symmetry;
  /// the metric of the parameter system, if it is given rather than induced
  /// by the position: a rank-2 Function on the parameter geometry in the
  /// parameter system, g_ij (both indices covariant) or its inverse g^ij (both
  /// contravariant, inverted at the points). Its volume factor is sqrt(det
  /// g), so a metric given for a coordinate system reduced by a symmetry must
  /// include the ignored direction. Null: the metric induced by the position
  std::shared_ptr<const Function<Real>> metric;
};

/// @brief A geometry whose cells are those of a parameter geometry and whose
///        coordinate system is that of a position Function over it.
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
  ///         GetCoordinateSystem()
  [[nodiscard]] const std::shared_ptr<const Function<Real>>& GetPosition()
    const noexcept
  {
    return position_;
  }

  /// @return the metric given at creation (EmbeddingOptions::metric), or
  ///         null when it is induced by the position
  [[nodiscard]] const std::shared_ptr<const Function<Real>>& GetMetric()
    const noexcept
  {
    return metric_;
  }

protected:
  [[nodiscard]] std::shared_ptr<BoundPoints> BindImpl(
    const BindRequest& request) const override;

  [[nodiscard]] Quadrature CreateQuadrature(
    const detail::QuadratureRule& rule) const override;

private:
  friend std::shared_ptr<const EmbeddedGeometry> CreateEmbeddedGeometry(
    std::shared_ptr<const Geometry> parameter_geometry,
    std::shared_ptr<const Function<Real>> position, EmbeddingOptions options);

  EmbeddedGeometry(std::shared_ptr<const Geometry> parameter,
                   std::shared_ptr<const Function<Real>> position,
                   std::shared_ptr<const Function<Real>> metric,
                   std::shared_ptr<const detail::EmbeddedSearch> search);

  std::shared_ptr<const Geometry> parameter_;
  std::shared_ptr<const Function<Real>> position_;
  std::shared_ptr<const Function<Real>> metric_;
  std::shared_ptr<const CoordinateSystem> system_;
  std::shared_ptr<const detail::EmbeddedSearch> search_;
};

/// @brief An embedded geometry over a structured parameter geometry.
///
/// Bind of physical points inverts the position: the images of the cells
/// are binned by their bounding boxes, and Gauss-Newton (Newton for a square
/// Jacobian) finds the parameter coordinates in the candidate cells of a
/// point's bin, modulo the periods of the coordinate system and the options'
/// symmetry. Bind also accepts points in the parameter system
/// (GetParameterSystem), located directly (no search). The BoundPoints carry
/// the parameter coordinates, so the evaluators of spaces on the parameter
/// geometry evaluate unchanged
/// (CreateEmbeddedSpace). The geometric factors at the points are J, the
/// Gradient of the position, and the given metric (EmbeddingOptions::metric)
/// or the one the position induces, g = J^T G J (G the coordinate system's
/// metric at the image); their volume factor is the physical volume per unit
/// of parameter volume. The geometry does not support
/// OutOfBoundsMode::NEAREST_BOUNDARY, and a Rebind searches anew.
/// @param parameter_geometry a uniform-grid geometry (the cells are boxes in
///        the parameter space), whose coordinate system is a parameter system
///        (csys::CreateParameterSystem) or at least not the position's
///        coordinate system
/// @param position a Function on `parameter_geometry` with a
///        ValueBasis::Coordinates basis, whose space evaluates Gradient
/// @param options the options of the search, the symmetry and the metric
/// @throws pcms_error if the geometry is not a 2D or 3D uniform grid, the
///         position is not a Function of coordinates on it with a Gradient,
///         the parameter system is the position's coordinate system, the
///         symmetry has an axis outside the coordinate system or a
///         non-positive length, or the metric is not a rank-2 Function on the
///         parameter geometry in the parameter system with indices of one
///         variance
[[nodiscard]] std::shared_ptr<const EmbeddedGeometry> CreateEmbeddedGeometry(
  std::shared_ptr<const Geometry> parameter_geometry,
  std::shared_ptr<const Function<Real>> position,
  EmbeddingOptions options = {});

/// @brief The space `space`, on the parameter geometry, as a space on the
///        embedded geometry built from it: same layout and values, evaluated
///        at BoundPoints of the embedded geometry.
///
/// Vectors and covectors whose basis is the parameter system are evaluated in
/// the components of the embedded geometry's coordinate system (the Pushforward
/// of the position's Jacobian: vectors by J, covectors by J^-T, which needs as
/// many parameters as coordinates). Only Value is supported, and the fields are
/// plain DOF holders data (the data of a SimpleFieldData), which is what
/// uniform-grid spaces use.
/// @throws pcms_error if `embedded` is not an embedded geometry of
///         `space`'s geometry
[[nodiscard]] std::shared_ptr<const FunctionSpace> CreateEmbeddedSpace(
  std::shared_ptr<const EmbeddedGeometry> embedded,
  std::shared_ptr<const FunctionSpace> space);

} // namespace pcms

#endif // PCMS_FIELD_EMBEDDING_HPP
