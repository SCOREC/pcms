#include "pcms/field/function_spaces.hpp"
#include "pcms/backends/point_cloud/evaluator.hpp"
#include "pcms/backends/point_cloud/geometry.hpp"
#include "pcms/field/simple_field_data.hpp"
#include "pcms/utility/assert.hpp"
#include "pcms/utility/detail/common.hpp"

namespace pcms
{
using backend::point_cloud::PointCloudEvaluatorFactory;
using backend::point_cloud::PointCloudGeometry;
} // namespace pcms

namespace pcms
{

std::shared_ptr<const FunctionSpace> CreatePolynomialReconstructionSpace(
  std::shared_ptr<const Geometry> geometry, MLSOptions options)
{
  if (geometry == nullptr) {
    throw pcms_error(
      "CreatePolynomialReconstructionSpace: the geometry must not be null");
  }
  auto cloud_geometry =
    std::dynamic_pointer_cast<const PointCloudGeometry>(geometry);
  if (cloud_geometry == nullptr) {
    throw pcms_error("CreatePolynomialReconstructionSpace: only point-cloud "
                     "geometries are supported");
  }
  auto layout = cloud_geometry->GetSiteLayout();
  auto eval_factory = std::make_shared<PointCloudEvaluatorFactory>(
    std::move(cloud_geometry), layout, options);
  return std::make_shared<FunctionSpace>(
    layout,
    [layout](Type t, ValueBasis basis) -> FieldDataVariant {
      return apply_to_type(t, [&](auto tag) -> FieldDataVariant {
        using T = typename decltype(tag)::type;
        if constexpr (!std::is_same_v<T, double>) {
          throw pcms_error("CreatePolynomialReconstructionSpace: only double "
                           "(Real) is supported");
        } else {
          return std::make_unique<SimpleFieldData<double>>(layout,
                                                           std::move(basis));
        }
      });
    },
    std::move(eval_factory));
}

} // namespace pcms
