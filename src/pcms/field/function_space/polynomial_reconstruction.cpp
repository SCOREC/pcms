#include "pcms/field/function_spaces.hpp"
#include "pcms/backends/point_cloud/evaluator.hpp"
#include "pcms/backends/point_cloud/point_cloud_geometry.hpp"
#include "pcms/field/detail/simple_field_data.hpp"
#include "pcms/utility/assert.hpp"
#include "pcms/utility/detail/common.hpp"

namespace pcms
{
using point_cloud::PointCloudGeometry;
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
  auto eval_factory = point_cloud::CreateMLSEvaluatorFactory(
    std::move(cloud_geometry), layout, options);
  return std::make_shared<FunctionSpace>(layout, eval_factory,
                                         eval_factory->GetGeometry());
}

} // namespace pcms
