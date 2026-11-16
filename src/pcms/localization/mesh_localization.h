#ifndef PCMS_LOCALIZATION_ADJACENCY_LOCALIZATION_H
#define PCMS_LOCALIZATION_ADJACENCY_LOCALIZATION_H

#include "pcms/localization/localization_factory.h"
#include "pcms/field/evaluator/mls_options.h"

#include <Omega_h_mesh.hpp>

namespace pcms
{

// LocalizationFactory implementation using mesh-adjacency BFS search.
//
// Stores a non-owning reference to the source mesh. The source mesh must
// outlive this factory.
//
// Build(CoordinateView) dispatch:
//   - source_entity_dim == VERT: uses two-mesh searchNeighbors (adjacency BFS).
//   - source_entity_dim != VERT: falls back to N^2 BuildPointCloudSupports.
//
class AdjacencyLocalizationFactory : public LocalizationFactory
{
public:
  AdjacencyLocalizationFactory(Omega_h::Mesh& source_mesh,
                               int source_entity_dim, MLSOptions options)
    : source_mesh_(source_mesh),
      source_entity_dim_(source_entity_dim),
      options_(options)
  {
  }

  SupportResults Build(
    CoordinateView<DeviceMemorySpace> target_coords) const override;
  SupportResults BuildSameMeshCentroidToVertex() const;

private:
  Omega_h::Mesh& source_mesh_;
  int source_entity_dim_;
  MLSOptions options_;
};

} // namespace pcms

#endif // PCMS_LOCALIZATION_ADJACENCY_LOCALIZATION_H
