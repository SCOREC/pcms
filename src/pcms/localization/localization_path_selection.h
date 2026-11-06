#ifndef PCMS_LOCALIZATION_PATH_SELECTION_H
#define PCMS_LOCALIZATION_PATH_SELECTION_H

#include "pcms/field/field_layout.h"
#include "pcms/utility/entity_types.h"

#include <optional>

namespace pcms::detail
{

enum class LocalizationPath
{
  PointCloudSupports,
  VertexAdjacencySearch,
  CentroidToVertexAdjacencySearch
};

// Returns the unique entity dimension represented by all DOF holders in the
// layout: the layout's GetDOFHolderEntityDim() when it has one DOF holder per
// local mesh entity of a dimension, or else, for a layout without ghosts, the
// dimension of its only populated entity block when its DOF-holder count
// matches that discretization entity count. Otherwise, returns std::nullopt.
inline std::optional<int> GetUniformDofHolderEntityDim(
  const FieldLayout& layout)
{
  auto disc = layout.GetDiscretization();
  if (disc == nullptr) {
    return std::nullopt;
  }

  // A layout with one DOF holder per local mesh entity of a dimension.
  if (const int mesh_entity_dim = layout.GetDOFHolderEntityDim();
      mesh_entity_dim >= 0) {
    return mesh_entity_dim;
  }

  const auto offsets = layout.GetEntOffsets();
  std::optional<int> entity_dim;
  for (int dim = 0; dim < ent_offsets_len - 1; ++dim) {
    const auto block_size = offsets[dim + 1] - offsets[dim];
    if (block_size > 0) {
      if (entity_dim.has_value()) {
        return std::nullopt;
      }
      entity_dim = dim;
    }
  }

  if (!entity_dim.has_value()) {
    return std::nullopt;
  }

  // Only a layout without ghosts covers exactly its discretization's entities.
  if (layout.HasGhosts() ||
      layout.GetNumDofHolder() != disc->GetNumEntities(*entity_dim)) {
    return std::nullopt;
  }

  return entity_dim;
}

inline LocalizationPath SelectLocalizationPath(
  const FieldLayout& source_layout, const FieldLayout* target_layout = nullptr)
{
  auto source_disc = source_layout.GetDiscretization();
  if (source_disc == nullptr) {
    return LocalizationPath::PointCloudSupports;
  }

  auto source_entity_dim = GetUniformDofHolderEntityDim(source_layout);
  if (source_entity_dim.has_value() && *source_entity_dim == Vertex) {
    return LocalizationPath::VertexAdjacencySearch;
  }

  if (target_layout == nullptr) {
    return LocalizationPath::PointCloudSupports;
  }

  auto target_disc = target_layout->GetDiscretization();
  auto target_entity_dim = GetUniformDofHolderEntityDim(*target_layout);
  // The same-mesh shortcut answers for every local target vertex, so it
  // applies only when those are exactly the owned query points.
  if (!target_layout->HasGhosts() && source_entity_dim.has_value() &&
      *source_entity_dim == Face && target_entity_dim.has_value() &&
      *target_entity_dim == Vertex && target_disc != nullptr &&
      source_disc->GetDimension() == 2 &&
      source_disc->SameEntities(*target_disc)) {
    return LocalizationPath::CentroidToVertexAdjacencySearch;
  }

  return LocalizationPath::PointCloudSupports;
}

} // namespace pcms::detail

#endif // PCMS_LOCALIZATION_PATH_SELECTION_H
