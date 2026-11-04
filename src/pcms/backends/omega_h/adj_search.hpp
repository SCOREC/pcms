#ifndef PCMS_BACKENDS_OMEGA_H_ADJ_SEARCH_HPP
#define PCMS_BACKENDS_OMEGA_H_ADJ_SEARCH_HPP

#include <pcms/localization/point_search.hpp>
#include <pcms/localization/mls_support_helpers.hpp>
#include <pcms/backends/omega_h/queue_visited.hpp>
#include <pcms/localization/support_results.hpp>
#include <pcms/utility/mesh_geometry.hpp>
#include <pcms/utility/print.hpp>

#include <Kokkos_MathematicalFunctions.hpp>

namespace pcms::backend::omega_h
{

static constexpr int max_dim = 3;

class FindSupports
{
private:
  Omega_h::Mesh& source_mesh;
  Omega_h::Mesh& target_mesh;

public:
  FindSupports(Omega_h::Mesh& source_mesh_, Omega_h::Mesh& target_mesh_)
    : source_mesh(source_mesh_), target_mesh(target_mesh_){};

  FindSupports(Omega_h::Mesh& mesh_) : source_mesh(mesh_), target_mesh(mesh_){};

  void adjBasedSearch(Omega_h::Write<Omega_h::LO>& supports_ptr,
                      Omega_h::Write<Omega_h::LO>& nSupports,
                      Omega_h::Write<Omega_h::LO>& support_idx,
                      Omega_h::Write<Omega_h::Real>& radii2,
                      bool is_build_csr_call);

  void adjBasedSearchCentroidNodes(Omega_h::Write<Omega_h::LO>& supports_ptr,
                                   Omega_h::Write<Omega_h::LO>& nSupports,
                                   Omega_h::Write<Omega_h::LO>& support_idx,
                                   Omega_h::Write<Omega_h::Real>& radii2,
                                   bool is_build_csr_call);
};

pcms::localization::SupportResults searchNeighbors(Omega_h::Mesh& source_mesh,
                               const Omega_h::Reals& target_coords,
                               Omega_h::LO n_targets,
                               Omega_h::Real& cutoffDistance,
                               Omega_h::LO min_req_support = 12,
                               Omega_h::LO max_allowed_support = 36,
                               bool adapt_radius = true);

pcms::localization::SupportResults searchNeighbors(Omega_h::Mesh& source_mesh,
                               Omega_h::Mesh& target_mesh,
                               Omega_h::Real& cutoffDistance,
                               Omega_h::LO min_req_support = 12,
                               Omega_h::LO max_allowed_support = 36,
                               bool adapt_radius = true);

pcms::localization::SupportResults searchNeighbors(Omega_h::Mesh& mesh,
                               Omega_h::Real cutoffDistance,
                               Omega_h::LO min_support = 12,
                               bool adapt_radius = true);
} // namespace pcms::backend::omega_h

#endif // PCMS_BACKENDS_OMEGA_H_ADJ_SEARCH_HPP
