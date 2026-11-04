#include "pcms/backends/omega_h/omega_h.hpp"
#include "pcms/backends/omega_h/discretization.hpp"
#include "pcms/utility/assert.hpp"

namespace pcms::backend::omega_h
{

std::shared_ptr<const Discretization> CreateDiscretization(
  Omega_h::Mesh& mesh, std::string global_id_name)
{
  return std::make_shared<OmegaHDiscretization>(mesh,
                                                std::move(global_id_name));
}

std::shared_ptr<const Discretization> CreateDiscretization(
  Omega_h::Mesh& mesh, Omega_h::Read<Omega_h::I8> owned_mask,
  std::string global_id_name)
{
  return std::make_shared<OmegaHDiscretization>(
    mesh, std::move(global_id_name), std::move(owned_mask));
}

Omega_h::Mesh& GetMesh(const Discretization& discretization)
{
  const auto* oh = dynamic_cast<const OmegaHDiscretization*>(&discretization);
  if (oh == nullptr) {
    throw pcms_error("pcms::backend::omega_h::GetMesh: the discretization is "
                     "not based on an Omega_h mesh");
  }
  return oh->GetMesh();
}

} // namespace pcms::backend::omega_h
