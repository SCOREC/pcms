#ifndef PCMS_GEOMETRY_DISTRIBUTED_GEOMETRY_HPP
#define PCMS_GEOMETRY_DISTRIBUTED_GEOMETRY_HPP

#include "pcms/geometry/bound_points.hpp"
#include "pcms/geometry/geometry.hpp"
#include "pcms/utility/memory_spaces.hpp"
#include "pcms/utility/types.hpp"

#include <Kokkos_Core.hpp>
#include <mpi.h>

#include <memory>
#include <vector>

namespace pcms
{

/** Where a query point lies in a partitioned domain. */
struct PointOwner
{
  /** The rank that answers for the point, -1 if no rank contains it. */
  int rank = -1;
  /**
   * The entity of `rank`'s local geometry containing the point, in that
   * rank's numbering, or -1 if unknown.
   */
  LO entity = -1;
};

/**
 * \brief The geometry of a partitioned domain: binds query points anywhere in
 * the whole decomposition, not only in this rank's part.
 *
 * Wraps the geometry of this rank's part. A kind of partition supplies only
 * LocateOwners: which rank owns each query point. Its Bind is collective:
 * every rank calls it, with its own points (possibly none).
 *
 * Fields are evaluated over it by wrapping the evaluator source of a space on
 * the local geometry, see CreateDistributedEvaluator.
 */
class DistributedGeometry : public Geometry
{
public:
  [[nodiscard]] const std::shared_ptr<const CoordinateSystem>&
  GetCoordinateSystem() const noexcept override
  {
    return local_->GetCoordinateSystem();
  }

  /** The geometry of this rank's part. */
  [[nodiscard]] const std::shared_ptr<const Geometry>& GetLocalGeometry()
    const noexcept
  {
    return local_;
  }

  /** The communicator of the partition. */
  [[nodiscard]] MPI_Comm GetComm() const noexcept { return comm_; }

  /** Hints naming the local geometry: its entities are the ones bound. */
  [[nodiscard]] bool TakesHints(const OwnerHints& owners,
                                LO num_points) const noexcept override
  {
    return owners.Names(*local_, num_points);
  }

  /**
   * The owner of each of `points` (one row per point, in the local
   * geometry's coordinate system). Collective.
   * \param owners claimed owners of `points` naming the local geometry, or
   *        empty; a kind may start its search from them (verified)
   * \return one owner per point
   */
  [[nodiscard]] virtual std::vector<PointOwner> LocateOwners(
    Kokkos::View<const Real**, DeviceMemorySpace> points,
    const OwnerHints& owners) const = 0;

  /** The dimension of the entities LocateOwners returns. */
  [[nodiscard]] virtual int OwnerEntityDim() const noexcept = 0;

protected:
  /** \throws pcms_error if `local` is null */
  DistributedGeometry(std::shared_ptr<const Geometry> local, MPI_Comm comm);

  [[nodiscard]] std::shared_ptr<BoundPoints> BindImpl(
    const BindRequest& request) const override;

private:
  std::shared_ptr<const Geometry> local_;
  MPI_Comm comm_;
};

} // namespace pcms

#endif // PCMS_GEOMETRY_DISTRIBUTED_GEOMETRY_HPP
