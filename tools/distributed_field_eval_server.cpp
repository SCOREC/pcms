/**
 * @file distributed_field_eval_server.cpp
 * @brief Standalone evaluation server for distributed field evaluation.
 *
 * Launched as its own single-rank `mpirun` job, separate from the field
 * application(s) that construct a matching
 * pcms::DistributedEvaluationChannel<T> as the Client role. See
 * pcms/src/pcms/coupler/distributed_evaluation_channel.h for the protocol this
 * implements.
 *
 * This server never touches field data; it only routes unresolved query points
 * to candidate owner ranks using a routing partition, then relays their answers
 * back. It does not read a partition file and holds no mesh: the routing
 * partition is assembled from the per-rank bounding boxes the field job
 * publishes at startup (PublishOwnedBounds), so it always matches the
 * distribution the field ranks actually evaluate against.
 *
 * Limitation: one server (and channel name) is needed per distinct field value
 * type / component count. This sketch only supports Real-valued,
 * single-component fields (pcms::DistributedEvaluationChannel<pcms::Real>);
 * extend main() for other PointEvaluatorVariant types as needed.
 */
#include <cstdio>
#include <cstdlib>
#include <memory>
#include <string>
#include <utility>
#include <mpi.h>
#include <redev.h>
#include "pcms/coupler/distributed_evaluation_channel.h"
#include "pcms/coupler/distributed_evaluation_partition.h"

int main(int argc, char** argv)
{
  MPI_Init(&argc, &argv);
  int rank = 0;
  int nproc = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nproc);
  if (nproc != 1) {
    if (rank == 0) {
      std::fprintf(stderr,
                   "distributed_field_eval_server must be launched with "
                   "exactly one rank (got %d)\n",
                   nproc);
    }
    MPI_Abort(MPI_COMM_WORLD, 1);
  }
  if (argc < 5) {
    std::fprintf(
      stderr,
      "Usage: %s <channel_name> <dim> <num_components> <num_rounds>\n"
      "  channel_name must match the name the field job passes to its\n"
      "  DistributedEvaluationChannel, and dim must match its field\n"
      "  dimension. The routing partition is NOT passed in: it is built\n"
      "  from the per-rank boxes the field job publishes at startup.\n"
      "  num_rounds is the total number of evaluation batch rounds to\n"
      "  service before exiting. It must equal the number of Evaluate()\n"
      "  calls the paired field job will make -- NOT the number of field\n"
      "  ranks: one round already fans out across every rank. The two\n"
      "  sides exchange this number and abort immediately on a mismatch\n"
      "  (see DistributedEvaluationChannel::CheckDeclaredRoundCount).\n",
      argv[0]);
    MPI_Abort(MPI_COMM_WORLD, 1);
  }
  const std::string name = argv[1];
  const int dim = std::atoi(argv[2]);
  const int num_components = std::atoi(argv[3]);
  const int num_rounds = std::atoi(argv[4]);

  // Everything that owns redev/ADIOS state is scoped so it is destroyed before
  // MPI_Finalize: AdiosChannel's destructor closes its engines, and closing an
  // ADIOS engine after MPI_Finalize trips "Attempting to use an MPI routine
  // after finalizing MPICH".
  {
    // redev requires a partition to construct a Server-role Redev instance,
    // and the channel handshake ships it to the clients. It must not be empty:
    // ClassPtn::DeserializeModelEntsAndRanks asserts the transferred buffer is
    // a multiple of its (dim, id, rank) stride, which a default-constructed
    // ClassPtn trips. Our actual routing partition is separate -- built from
    // the field job's published boxes below -- so one trivial model entity
    // suffices here.
    redev::ClassPtn dummy_partition(MPI_COMM_WORLD, redev::LOs{0},
                                    redev::ClassPtn::ModelEntVec{{0, 0}});
    redev::Redev redev(MPI_COMM_WORLD,
                       redev::Partition{std::move(dummy_partition)},
                       redev::ProcessType::Server);
    adios2::Params params{{"Streaming", "On"}, {"OpenTimeoutSecs", "60"}};
    auto channel =
      redev.CreateAdiosChannel(name, params, redev::TransportType::BP4);

    pcms::DistributedEvaluationChannel<pcms::Real> eval_channel(
      redev, channel, MPI_COMM_WORLD, name, dim, num_components,
      /*expected_rounds=*/num_rounds);

    // Assemble the routing partition from the boxes the field ranks published.
    std::unique_ptr<pcms::EvaluationPartitionQuery> partition =
      eval_channel.ReceiveRoutingPartition();

    for (int round = 0; round < num_rounds; ++round) {
      eval_channel.RunServerRound(*partition);
    }
  }

  MPI_Finalize();
  return 0;
}
