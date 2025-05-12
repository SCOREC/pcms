#ifndef PCMS_INTERPOLATOR_LOGGER_HPP
#define PCMS_INTERPOLATOR_LOGGER_HPP

#include <pcms/interpolator/pcms_interpolator_aliases.hpp>
#include <Kokkos_Printf.hpp>

namespace pcms
{
enum class LogLevel
{
  INFO,
  WARNING,
  ERROR,
  DEBUG
};

class Logger
{
private:
  size_t selected_league_rank_;

public:
  KOKKOS_INLINE_FUNCTION
  explicit Logger(size_t selected_league_rank)
    : selected_league_rank_(selected_league_rank)
  {
  }

  // log level info
  template <typename... Args>
  KOKKOS_INLINE_FUNCTION void log(const member_type& team, const LogLevel level,
                                  const char* fmt, Args... args)
  {
    if (team.league_rank() == selected_league_rank_) {
      Kokkos::single(Kokkos::PerTeam(team), [&]() {
        Kokkos::printf("[%s] (League %d) ", logLevelToString(level),
                       team.league_rank());
        Kokkos::printf(fmt, args...);
        Kokkos::printf("\n");
      });
    }
  }

  KOKKOS_INLINE_FUNCTION
  void logStruct(const member_type& team, const LogLevel level, const Coord& p,
                 const char* name)
  {

    if (team.league_rank() == selected_league_rank_) {
      Kokkos::single(Kokkos::PerTeam(team), [&]() {
        Kokkos::printf("[%s] (League %d) %s: \n", logLevelToString(level),
                       team.league_rank(), name);
        Kokkos::printf("(%12.6f, %12.6f)", p.x, p.y);
        Kokkos::printf("\n");
      });
    }
  }

  // log array
  KOKKOS_INLINE_FUNCTION
  void logArray(const member_type& team, const LogLevel level,
                const double* array, const int size, const char* name)
  {

    if (team.league_rank() == selected_league_rank_) {
      Kokkos::single(Kokkos::PerTeam(team), [&]() {
        Kokkos::printf("[%s] (League %d) %s: \n", logLevelToString(level),
                       team.league_rank(), name);

        for (int i = 0; i < size; ++i) {
          Kokkos::printf("%12.6f\n", array[i]);
        }
        Kokkos::printf("\n");
      });
    }
  }
  // log scratch vector
  KOKKOS_INLINE_FUNCTION
  void logVector(const member_type& team, const LogLevel level,
                 const ScratchVecView& vector, const char* name) const
  {
    if (team.league_rank() == selected_league_rank_) {
      Kokkos::single(Kokkos::PerTeam(team), [&]() {
        Kokkos::printf("[%s] (League %d) %s: \n", logLevelToString(level),
                       team.league_rank(), name);
        for (int i = 0; i < vector.size(); ++i) {
          Kokkos::printf("%12.6f\n", vector(i));
        }
        Kokkos::printf("\n");
      });
    }
  }

  // log scratch matrix
  KOKKOS_INLINE_FUNCTION
  void logMatrix(const member_type& team, const LogLevel level,
                 const ScratchMatView& matrix, const char* name) const
  {
    if (team.league_rank() == selected_league_rank_) {
      Kokkos::single(Kokkos::PerTeam(team), [&]() {
        Kokkos::printf("[%s] (League %d) %s: \n", logLevelToString(level),
                       team.league_rank(), name);
        for (int i = 0; i < matrix.extent(0); ++i) {
          for (int j = 0; j < matrix.extent(1); ++j) {
            Kokkos::printf("%20.8f", matrix(i, j));
          }
          Kokkos::printf("\n");
        }

        Kokkos::printf("\n");
      });
    }
  }

  // log scalar
  KOKKOS_INLINE_FUNCTION
  void logScalar(const member_type& team, const LogLevel level,
                 const double value, const char* name) const
  {
    if (team.league_rank() == selected_league_rank_) {
      Kokkos::single(Kokkos::PerTeam(team), [&]() {
        Kokkos::printf("[%s] (League %d) %s: \n", logLevelToString(level),
                       team.league_rank(), name);
        Kokkos::printf("%12.6f", value);
        Kokkos::printf("\n");
      });
    }
  }

private:
  KOKKOS_INLINE_FUNCTION
  const char* logLevelToString(LogLevel loglevel) const
  {
    switch (loglevel) {
      case LogLevel::INFO: return "INFO";

      case LogLevel::WARNING: return "WARNING";

      case LogLevel::ERROR: return "ERROR";

      case LogLevel::DEBUG: return "DEBUG";
      default: return "UNKNOWN";
    }
  }
};
} // end namespace pcms
#endif
