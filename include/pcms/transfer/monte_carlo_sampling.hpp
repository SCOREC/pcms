#ifndef PCMS_TRANSFER_MONTE_CARLO_SAMPLING_HPP
#define PCMS_TRANSFER_MONTE_CARLO_SAMPLING_HPP

namespace pcms
{

// Strategy for generating sample points on each target element.
enum class MonteCarloSampling
{
  UniformRandom, // independent uniform draws per element (XorShift64 pool)
};

} // namespace pcms

#endif // PCMS_TRANSFER_MONTE_CARLO_SAMPLING_HPP
