#ifndef PCMS_COUPLING_FIELD_EVALUATION_METHODS_H
#define PCMS_COUPLING_FIELD_EVALUATION_METHODS_H
namespace pcms
{

template <int o>
struct Lagrange
{
  static constexpr int order{o};
};
// variable order lagrange interpolation
template <>
struct Lagrange<0>
{
  explicit Lagrange(int o) : order{o} {};
  int order;
};

struct NearestNeighbor
{};

struct Copy
{};

} // namespace pcms

#endif // PCMS_COUPLING_FIELD_EVALUATION_METHODS_H
