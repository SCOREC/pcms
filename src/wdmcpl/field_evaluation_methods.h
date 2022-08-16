#ifndef WDM_COUPLING_FIELD_EVALUATION_METHODS_H
#define WDM_COUPLING_FIELD_EVALUATION_METHODS_H
namespace wdmcpl
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
  explicit Lagrange(int order) : order{order} {};
  int order;
};

struct NearestNeighbor
{};

} // namespace wdmcpl

#endif // WDM_COUPLING_FIELD_EVALUATION_METHODS_H
