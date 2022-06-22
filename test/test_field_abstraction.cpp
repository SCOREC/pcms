#include <catch2/catch.hpp>
#include <wdmcpl/types.h>
//#include <wdmcpl/coordinate_systems.h>
#include <wdmcpl/coordinate_transform.h>
using wdmcpl::Cartesian;
using wdmcpl::Coordinate;
using wdmcpl::CoordinateTransform;
using wdmcpl::Cylindrical;
using wdmcpl::DataTransform;
using wdmcpl::GO;
using wdmcpl::LO;
using wdmcpl::Real;
using wdmcpl::ScalarData;

namespace wdmcpl
{
template <typename InputCoordinateIterator, typename OutputCoordinateIterator>
constexpr OutputCoordinateIterator CoordinateTransform(
  InputCoordinateIterator first, InputCoordinateIterator last,
  OutputCoordinateIterator out) noexcept
{
  using FromCoordinateSystem =
    typename InputCoordinateIterator::CoordinateSystem;
  using ToCoordinateSystem =
    typename OutputCoordinateIterator::CoordinateSystem;
  // implementation that doesn't require construction
  if constexpr (std::is_same_v<FromCoordinateSystem, ToCoordinateSystem>) {
    return last;
  }
  for (auto it = first; it != last; ++it) {
    // Coordinate transformation functions
    // if the InputIterator::value_type is a Coordinate
    *out = CoordinateTransform<ToCoordinateSystem, FromCoordinateSystem>(*it);
    // if InputIteraterator::value_type is a numeric type
    // *out = detail::CoordinateTransform<ToCoordinateSystem, FromCoordinateSystem>(*it);
    ++out;
  }
  // if output iterator and input iterator are
  // random access iterator parallelize w/ parallel_for
  //parallel_for(std::distance(last,first),  [](size_t i) {
  //  auto it = std::next(first, i); // only efficient with random access iterator
  //  out_it = std::next(out, i);
  //  *out_it = CoordinateTransform<ToCoordinateSystem, FromCoordinateSystem>(*it);
  // });
  // return std::next(out, std::distance(first,last));
}
} // namespace wdmcpl

// field evaluation methods:
// lagrange
// bspline
// SPR patch recovery
// mean element value
// max element value?

struct Lagrange
{
  Lagrange() = default;
  Lagrange(int order) : order(order) {}
  int order{1};
};
/* do we want to be able to swap between compile time and runtime order?
template <>
struct Lagrange <0>{
  Lagrange(int order) : order(order) {}
  int order;
};
*/
struct Node;
struct NodeIterator
{
  // dummy iterator for now auto for loop works...
  Node* begin();
  Node* end();
  const Node* begin() const;
  const Node* end() const;
};
/*
class FieldConcept
{
public:
  // these types can be changed in implementation. Other components of interface
remain fixed using CoordinateSystem = Cartesian; using DataT = Real; using
NodeHandleT = LO; using CoordT = Real;

  // methods with return types can return by value or const and non-const
reference template <typename MethodT> ScalarData<DataT> evaluate(MethodT method,
Coordinate<CoordinateSystem, CoordT> coordinate) const; template <typename
MethodT> ScalarData<std::vector<DataT>> evaluate(MethodT method, const
Coordinate<CoordinateSystem, std::vector<CoordT>>& coordinates) const;

  // NodeIterator must abide by InputIterator concept
  NodeIterator nodes();

  // set functions for higher order
  void set(const NodeHandleT & node, const ScalarData<DataT>& data);
  void set(const std::vector<NodeHandleT>& ids, const
ScalarData<std::vector<DataT>>& data);

  Coordinate<CoordinateSystem, CoordT> get_coordinate(const Node& nd) const;
  Coordinate<CoordinateSystem, std::vector<CoordT>> get_nodal_coordinates();

  std::vector<NodeHandleT> get_node_handles();
};
*/

template <typename CoordinateSystemT = Cartesian>
class SampleField
{
public:
  using CoordinateSystem = CoordinateSystemT;
  // using StorageFormat = ScalarData<std::vector<Real>>;
  using StorageFormat = std::vector<Real>;
  // if we never need to pass state into evaluate
  // method can become an enum type this would have the beneficial
  // aspect of warning for unhandled method types. Although less flexible
  // template <typename MethodT>
  // ScalarData<Real> evaluate(MethodT method,
  //                          Coordinate<CoordinateSystem, Real> coordinate)
  [[nodiscard]] ScalarData<Real> evaluate(
    Lagrange method, Coordinate<CoordinateSystem, Real> coordinate) const
  {
    // transform coordinate into current coordinate system
    return {};
  }
  // template <typename MethodT>
  // ScalarData<std::vector<Real>> evaluate(
  //  MethodT method,
  //  Coordinate<CoordinateSystem, std::vector<Real> > coordinates)

  template <typename MethodT, typename InputIt, typename OutputIt>
  constexpr OutputIt evaluate(MethodT method, InputIt first, InputIt last,
                              OutputIt out_first)
  {
    // transform coordinate into current coordinate system
    return {};
  }

  struct Node
  {
    GO id;
    // Coordinate<CoordinateSystem> coordinate;
  };
  Coordinate<CoordinateSystem, Real> get_coordinate(const Node& nd) const
  {
    return {};
  }

  struct NodeIterator
  {
    // dummy iterator for now auto for loop works...
    Node* begin() { return {}; }
    Node* end() { return {}; }
    const Node* begin() const { return {}; }
    const Node* end() const { return {}; }
  };
  NodeIterator nodes() { return {}; }
  void set(const Node& node, ScalarData<Real> result) {}
  void set(const std::vector<LO>& ids, ScalarData<std::vector<Real>> result) {}
  Coordinate<CoordinateSystem, std::vector<Real>> get_nodal_coordinates()
  {
    return {};
  }
  std::vector<LO> get_node_handles() { return {}; }
};

template <typename T>
struct HasBatchMethods : std::true_type
{
};

// InputIterator and output iterator need to meet the requirements of
// LegacyInput iterator such that std::next and std::distance work
template <typename Field, typename Method, typename InputIt, typename OutputIt>
OutputIt evaluate(Field field, Method method, InputIt first, InputIt last,
                  OutputIt out_first);

template <typename SourceField, typename TargetField, typename Method>
void transfer_field_data(const Method& method, const SourceField& source_field,
                         TargetField& target_field)
{

  using std::begin;
  using std::end;
  auto coordinates = target_field.get_nodal_coordinates(); // name???
  // in place transform coordinates
  CoordinateTransform<typename SourceField::CoordinateSystem>(begin(coordinates), end(coordinates), begin(transformed_coordinates));
  auto node_ids = target_field.get_node_handles(); // note handles can be a list of ids or any
  typename TargetField::StorageType results{std::distance(begin(coordinates), end(coordinates))};
  evaluate(source_field, method, begin(coordinates), end(coordinates),begin(results));
  // Interface forces a copy
  DataTransform<typename TargetField::CoordinateSystem>(begin(results), end(results), begin(transformed_results));
  target_field.set(node_ids, target_result);
}
/*
// templates left unconstrained for slideware
template <typename SourceField, typename TargetField, typename Method>
void transfer_field_data(const Method& method, const SourceField& source_field,
                     TargetField& target_field)
{
using std::begin;
using std::end;
// if the meshes have API for data parallel operation
if constexpr (HasBatchMethods<SourceField>::value &&
            HasBatchMethods<TargetField>::value) {
auto coordinates = target_field.get_nodal_coordinates(); // name???
auto transformed_coordinates =
  CoordinateTransform<typename SourceField::CoordinateSystem>(coordinates);
auto node_ids =
  target_field
    .get_node_handles(); // note handles can be a list of ids or any
//auto result = source_field.evaluate(method, transformed_coordinates);
// NEED TAGGED ITERATOR...
ScalarData<std::vector<Real>> // TODO make type generic
auto result = source_field.evaluate(method,
                                    begin(transformed_coordinates),
                                    end(transformed_coordinates),);

auto&& target_result =
  DataTransform<typename TargetField::CoordinateSystem>(result);
target_field.set(node_ids, target_result);

}
else {
for (const auto& node : target_field.nodes()) {
  auto coordinate =
    CoordinateTransform<typename SourceField::CoordinateSystem>(
      target_field.get_coordinate(node));
  auto&& result = source_field.evaluate(method, coordinate);
  auto&& target_result =
    DataTransform<typename TargetField::CoordinateSystem>(result);
  target_field.set(node, target_result);
}
}
}
*/

TEST_CASE("evaluate scalar field")
{
  SampleField<Cartesian> target_field;
  SampleField<Cylindrical> source_field;
  Lagrange lagrange_quadratic{2};
  transfer_field_data(lagrange_quadratic, source_field, target_field);
}

/*
TEST_CASE("batched evaluate field")
{
  SampleField target_field;
  SampleField<Cylindrical> source_field;
  Lagrange method{2};
  auto coordinates = target_field.get_nodal_coordinates(); // name???
  auto transformed_coordinates =
    CoordinateTransform<decltype(source_field)::CoordinateSystem>(coordinates);
  auto node_ids =
    target_field.get_node_handles(); // note handles can be a list of ids or any
  // auto result = source_field.evaluate(method, coordinates);
  auto result = source_field.evaluate(method, transformed_coordinates);

  auto&& target_result =
    DataTransform<decltype(target_field)::CoordinateSystem>(result);
  target_field.set(node_ids, target_result);
}
*/