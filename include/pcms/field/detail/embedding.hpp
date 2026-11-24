#ifndef PCMS_FIELD_DETAIL_EMBEDDING_HPP
#define PCMS_FIELD_DETAIL_EMBEDDING_HPP

#include <cstddef>

namespace pcms::detail
{

// the search of an embedded geometry (embedding.cpp)
struct EmbeddedSearch;

/// @brief The most sample points of the cell images the search of an embedded
///        geometry evaluates at once
inline std::size_t EmbeddingSampleBatch = std::size_t{1} << 20;

} // namespace pcms::detail

#endif // PCMS_FIELD_DETAIL_EMBEDDING_HPP
