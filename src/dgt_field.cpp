#include <fmt/core.h>

#include "dgt_field.hpp"

namespace dgt {

Field::Field(
    std::string const& name,
    int const num_blocks,
    int const extent0,
    int const extent1,
    int const extent2)
{
  using Kokkos::view_alloc;
  using Kokkos::WithoutInitializing;

  m_name = name;
  m_storage.resize(num_blocks);
  m_accessor = accessor_t(name, num_blocks);
  for (int block = 0; block < num_blocks; ++block) {
    auto const block_name = fmt::format("{}[{}]", block, name);
    auto const alloc = view_alloc(block_name, WithoutInitializing);
    m_storage[block] = view_t(alloc, extent0, extent1, extent2);
    m_accessor[block] = uview_t(m_storage[block].data(), extent0, extent1, extent2);
  }
}

}
