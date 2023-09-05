#include <fmt/core.h>

#include "dgt_field.hpp"

namespace dgt {

ModalField::ModalField(
    std::string const& name,
    int const num_blocks,
    int const num_cells,
    int const num_eqs,
    int const num_modes)
{
  using Kokkos::view_alloc;
  using Kokkos::WithoutInitializing;
  auto const outer_alloc = view_alloc(name, WithoutInitializing);
  m_accessor = accessor_t(outer_alloc, num_blocks);
  for (int block = 0; block < num_blocks; ++block) {
    auto const bname = fmt::format("{}[{}]", block, name);
    auto const inner_alloc = view_alloc(bname, WithoutInitializing);
    new (&m_accessor[block]) view_t(inner_alloc, num_cells, num_eqs, num_modes);
  }
}

ModalField::~ModalField()
{
  Kokkos::fence();
  for (std::size_t block = 0; block < m_accessor.extent(0); ++block) {
    m_accessor[block].~view_t();
  }
  m_accessor = accessor_t();
}

}
