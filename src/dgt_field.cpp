#include <fmt/core.h>

#include "dgt_field.hpp"

namespace dgt {

DGField::DGField(
    std::string const& name,
    int const num_blocks,
    int const num_cells,
    int const num_eqs,
    BasisInfo const& basis)
{
  using Kokkos::view_alloc;
  using Kokkos::WithoutInitializing;
  m_name = name;
  m_basis = basis;
  m_storage.resize(num_blocks);
  m_accessor = accessor_t(name, num_blocks);
  int const nmodes = num_modes(basis.dim, basis.p, basis.tensor);
  for (int block = 0; block < num_blocks; ++block) {
    auto const block_name = fmt::format("{}[{}]", block, name);
    auto const alloc = view_alloc(block_name, WithoutInitializing);
    m_storage[block] = view_t(alloc, num_cells, num_eqs, nmodes);
    m_accessor[block] = uview_t(m_storage[block].data(), num_cells, num_eqs, nmodes);
  }
}

DGField::~DGField()
{
}

}
