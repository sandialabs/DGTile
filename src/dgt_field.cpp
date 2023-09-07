#include <fmt/core.h>

#include "dgt_field.hpp"

namespace dgt {
namespace field {

template <class T>
void CellBase<T>::allocate(
    int const num_blocks,
    int const extent0,
    int const extent1,
    int const extent2)
{
  using Kokkos::view_alloc;
  using Kokkos::WithoutInitializing;
  m_accessor = accessor_t(m_name, num_blocks);
  for (int block = 0; block < num_blocks; ++block) {
    auto const block_name = fmt::format("{}[{}]", m_name, block);
    auto const alloc = view_alloc(block_name, WithoutInitializing);
    m_storage[block] = view_t(alloc, extent0, extent1, extent2);
    m_accessor[block] = uview_t(m_storage[block].data(), extent0, extent1, extent2);
  }
}

template class CellBase<real**>;
template class CellBase<real***>;

Modal::Modal(
    std::string const& name,
    Grid3 const& cell_grid,
    int const num_blocks,
    int const num_comps,
    int const num_modes)
{
  m_name = name;
  int const dim = infer_dimension(cell_grid);
  int const num_cells = generalize(dim, cell_grid).size();
  allocate(num_blocks, num_cells, num_comps, num_modes);
}

}
}
