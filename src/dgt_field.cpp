#include <fmt/core.h>

#include "dgt_cartesian.hpp"
#include "dgt_field.hpp"

namespace dgt {

using Kokkos::view_alloc;
using Kokkos::WithoutInitializing;

template <class T>
void Field<T>::create(
    std::string const& name,
    int const num_blocks,
    int const extent0,
    int const,
    int const)
{
  m_storage.resize(num_blocks);
}


//  storage.resize(num_blocks);
//  accessor = typename BaseT::accessor_t(name, num_blocks);
//  for (int block = 0; block < num_blocks; ++block) {
//    auto const block_name = fmt::format("{}[{}]", name, block);
//    auto const alloc = view_alloc(block_name, WithoutInitializing);
//    if (view_dim == 2) {
//      storage[block] = typename BaseT::view_t(alloc, n0, n1);
//      accessor[block] = typename BaseT::uview_t(storage[block].data(), n0, n1);
//    }
//    if (view_dim == 3) {
//      storage[block] = typename BaseT::view_t(alloc, n0, n1, n2);
//      accessor[block] = typename BaseT::uview_t(storage[block].data(), n0, n1, n2);
//    }
//  }
//}

}
