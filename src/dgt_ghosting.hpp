#pragma once

#include <Kokkos_DualView.hpp>

#include <mpicpp.hpp>

#include "dgt_subgrid3.hpp"
#include "dgt_tree.hpp"
#include "dgt_view.hpp"

namespace dgt {
namespace ghosting {

class PackData
{

  private:

    template <class T>
    using DualView = typename Kokkos::DualView<T>;

    int m_num_messages = 0;
    int m_num_cells = 0;
    DualView<int*> m_blocks;
    DualView<int*> m_offsets;
    DualView<Subgrid3*> m_subgrids;
    HostPinnedView<real*> m_values;

  public:

    PackData() = default;

    void build(
        Grid3 const& cell_grid,
        tree::Adjacencies const& adjs,
        tree::OwnedLeaves const& leaves,
        int const num_max_eqs,
        int const num_modes);

};

}
}
