#pragma once

#include <Kokkos_DualView.hpp>

#include <mpicpp.hpp>

#include "dgt_field.hpp"
#include "dgt_subgrid3.hpp"
#include "dgt_tree.hpp"
#include "dgt_view.hpp"

namespace dgt {
namespace ghosting {

class Packing
{

  private:

    template <class T>
    using DualView = typename Kokkos::DualView<T>;

    int m_num_messages = 0;
    int m_num_cells = 0;
    Grid3 m_cell_grid = {0,0,0};
    DualView<Subgrid3*> m_subgrids;
    DualView<int*> m_cell_offsets;
    DualView<int*> m_block_offsets;
    HostPinnedRightView<real***> m_values;

  public:

    Packing() = default;

    void build(
        Grid3 const& cell_grid,
        tree::Adjacencies const& adjs,
        tree::OwnedLeaves const& leaves,
        int const num_max_eqs,
        int const num_modes);

    void pack(Field<real***> const& field);

};

}
}
