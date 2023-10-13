#include <dgt_cartesian.hpp>
#include <dgt_ghosting.hpp>

namespace dgt {
namespace ghosting {

static int count_messages(
    tree::Adjacencies const& adjs,
    tree::OwnedLeaves const& leaves)
{
  int num_messages = 0;
  for (auto const leaf_id : leaves) {
    num_messages += adjs.at(leaf_id).size();
  }
  return num_messages;
}

void Packing::build(
    Grid3 const& cell_grid,
    tree::Adjacencies const& adjs,
    tree::OwnedLeaves const& leaves,
    int const num_max_eqs,
    int const num_modes)
{
  m_num_messages = count_messages(adjs, leaves);
  m_blocks = DualView<int*>("Packing::blocks", m_num_messages);
  m_offsets = DualView<int*>("Packing::offsets", m_num_messages);
  m_subgrids = DualView<Subgrid3*>("Packing::subgrids", m_num_messages);
  m_num_cells = 0;
  int msg_idx = 0;
  int local_block_idx = 0;
  int const dim = infer_dimension(cell_grid);
  for (auto const leaf_id : leaves) {
    for (auto const& adj : adjs.at(leaf_id)) {
      Subgrid3 subgrid;
      if (adj.level_offset == 0) {
        subgrid = get_cells(OWNED, cell_grid, adj.ijk_offset);
      } else if (adj.level_offset == -1) {
        subgrid = get_fine_to_coarse_cells(OWNED, cell_grid, adj.ijk_offset);
      } else if (adj.level_offset == 1) {
        subgrid = get_coarse_to_fine_cells(OWNED, cell_grid, adj.ijk_offset);
      } else {
        throw std::runtime_error("Packing::build - invalid adjacencies");
      }
      m_blocks.h_view[msg_idx] = local_block_idx;
      m_subgrids.h_view[msg_idx] = subgrid;
      m_offsets.h_view[msg_idx] = m_num_cells;
      m_num_cells += generalize(dim, subgrid).size();
      msg_idx++;
    }
    local_block_idx++;
  }
  m_values = HostPinnedRightView<real***>(
      "Packing::values", m_num_cells, num_max_eqs, num_modes);
  Kokkos::deep_copy(m_blocks.d_view, m_blocks.h_view);
  Kokkos::deep_copy(m_offsets.d_view, m_offsets.h_view);
  Kokkos::deep_copy(m_subgrids.d_view, m_subgrids.d_view);
  Kokkos::fence();
}

}
}
