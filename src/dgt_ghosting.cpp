#include <dgt_cartesian.hpp>
#include <dgt_ghosting.hpp>

#include <dgt_print.hpp> // debug

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

void PackData::build(
    Grid3 const& cell_grid,
    tree::Adjacencies const& adjs,
    tree::OwnedLeaves const& leaves,
    int const num_max_eqs,
    int const num_modes)
{
  m_num_messages = count_messages(adjs, leaves);
  m_blocks = DualView<int*>("Pack::blocks", m_num_messages);
  m_offsets = DualView<int*>("Pack::offsets", m_num_messages);
  m_subgrids = DualView<Subgrid3*>("Pack::subgrids", m_num_messages);
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
        //TODO: modify to 4x4x4 grid
        Vec3<int> const child_ijk = Vec3<int>::zero();
        subgrid = get_coarse_to_fine_cells(
            OWNED, cell_grid, child_ijk, adj.ijk_offset);
      } else {
        throw std::runtime_error("PackData::build - invalid adjacencies");
      }
      m_blocks.h_view[msg_idx] = local_block_idx;
      m_subgrids.h_view[msg_idx] = subgrid;
      m_offsets.h_view[msg_idx] = m_num_cells;
      m_num_cells += generalize(dim, subgrid).size();
      msg_idx++;
    }
    local_block_idx++;
  }
}

}
}
