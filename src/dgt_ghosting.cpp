#include <dgt_cartesian.hpp>
#include <dgt_for_each.hpp>
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
  int const num_blocks = leaves.size();
  m_num_messages = count_messages(adjs, leaves);
  m_cell_offsets = DualView<int*>("Packing::offsets", m_num_messages);
  m_block_offsets = DualView<int*>("Packing::block_offsets", num_blocks + 1);
  m_subgrids = DualView<Subgrid3*>("Packing::subgrids", m_num_messages);
  m_num_cells = 0;
  int msg_idx = 0;
  int block_idx = 0;
  int const dim = infer_dimension(cell_grid);
  for (auto const leaf_id : leaves) {
    m_block_offsets.h_view[block_idx] = msg_idx;
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
      m_subgrids.h_view[msg_idx] = subgrid;
      m_cell_offsets.h_view[msg_idx] = m_num_cells;
      m_num_cells += generalize(dim, subgrid).size();
      msg_idx++;
    }
    m_block_offsets.h_view[block_idx+1] = msg_idx;
    block_idx++;
  }
  m_values = HostPinnedRightView<real***>(
      "Packing::values", m_num_cells, num_max_eqs, num_modes);
  Kokkos::deep_copy(m_subgrids.d_view, m_subgrids.d_view);
  Kokkos::deep_copy(m_cell_offsets.d_view, m_cell_offsets.h_view);
  Kokkos::deep_copy(m_block_offsets.d_view, m_block_offsets.h_view);
  Kokkos::fence();
}

void Packing::pack(Field<real***> const& field)
{
  auto const U = field.get();
  auto const subgrids = m_subgrids.d_view;
  auto const offsets_c = m_cell_offsets.d_view;
  auto const offsets_b = m_block_offsets.d_view;
  auto values = m_values;
  int const num_blocks = U.size();
  int const num_eqs = U[0].extent(1);
  int const num_modes = U[0].extent(2);
  auto functor = [=] DGT_DEVICE (
      int const block,
      Vec3<int> const& cell_ijk) DGT_ALWAYS_INLINE
  {
    int const cell = m_cell_grid.index(cell_ijk);
    int const start = offsets_b[block];
    int const end = offsets_b[block+1];
    for (int msg_idx = start; msg_idx < end; ++msg_idx) {
      Subgrid3 const s = subgrids[msg_idx];
      if (!s.contains(cell_ijk)) return;
      int const packed_offset = offsets_c[msg_idx];
      int const local_cell = s.index(cell_ijk);
      int const packed_cell = packed_offset + local_cell;
      for (int eq = 0; eq < num_eqs; ++eq) {
        for (int mode = 0; mode < num_modes; ++mode) {
          values(packed_cell, eq, mode) = U[block](cell, eq, mode);
        }
      }
    }
  };
  for_each("pack", num_blocks, m_cell_grid, functor);
}

}
}
