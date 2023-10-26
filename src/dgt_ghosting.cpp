#include <dgt_cartesian.hpp>
#include <dgt_for_each.hpp>
#include <dgt_ghosting.hpp>
#include <dgt_mesh.hpp>
#include <dgt_partitioning.hpp>

namespace dgt {

#if 0
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
#endif

static int count_messages(
    tree::ZLeaves const& owned_leaves,
    tree::Adjacencies const& owned_adjs)
{
  int num_messages = 0;
  for (auto const leaf_id : owned_leaves) {
    num_messages += owned_adjs.at(leaf_id).size();
  }
  return num_messages;
}

void Ghosting::build_helper_views(Mesh const& mesh)
{
  int msg = 0;
  int block = 0;
  m_num_cells = 0;
  int const dim = mesh.dim();
  auto const& leaves = mesh.owned_leaves();
  auto const& adjs = mesh.owned_adjacencies();
  int const num_blocks = int(leaves.size());
  m_msgs_per_block = DualView<int*>("Ghosting::msgs_per_block", num_blocks+1);
  m_subgrids[SEND] = DualView<Subgrid3*>("Ghosting::subgrids[SEND]", m_num_messages);
  m_subgrids[RECV] = DualView<Subgrid3*>("Ghosting::subgrids[RECV]", m_num_messages);
  m_buffer_offsets[SEND] = DualView<int*>("Ghosting::buffer_offsets[SEND]", m_num_messages);
  m_buffer_offsets[RECV] = DualView<int*>("Ghosting::buffer_offsets[RECV]", m_num_messages);
  for (auto const leaf_id : leaves) {
    m_msgs_per_block.h_view[block] = msg;
    for (auto const& adj : adjs.at(leaf_id)) {
      Subgrid3 owned_subgrid;
      Subgrid3 ghost_subgrid;
      Vec3<std::int8_t> const meta_ijk = adj.meta_ijk;
      if (adj.level_diff == 0) {
        owned_subgrid = get_cells(OWNED, m_cell_grid, meta_ijk);
        ghost_subgrid = get_cells(GHOST, m_cell_grid, meta_ijk);
      } else if (adj.level_diff == 1) {
        owned_subgrid = get_coarse_to_fine_cells(OWNED, m_cell_grid, meta_ijk);
        ghost_subgrid = get_coarse_to_fine_cells(GHOST, m_cell_grid, meta_ijk);
      } else if (adj.level_diff == -1) {
        owned_subgrid = get_fine_to_coarse_cells(OWNED, m_cell_grid, meta_ijk);
        ghost_subgrid = get_fine_to_coarse_cells(GHOST, m_cell_grid, meta_ijk);
      }
      m_subgrids[SEND].h_view[msg] = owned_subgrid;
      m_subgrids[RECV].h_view[msg] = ghost_subgrid;
      m_num_cells += generalize(dim, ghost_subgrid).size();
      msg++;
    }
    m_msgs_per_block.h_view[block+1] = msg;
    block++;
  }
  Kokkos::deep_copy(m_msgs_per_block.d_view, m_msgs_per_block.h_view);
  Kokkos::deep_copy(m_subgrids[SEND].d_view, m_subgrids[SEND].h_view);
  Kokkos::deep_copy(m_subgrids[RECV].d_view, m_subgrids[RECV].h_view);
  Kokkos::deep_copy(m_buffer_offsets[SEND].d_view, m_buffer_offsets[SEND].h_view);
  Kokkos::deep_copy(m_buffer_offsets[RECV].d_view, m_buffer_offsets[RECV].h_view);
}

void Ghosting::build_buffers(Mesh const& mesh)
{
  int max_eqs = 0;
  int const nmodes = mesh.basis().num_modes;
  for (auto const& d : mesh.get_modal_descriptors()) {
    max_eqs = std::max(max_eqs, d.num_comps);
  }
  m_buffer[SEND] = HostPinnedRightView<real***>(
      "Ghosting::buffer[SEND]", m_num_cells, max_eqs, nmodes);
  m_buffer[RECV] = HostPinnedRightView<real***>(
      "Ghosting::buffer[RECV]", m_num_cells, max_eqs, nmodes);
}

static int get_tag(
    int const local_block,
    int const nlocal_blocks,
    Vec3<std::int8_t> const& meta_ijk)
{
  Vec3<int> const ijk(meta_ijk.x(), meta_ijk.y(), meta_ijk.z());
  int const meta_idx = fine_meta_grid.index(ijk);
  int const tag = meta_idx + local_block * nlocal_blocks;
  return tag;
}

static int get_num_blocks(
    mpicpp::comm const* comm,
    tree::ZLeaves const& leaves,
    int const rank)
{
  int const nleaves = int(leaves.size());
  int const nranks = comm->size();
  return linear_partitioning::get_num_local(nleaves, nranks, rank);
}

void Ghosting::build_messages(Mesh const& mesh)
{
  int msg = 0;
  int block = 0;
  mpicpp::comm const* comm = mesh.comm();
  auto const& leaves = mesh.owned_leaves();
  auto const& adjs = mesh.owned_adjacencies();
  auto const& partitioning = mesh.partitioning();
  m_messages[SEND].resize(m_num_messages);
  m_messages[RECV].resize(m_num_messages);
  for (auto const leaf_id : leaves) {
    for (auto const& adj : adjs.at(leaf_id)) {
      tree::ID const adj_leaf_id = adj.neighbor;
      int const adj_rank = partitioning.at(adj_leaf_id).rank;
      int const adj_block = partitioning.at(adj_leaf_id).block;
      int const num_adj_rank_blocks = get_num_blocks(comm, leaves, adj_rank);
      m_messages[SEND][msg].rank = adj_rank;
      m_messages[RECV][msg].rank = adj_rank;
      //TODO: figure out what the tag is
      //TODO: figure out what the size of the message is
      //TODO: point to the appropriate place in the buffer view
      msg++;
    }
    block++;
  }
}

void Ghosting::build(Mesh const& mesh)
{
  auto const& leaves = mesh.owned_leaves();
  auto const& adjs = mesh.owned_adjacencies();
  m_cell_grid = mesh.cell_grid();
  m_num_messages = count_messages(leaves, adjs);
  build_helper_views(mesh);
  build_buffers(mesh);
  build_messages(mesh);
}

}
