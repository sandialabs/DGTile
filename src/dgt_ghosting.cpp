#include "dgt_cartesian.hpp"
#include "dgt_for_each.hpp"
#include "dgt_ghosting.hpp"
#include "dgt_partitioning.hpp"
#include "dgt_tree.hpp"

namespace dgt {

static int count_messages(Mesh const& mesh)
{
  int num_msg = 0;
  auto const& leaves = mesh.owned_leaves();
  auto const& adjs = mesh.owned_adjacencies();
  for (tree::ID const global_id : leaves) {
    num_msg += adjs.at(global_id).size();
  }
  return num_msg;
}

static int get_max_eqs(Mesh const& mesh)
{
  int result = 0;
  for (auto const& f : mesh.get_modal_descriptors()) {
    result = std::max(result, f.num_comps);
  }
  return result;
}

template <class T>
void copy_to_device(T const& v)
{
  Kokkos::deep_copy(v.d_view, v.h_view);
}

void Ghosting::build_views(Mesh const& mesh)
{
  auto const& leaves = mesh.owned_leaves();
  auto const& adjs = mesh.owned_adjacencies();
  int const dim = mesh.dim();
  int const num_blocks = int(leaves.size());
  m_cell_grid = mesh.cell_grid();
  m_num_msgs = count_messages(mesh);
  m_max_eqs = get_max_eqs(mesh);
  m_num_modes = mesh.basis().num_modes;
  m_num_blocks = int(leaves.size());
  m_block_offsets = DualView<int*>("Ghosting::block_offsets", num_blocks+1);
  m_buffer_offsets = DualView<int*>("Ghosting::buffer_offsets", m_num_msgs+1);
  m_adjacencies = DualView<tree::Adjacent*>("Ghosting::adjacencies", m_num_msgs);
  m_subgrids[SEND] = DualView<Subgrid3*>("Ghosting::subgrids[SEND]", m_num_msgs);
  m_subgrids[RECV] = DualView<Subgrid3*>("Ghosting::subgrids[RECV]", m_num_msgs);
  int msg = 0;
  int num_buffer_cells = 0;
  m_block_offsets.h_view[0] = 0;
  m_buffer_offsets.h_view[0] = 0;
  for (int block = 0; block < m_num_blocks; ++block) {
    tree::ID const global_id = leaves[block];
    for (auto const& adj : adjs.at(global_id)) {
      int const axis = adj.axis;
      int const dir = adj.dir;
      int const kind = adj.kind;
      int const child = adj.which_child;
      Subgrid3 const owned_cells = generalize(dim, get_cells(
          m_cell_grid, OWNED, kind, axis, dir, child));
      Subgrid3 const ghost_cells = generalize(dim, get_cells(
          m_cell_grid, GHOST, kind, axis, dir, child));
      num_buffer_cells += ghost_cells.size();
      m_subgrids[SEND].h_view[msg] = owned_cells;
      m_subgrids[RECV].h_view[msg] = ghost_cells;
      m_adjacencies.h_view[msg] = adj;
      m_buffer_offsets.h_view[msg+1] = num_buffer_cells;
      msg++;
    }
    m_block_offsets.h_view[block+1] = msg;
  }
  m_buffers[SEND] = buffer_t(
      "Ghosting::buffers[SEND]", num_buffer_cells, m_max_eqs, m_num_modes);
  m_buffers[RECV] = buffer_t(
      "Ghosting::buffers[RECV]", num_buffer_cells, m_max_eqs, m_num_modes);
  copy_to_device(m_block_offsets);
  copy_to_device(m_buffer_offsets);
  copy_to_device(m_adjacencies);
  copy_to_device(m_subgrids[SEND]);
  copy_to_device(m_subgrids[RECV]);
}

static int get_tag(
    int const block,
    int const axis,
    int const dir,
    int const which_child)
{
  static constexpr int num_border = DIMENSIONS * DIRECTIONS;
  static constexpr int num_child = child_grid.size();
  int const border = axis * DIRECTIONS + dir;
  return (block * num_border + border) * num_child + which_child;
}

static int invert_child(int const axis, int const child)
{
  Vec3<int> const child_ijk = get_child_ijk(child);
  Vec3<int> ichild_ijk = child_ijk;
  ichild_ijk[axis] = (child_ijk[axis] == 0) ? 1 : 0;
  int const ichild = get_which_child(ichild_ijk);
  return ichild;
}

void Ghosting::build_messages(Mesh const& mesh)
{
  using namespace linear_partitioning;
  int const num_ranks = mesh.comm()->size();
  auto const& leaves = mesh.owned_leaves();
  auto const& adjs = mesh.owned_adjacencies();
  auto const inv_z_leaves = mesh.inv_z_leaves();
  m_messages[SEND].resize(m_num_msgs);
  m_messages[RECV].resize(m_num_msgs);
  int msg = 0;
  for (int block = 0; block < int(leaves.size()); ++block) {
    tree::ID const global_id = leaves[block];
    for (auto const& adj : adjs.at(global_id)) {
      tree::ID const adj_global_id = adj.neighbor;
      PartInfo const I = get_part_info(num_ranks, adj_global_id, inv_z_leaves);
      int const axis = adj.axis;
      int const dir = adj.dir;
      int const child = adj.which_child;
      int const idir = invert_dir(dir);
      int const ichild = invert_child(axis, child);
      int const buffer_start = m_buffer_offsets.h_view[msg];
      int const buffer_end = m_buffer_offsets.h_view[msg+1];
      int const num_cells = buffer_end - buffer_start;
      int const size = num_cells * m_max_eqs * m_num_modes;
      m_messages[SEND][msg].rank = I.rank;
      m_messages[SEND][msg].tag = get_tag(block, axis, dir, child);
      m_messages[SEND][msg].data = &(m_buffers[SEND](buffer_start, 0, 0));
      m_messages[SEND][msg].size = size;
      m_messages[RECV][msg].rank = I.rank;
      m_messages[RECV][msg].tag = get_tag(I.block, axis, idir, ichild);
      m_messages[RECV][msg].data = &(m_buffers[RECV](buffer_start, 0, 0));
      m_messages[RECV][msg].size = size;
      msg++;
    }
  }
}

void Ghosting::build(Mesh const& mesh)
{
  build_views(mesh);
  build_messages(mesh);
}

void Ghosting::begin_transfer(
    Field<real***> const& U,
    Basis<View> const& B,
    int const start_eq,
    int const end_eq)
{
  pack(U, B, start_eq, end_eq);
}

void Ghosting::pack(
    Field<real***> const& U_field,
    Basis<View> const& B,
    int const eq_start,
    int const eq_end)
{
  int const num_modes = m_num_modes;
  Grid3 const cell_grid = m_cell_grid;
  auto const U = U_field.get();
  auto const subgrids = m_subgrids[SEND].d_view;
  auto const block_offsets = m_block_offsets.d_view;
  auto const buffer_offsets = m_buffer_offsets.d_view;
  auto const adjacencies = m_adjacencies.d_view;
  auto buffer = m_buffers[SEND];
  auto functor = [=] DGT_DEVICE (
      int const block,
      Vec3<int> const& cell_ijk) DGT_ALWAYS_INLINE
  {
    int const msg_begin = block_offsets[block];
    int const msg_end = block_offsets[block+1];
    for (int msg = msg_begin; msg < msg_end; ++msg) {
      Subgrid3 const& subgrid = subgrids[msg];
      if (!subgrid.contains(cell_ijk)) continue;
      int const buffer_start = buffer_offsets[msg];
      int const buffer_idx = buffer_start + subgrid.index(cell_ijk);
      int const cell = cell_grid.index(cell_ijk);
      tree::Adjacent const& adj = adjacencies[msg];
      if (adj.kind == tree::EQUAL) {
        for (int eq = eq_start; eq < eq_end; ++eq) {
          for (int mode = 0; mode < num_modes; ++mode) {
            buffer(buffer_idx, eq, mode) = U[block](cell, eq, mode);
          }
        }
      }
      else {
        throw std::runtime_error("oops\n"); //TODO: do AMR
      }
    }
  };
  for_each("Ghosting::pack", m_num_blocks, m_cell_grid, functor);
  (void)B; //TODO: used for restricting
}

}
