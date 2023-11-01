#include "dgt_cartesian.hpp"
#include "dgt_ghosting.hpp"
#include "dgt_mesh.hpp"
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
  int const num_msgs = count_messages(mesh);
  Grid3 const grid = mesh.cell_grid();
  m_max_eqs = get_max_eqs(mesh);
  m_num_modes = mesh.basis().num_modes;
  m_msgs_per_block = DualView<int*>("Ghosting::msgs_per_block", num_blocks+1);
  m_buffer_offsets = DualView<int*>("Ghosting::buffer_offsets", num_msgs+1);
  m_adjacencies = DualView<tree::Adjacent*>("Ghosting::adjacencies", num_msgs);
  m_subgrids[SEND] = DualView<Subgrid3*>("Ghosting::subgrids[SEND]", num_msgs);
  m_subgrids[RECV] = DualView<Subgrid3*>("Ghosting::subgrids[RECV]", num_msgs);
  int msg = 0;
  int num_buffer_cells = 0;
  m_msgs_per_block.h_view[0] = 0;
  m_buffer_offsets.h_view[0] = 0;
  for (int block = 0; block < int(leaves.size()); ++block) {
    tree::ID const global_id = leaves[block];
    for (auto const& adj : adjs.at(global_id)) {
      int const axis = adj.axis;
      int const dir = adj.dir;
      int const kind = adj.kind;
      int const child = adj.which_child;
      Subgrid3 const owned_cells = get_cells(grid, OWNED, kind, axis, dir, child);
      Subgrid3 const ghost_cells = get_cells(grid, GHOST, kind, axis, dir, child);
      num_buffer_cells += generalize(dim, ghost_cells).size();
      m_subgrids[SEND].h_view[msg] = owned_cells;
      m_subgrids[RECV].h_view[msg] = ghost_cells;
      m_adjacencies.h_view[msg] = adj;
      m_buffer_offsets.h_view[msg+1] = num_buffer_cells;
      msg++;
    }
    m_msgs_per_block.h_view[block+1] = msg;
  }
  m_buffers[SEND] = buffer_t(
      "Ghosting::buffers[SEND]", num_buffer_cells, m_max_eqs, m_num_modes);
  m_buffers[RECV] = buffer_t(
      "Ghosting::buffers[RECV]", num_buffer_cells, m_max_eqs, m_num_modes);
  copy_to_device(m_msgs_per_block);
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

static int invert_child(int const child)
{
  return child;
}

void Ghosting::build_messages(Mesh const& mesh)
{
  using namespace linear_partitioning;
  int const num_ranks = mesh.comm()->size();
  auto const& leaves = mesh.owned_leaves();
  auto const& adjs = mesh.owned_adjacencies();
  auto const inv_z_leaves = mesh.inv_z_leaves();
  int const num_msg = m_adjacencies.h_view.size();
  m_messages[SEND].resize(num_msg);
  m_messages[RECV].resize(num_msg);
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
      int const ichild = invert_child(child);
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

}
