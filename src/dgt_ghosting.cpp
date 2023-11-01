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

void Ghosting::build_messages(Mesh const& mesh)
{
  using namespace linear_partitioning;
  int const num_ranks = mesh.comm()->size();
  auto const& leaves = mesh.owned_leaves();
  auto const& adjs = mesh.owned_adjacencies();
  auto const& inv_zleaves = mesh.inv_z_leaves();
  int const num_msg = count_messages(mesh);
  m_messages[SEND].resize(num_msg);
  m_messages[RECV].resize(num_msg);
  int msg = 0;
  for (int block = 0; block < int(leaves.size()); ++block) {
    tree::ID const global_id = leaves[block];
    for (auto const& adj : adjs.at(global_id)) {
      tree::ID const adj_global_id = adj.neighbor;
      PartInfo const I = get_part_info(num_ranks, adj_global_id, inv_zleaves);
      int const tag = get_tag(I.block, adj.axis, adj.dir, adj.which_child);

      std::cout << tag << "\n";

      m_messages[SEND][msg].rank = I.rank;
      m_messages[SEND][msg].tag = tag;
      m_messages[SEND][msg].size = 0;
      m_messages[SEND][msg].data = nullptr;
      msg++;
    }
  }
}

void Ghosting::build(Mesh const& mesh)
{
  build_messages(mesh);
}


}
