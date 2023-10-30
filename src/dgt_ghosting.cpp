#include "dgt_ghosting.hpp"
#include "dgt_mesh.hpp"

namespace dgt {

static int count_messages(Mesh const& mesh)
{
  int num_msg = 0;
  auto const& leaves = mesh.owned_leaves();
  auto const& adjs = mesh.owned_adjacencies();
  for (tree::ID const leaf_id : mesh.owned_leaves()) {
    num_msg += adjs.at(leaf_id).size();
  }
  return num_msg;
}

void Ghosting::build(Mesh const& mesh)
{
  m_num_messages = count_messages(mesh);
  (void)mesh;
}

}
