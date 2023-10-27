#include <dgt_tree.hpp>

#include <gtest/gtest.h>

using namespace dgt;
using namespace dgt::tree;

namespace dgt {
namespace tree {

static ZLeaves get_example_amr(int const dim)
{
  ZLeaves z_leaves;
  if (dim == 1) z_leaves = {3,4,2};
  if (dim == 2) z_leaves = {5,6,9,10,2};
  if (dim == 3) z_leaves = {9,10,13,14,25,26,29,30,2};
  return z_leaves;
}

TEST(tree_amr, get_level_2D)
{
  int const dim = 2;
  ZLeaves const z_leaves = get_example_amr(dim);
  EXPECT_EQ(get_min_level(dim, z_leaves), 1);
  EXPECT_EQ(get_max_level(dim, z_leaves), 2);
}

//TODO: test base_pt in here..

TEST(tree_amr, get_adjacencies_2D)
{
  int const dim = 2;
  ZLeaves const z_leaves = get_example_amr(dim);
  Leaves leaves;
  for (auto const leaf_id : z_leaves) { leaves.insert(leaf_id); }
  Point const base_pt = get_base_point(dim, z_leaves);
  Vec3<bool> const periodic(false, false, false);
  Adjacencies const adjs = get_adjacencies(
      dim, z_leaves, leaves, base_pt, periodic);
}

}
}
