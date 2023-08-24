#include <dgt_tree.hpp>

#include <gtest/gtest.h>

using namespace dgt;
using namespace dgt::tree;

namespace dgt {
namespace tree {

static bool operator==(Point const& a, Point const& b)
{
  return
    (a.level == b.level) &&
    (a.ijk == b.ijk);
}

}
}

static Leaves refine_zleaf(
    int const dim,
    Leaves const& leaves,
    int const index)
{
  ZLeaves const z_leaves = order(dim, leaves);
  Marks marks(z_leaves.size(), REMAIN);
  marks[index] = REFINE;
  return modify(dim, z_leaves, marks);
}

static Leaves get_example_refined(int const dim)
{
  Leaves const leaves = create(dim, {3,2,1});
  return refine_zleaf(dim, leaves, 0);
}

static Leaves get_example_refined2(int const dim)
{
  Leaves const leaves = get_example_refined(dim);
  return refine_zleaf(dim, leaves, 3);
}

TEST(tree_amr, modify_1D)
{
  int const dim = 1;
  Leaves const leaves = get_example_refined(dim);
  EXPECT_EQ(leaves.size(), 4);
  EXPECT_EQ(leaves.count(ID(3)), 0);
  EXPECT_EQ(leaves.count(ID(4)), 1);
  EXPECT_EQ(leaves.count(ID(5)), 1);
  EXPECT_EQ(leaves.count(ID(6)), 0);
  EXPECT_EQ(leaves.count(ID(7)), 1);
  EXPECT_EQ(leaves.count(ID(8)), 1);
  EXPECT_EQ(leaves.count(ID(9)), 0);
}

TEST(tree_amr, modify_2D)
{
  int const dim = 2;
  Leaves const leaves = get_example_refined(dim);
  EXPECT_EQ(leaves.size(), 9);
  EXPECT_EQ(leaves.count(ID(5)), 0);
  EXPECT_EQ(leaves.count(ID(6)), 1);
  EXPECT_EQ(leaves.count(ID(7)), 1);
  EXPECT_EQ(leaves.count(ID(8)), 0);
  EXPECT_EQ(leaves.count(ID(9)), 1);
  EXPECT_EQ(leaves.count(ID(10)), 1);
  EXPECT_EQ(leaves.count(ID(11)), 1);
  EXPECT_EQ(leaves.count(ID(12)), 0);
  EXPECT_EQ(leaves.count(ID(20)), 0);
  EXPECT_EQ(leaves.count(ID(21)), 1);
  EXPECT_EQ(leaves.count(ID(22)), 1);
  EXPECT_EQ(leaves.count(ID(29)), 1);
  EXPECT_EQ(leaves.count(ID(30)), 1);
  EXPECT_EQ(leaves.count(ID(31)), 0);
}

TEST(tree_amr, modify_3D)
{
  int const dim = 3;
  Leaves const leaves = get_example_refined(dim);
  EXPECT_EQ(leaves.size(), 13);
  EXPECT_EQ(leaves.count(ID(10)), 1);
  EXPECT_EQ(leaves.count(ID(11)), 1);
  EXPECT_EQ(leaves.count(ID(13)), 1);
  EXPECT_EQ(leaves.count(ID(14)), 1);
  EXPECT_EQ(leaves.count(ID(15)), 1);
  EXPECT_EQ(leaves.count(ID(73)), 1);
  EXPECT_EQ(leaves.count(ID(74)), 1);
  EXPECT_EQ(leaves.count(ID(81)), 1);
  EXPECT_EQ(leaves.count(ID(82)), 1);
  EXPECT_EQ(leaves.count(ID(137)), 1);
  EXPECT_EQ(leaves.count(ID(138)), 1);
  EXPECT_EQ(leaves.count(ID(145)), 1);
  EXPECT_EQ(leaves.count(ID(146)), 1);
}

TEST(tree_amr, get_min_level_1D)
{
  int const dim = 1;
  Leaves const leaves = get_example_refined(dim);
  ZLeaves const z_leaves = order(dim, leaves);
  EXPECT_EQ(get_min_level(dim, leaves), 2);
  EXPECT_EQ(get_min_level(dim, z_leaves), 2);
}

TEST(tree_amr, get_min_level_2D)
{
  int const dim = 2;
  Leaves const leaves = get_example_refined(dim);
  ZLeaves const z_leaves = order(dim, leaves);
  EXPECT_EQ(get_min_level(dim, leaves), 2);
  EXPECT_EQ(get_min_level(dim, z_leaves), 2);
}

TEST(tree_amr, get_min_level_3D)
{
  int const dim = 3;
  Leaves const leaves = get_example_refined(dim);
  ZLeaves const z_leaves = order(dim, leaves);
  EXPECT_EQ(get_min_level(dim, leaves), 2);
  EXPECT_EQ(get_min_level(dim, z_leaves), 2);
}

TEST(tree_amr, get_max_level_1D)
{
  int const dim = 1;
  Leaves const leaves = get_example_refined(dim);
  ZLeaves const z_leaves = order(dim, leaves);
  EXPECT_EQ(get_max_level(dim, leaves), 3);
  EXPECT_EQ(get_max_level(dim, z_leaves), 3);
}

TEST(tree_amr, get_max_level_2D)
{
  int const dim = 2;
  Leaves const leaves = get_example_refined(dim);
  ZLeaves const z_leaves = order(dim, leaves);
  EXPECT_EQ(get_max_level(dim, leaves), 3);
  EXPECT_EQ(get_max_level(dim, z_leaves), 3);
}

TEST(tree_amr, get_max_level_3D)
{
  int const dim = 3;
  Leaves const leaves = get_example_refined(dim);
  ZLeaves const z_leaves = order(dim, leaves);
  EXPECT_EQ(get_max_level(dim, leaves), 3);
  EXPECT_EQ(get_max_level(dim, z_leaves), 3);
}

TEST(tree_amr, get_base_point_1D)
{
  int const dim = 1;
  Leaves const leaves = get_example_refined(dim);
  ZLeaves const z_leaves = order(dim, leaves);
  EXPECT_EQ(get_base_point(dim, leaves), Point(2, {3,0,0}));
  EXPECT_EQ(get_base_point(dim, z_leaves), Point(2, {3,0,0}));
}

TEST(tree_amr, get_base_point_2D)
{
  int const dim = 2;
  Leaves const leaves = get_example_refined(dim);
  ZLeaves const z_leaves = order(dim, leaves);
  EXPECT_EQ(get_base_point(dim, leaves), Point(2, {3,2,0}));
  EXPECT_EQ(get_base_point(dim, z_leaves), Point(2, {3,2,0}));
}

TEST(tree_amr, get_base_point_3D)
{
  int const dim = 3;
  Leaves const leaves = get_example_refined(dim);
  ZLeaves const z_leaves = order(dim, leaves);
  EXPECT_EQ(get_base_point(dim, leaves), Point(2, {3,2,1}));
  EXPECT_EQ(get_base_point(dim, z_leaves), Point(2, {3,2,1}));
}

TEST(tree_amr, get_adjacencies_failure)
{
  int const dim = 2;
  Leaves const leaves = get_example_refined2(dim);
  ZLeaves const z_leaves = order(dim, leaves);
  Point const base_pt = get_base_point(dim, z_leaves);
  Periodic const periodic(true, false, false);
  EXPECT_THROW((void)get_adjacencies(dim, leaves, base_pt, periodic), std::runtime_error);
}

TEST(tree_amr, balance_2D)
{
  int const dim = 2;
  Periodic const periodic(true, false, false);
  Box3<real> const domain({0.,0.,0.}, {3.,2.,0.});
  Leaves const leaves = get_example_refined2(dim);
  ZLeaves const z_leaves = order(dim, leaves);
  Point const base_pt = get_base_point(dim, leaves);
  Leaves const balanced_leaves = balance(dim, leaves, base_pt, periodic);
  ZLeaves const balanced_z_leaves = order(dim, balanced_leaves);

  write_vtu(dim, "debug_before", z_leaves, domain);
  write_vtu(dim, "debug_after", balanced_z_leaves, domain);

}

TEST(tree_amr, write_non_uniform_1D)
{
  int const dim = 1;
  Leaves const leaves = get_example_refined(dim);
  ZLeaves const z_leaves = order(dim, leaves);
  Box3<real> const domain({0.,0.,0.}, {3.,0.,0.});
  write_vtu(dim, "out_tree_amr_1D", z_leaves, domain);
}

TEST(tree_amr, write_non_uniform_2D)
{
  int const dim = 2;
  Leaves const leaves = get_example_refined(dim);
  ZLeaves const z_leaves = order(dim, leaves);
  Box3<real> const domain({0.,0.,0.}, {3.,2.,0.});
  write_vtu(dim, "out_tree_amr_2D", z_leaves, domain);
}

TEST(tree_amr, write_non_uniform_3D)
{
  int const dim = 3;
  Leaves const leaves = get_example_refined(dim);
  ZLeaves const z_leaves = order(dim, leaves);
  Box3<real> const domain({0.,0.,0.}, {3.,2.,1.});
  write_vtu(dim, "out_tree_amr_3D", z_leaves, domain);
}
