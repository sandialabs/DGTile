#include <dgt_tree.hpp>

#include <gtest/gtest.h>

using namespace dgt;
using namespace dgt::tree;

TEST(tree, point_construction)
{
  Point const p(2, {1,2,3});
  EXPECT_EQ(p.level, 2);
  EXPECT_EQ(p.ijk, Vec3<int>(1,2,3));
}

TEST(tree, child_grid)
{
  EXPECT_EQ(child_grid.extents().x(), 2);
  EXPECT_EQ(child_grid.extents().y(), 2);
  EXPECT_EQ(child_grid.extents().z(), 2);
}

TEST(tree, adj_grid)
{
  EXPECT_EQ(adj_grid.extents().x(), 3);
  EXPECT_EQ(adj_grid.extents().y(), 3);
  EXPECT_EQ(adj_grid.extents().z(), 3);
}

TEST(tree, fine_adj_grid)
{
  EXPECT_EQ(fine_adj_grid.extents().x(), 4);
  EXPECT_EQ(fine_adj_grid.extents().y(), 4);
  EXPECT_EQ(fine_adj_grid.extents().z(), 4);
}

TEST(tree, get_level_offset_1D)
{
  int const dim = 1;
  EXPECT_EQ(get_level_offset(dim, 0), ID(0));
  EXPECT_EQ(get_level_offset(dim, 1), ID(1));
  EXPECT_EQ(get_level_offset(dim, 2), ID(3));
  EXPECT_EQ(get_level_offset(dim, 3), ID(7));
  EXPECT_EQ(get_level_offset(dim, 4), ID(15));
  EXPECT_EQ(get_level_offset(dim, 5), ID(31));
  EXPECT_EQ(get_level_offset(dim, 6), ID(63));
  EXPECT_EQ(get_level_offset(dim, 7), ID(127));
  EXPECT_EQ(get_level_offset(dim, 8), ID(255));
  EXPECT_EQ(get_level_offset(dim, 9), ID(511));
  EXPECT_EQ(get_level_offset(dim, 10), ID(1023));
  EXPECT_EQ(get_level_offset(dim, 11), ID(2047));
  EXPECT_EQ(get_level_offset(dim, 12), ID(4095));
}

TEST(tree, get_level_offset_2D)
{
  int const dim = 2;
  EXPECT_EQ(get_level_offset(dim, 0), ID(0));
  EXPECT_EQ(get_level_offset(dim, 1), ID(1));
  EXPECT_EQ(get_level_offset(dim, 2), ID(5));
  EXPECT_EQ(get_level_offset(dim, 3), ID(21));
  EXPECT_EQ(get_level_offset(dim, 4), ID(85));
  EXPECT_EQ(get_level_offset(dim, 5), ID(341));
  EXPECT_EQ(get_level_offset(dim, 6), ID(1365));
  EXPECT_EQ(get_level_offset(dim, 7), ID(5461));
  EXPECT_EQ(get_level_offset(dim, 8), ID(21845));
  EXPECT_EQ(get_level_offset(dim, 9), ID(87381));
  EXPECT_EQ(get_level_offset(dim, 10), ID(349525));
  EXPECT_EQ(get_level_offset(dim, 11), ID(1398101));
  EXPECT_EQ(get_level_offset(dim, 12), ID(5592405));
}

TEST(tree, get_level_offset_3D)
{
  int const dim = 3;
  EXPECT_EQ(get_level_offset(dim, 0), ID(0));
  EXPECT_EQ(get_level_offset(dim, 1), ID(1));
  EXPECT_EQ(get_level_offset(dim, 2), ID(9));
  EXPECT_EQ(get_level_offset(dim, 3), ID(73));
  EXPECT_EQ(get_level_offset(dim, 4), ID(585));
  EXPECT_EQ(get_level_offset(dim, 5), ID(4681));
  EXPECT_EQ(get_level_offset(dim, 6), ID(37449));
  EXPECT_EQ(get_level_offset(dim, 7), ID(299593));
  EXPECT_EQ(get_level_offset(dim, 8), ID(2396745));
  EXPECT_EQ(get_level_offset(dim, 9), ID(19173961));
  EXPECT_EQ(get_level_offset(dim, 10), ID(153391689));
  EXPECT_EQ(get_level_offset(dim, 11), ID(1227133513));
  EXPECT_EQ(get_level_offset(dim, 12), ID(9817068105));
}

TEST(tree, get_level_id_1D)
{
  int const dim = 1;
  EXPECT_EQ(get_level_id(dim, {0, {0,0,0}}), ID(0));
  EXPECT_EQ(get_level_id(dim, {1, {0,0,0}}), ID(0));
  EXPECT_EQ(get_level_id(dim, {1, {1,0,0}}), ID(1));
  EXPECT_EQ(get_level_id(dim, {2, {0,0,0}}), ID(0));
  EXPECT_EQ(get_level_id(dim, {2, {1,0,0}}), ID(1));
  EXPECT_EQ(get_level_id(dim, {2, {2,0,0}}), ID(2));
  EXPECT_EQ(get_level_id(dim, {2, {3,0,0}}), ID(3));
  EXPECT_EQ(get_level_id(dim, {3, {0,0,0}}), ID(0));
  EXPECT_EQ(get_level_id(dim, {3, {1,0,0}}), ID(1));
  EXPECT_EQ(get_level_id(dim, {3, {7,0,0}}), ID(7));
}

TEST(tree, get_level_id_2D)
{
  int const dim = 2;
  EXPECT_EQ(get_level_id(dim, {0, {0,0,0}}), ID(0));
  EXPECT_EQ(get_level_id(dim, {1, {0,0,0}}), ID(0));
  EXPECT_EQ(get_level_id(dim, {1, {1,0,0}}), ID(1));
  EXPECT_EQ(get_level_id(dim, {1, {0,1,0}}), ID(2));
  EXPECT_EQ(get_level_id(dim, {1, {1,1,0}}), ID(3));
  EXPECT_EQ(get_level_id(dim, {2, {0,0,0}}), ID(0));
  EXPECT_EQ(get_level_id(dim, {2, {1,0,0}}), ID(1));
  EXPECT_EQ(get_level_id(dim, {2, {2,3,0}}), ID(14));
  EXPECT_EQ(get_level_id(dim, {2, {3,3,0}}), ID(15));
}

TEST(tree, get_level_id_3D)
{
  int const dim = 3;
  EXPECT_EQ(get_level_id(dim, {0, {0,0,0}}), ID(0));
  EXPECT_EQ(get_level_id(dim, {1, {0,0,0}}), ID(0));
  EXPECT_EQ(get_level_id(dim, {1, {1,0,0}}), ID(1));
  EXPECT_EQ(get_level_id(dim, {1, {0,1,0}}), ID(2));
  EXPECT_EQ(get_level_id(dim, {1, {1,1,0}}), ID(3));
  EXPECT_EQ(get_level_id(dim, {1, {0,0,1}}), ID(4));
  EXPECT_EQ(get_level_id(dim, {1, {1,0,1}}), ID(5));
  EXPECT_EQ(get_level_id(dim, {1, {0,1,1}}), ID(6));
  EXPECT_EQ(get_level_id(dim, {1, {1,1,1}}), ID(7));
  EXPECT_EQ(get_level_id(dim, {2, {0,0,0}}), ID(0));
  EXPECT_EQ(get_level_id(dim, {2, {1,0,0}}), ID(1));
  EXPECT_EQ(get_level_id(dim, {2, {2,3,3}}), ID(62));
  EXPECT_EQ(get_level_id(dim, {2, {3,3,3}}), ID(63));
}
