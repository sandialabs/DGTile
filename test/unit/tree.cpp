#include <dgt_tree.hpp>

#include <gtest/gtest.h>

#include <dgt_print.hpp>

using namespace dgt;
using namespace dgt::tree;

TEST(tree, mark_definitions)
{
  EXPECT_EQ(DEREFINE, -1);
  EXPECT_EQ(REMAIN, 0);
  EXPECT_EQ(REFINE, 1);
}

TEST(tree, adjacency_definitions)
{
  EXPECT_EQ(FINE_TO_COARSE, -1);
  EXPECT_EQ(EQUAL, 0);
  EXPECT_EQ(COARSE_TO_FINE, 1);
}

TEST(tree, point_construction)
{
  Point const p(2, {1,2,3});
  EXPECT_EQ(p.level, 2);
  EXPECT_EQ(p.ijk, Vec3<int>(1,2,3));
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

TEST(tree, create_1D)
{
  Leaves const leaves = create(1, {3,0,0});
  EXPECT_EQ(leaves.size(), 3);
  EXPECT_EQ(leaves.count(ID(0)), 0);
  EXPECT_EQ(leaves.count(ID(1)), 0);
  EXPECT_EQ(leaves.count(ID(2)), 0);
  EXPECT_EQ(leaves.count(ID(3)), 1);
  EXPECT_EQ(leaves.count(ID(4)), 1);
  EXPECT_EQ(leaves.count(ID(5)), 1);
  EXPECT_EQ(leaves.count(ID(6)), 0);
}

TEST(tree, create_2D)
{
  Leaves const leaves = create(2, {3,2,0});
  EXPECT_EQ(leaves.size(), 6);
  EXPECT_EQ(leaves.count(ID(0)), 0);
  EXPECT_EQ(leaves.count(ID(1)), 0);
  EXPECT_EQ(leaves.count(ID(1)), 0);
  EXPECT_EQ(leaves.count(ID(5)), 1);
  EXPECT_EQ(leaves.count(ID(6)), 1);
  EXPECT_EQ(leaves.count(ID(7)), 1);
  EXPECT_EQ(leaves.count(ID(8)), 0);
  EXPECT_EQ(leaves.count(ID(9)), 1);
  EXPECT_EQ(leaves.count(ID(10)), 1);
  EXPECT_EQ(leaves.count(ID(11)), 1);
  EXPECT_EQ(leaves.count(ID(12)), 0);
}

TEST(tree, create_3D)
{
  Leaves const leaves = create(3, {3,2,1});
  EXPECT_EQ(leaves.size(), 6);
  EXPECT_EQ(leaves.count(ID(8)), 0);
  EXPECT_EQ(leaves.count(ID(9)), 1);
  EXPECT_EQ(leaves.count(ID(10)), 1);
  EXPECT_EQ(leaves.count(ID(11)), 1);
  EXPECT_EQ(leaves.count(ID(12)), 0);
  EXPECT_EQ(leaves.count(ID(13)), 1);
  EXPECT_EQ(leaves.count(ID(14)), 1);
  EXPECT_EQ(leaves.count(ID(15)), 1);
}
