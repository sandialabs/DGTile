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

TEST(tree, marking_definitions)
{
  EXPECT_EQ(DEREFINE, -1);
  EXPECT_EQ(REMAIN, 0);
  EXPECT_EQ(REFINE, 1);
}

TEST(tree, adjacency_definitions)
{
  EXPECT_EQ(COARSE_TO_FINE, -1);
  EXPECT_EQ(EQUAL, 0);
  EXPECT_EQ(FINE_TO_COARSE, 1);
}

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

TEST(tree, meta_grid)
{
  EXPECT_EQ(meta_grid.extents().x(), 3);
  EXPECT_EQ(meta_grid.extents().y(), 3);
  EXPECT_EQ(meta_grid.extents().z(), 3);
}

TEST(tree, fine_meta_grid)
{
  EXPECT_EQ(fine_meta_grid.extents().x(), 4);
  EXPECT_EQ(fine_meta_grid.extents().y(), 4);
  EXPECT_EQ(fine_meta_grid.extents().z(), 4);
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

TEST(tree, order_1D)
{
  int const dim = 1;
  Leaves const leaves = create(dim, {3,0,0});
  ZLeaves const z_leaves = order(dim, leaves);
  EXPECT_EQ(z_leaves.size(), 3);
  EXPECT_EQ(z_leaves[0], ID(3));
  EXPECT_EQ(z_leaves[1], ID(4));
  EXPECT_EQ(z_leaves[2], ID(5));
}

TEST(tree, order_2D)
{
  int const dim = 2;
  Leaves const leaves = create(dim, {3,2,0});
  ZLeaves const z_leaves = order(dim, leaves);
  EXPECT_EQ(z_leaves.size(), 6);
  EXPECT_EQ(z_leaves[0], ID(5));
  EXPECT_EQ(z_leaves[1], ID(6));
  EXPECT_EQ(z_leaves[2], ID(9));
  EXPECT_EQ(z_leaves[3], ID(10));
  EXPECT_EQ(z_leaves[4], ID(7));
  EXPECT_EQ(z_leaves[5], ID(11));
}

TEST(tree, order_3D)
{
  int const dim = 3;
  Leaves const leaves = create(dim, {3,2,1});
  ZLeaves const z_leaves = order(dim, leaves);
  ASSERT_EQ(z_leaves.size(), 6);
  EXPECT_EQ(z_leaves[0], ID(9));
  EXPECT_EQ(z_leaves[1], ID(10));
  EXPECT_EQ(z_leaves[2], ID(13));
  EXPECT_EQ(z_leaves[3], ID(14));
  EXPECT_EQ(z_leaves[4], ID(11));
  EXPECT_EQ(z_leaves[5], ID(15));
}

static Leaves get_example_refined(int const dim)
{
  Leaves const leaves = create(dim, {3,2,1});
  ZLeaves const z_leaves = order(dim, leaves);
  Marks marks(z_leaves.size(), REMAIN);
  marks[0] = REFINE;
  Leaves const new_leaves = modify(dim, z_leaves, marks);
  return new_leaves;
}

TEST(tree, modify_1D)
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

TEST(tree, modify_2D)
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

TEST(tree, modify_3D)
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

TEST(tree, get_min_level_uniform_1D)
{
  int const dim = 1;
  Leaves const leaves1 = create(dim, {1,0,0});
  Leaves const leaves2 = create(dim, {3,0,0});
  Leaves const leaves3 = create(dim, {22,0,0});
  ZLeaves const z_leaves1 = order(dim, leaves1);
  ZLeaves const z_leaves2 = order(dim, leaves2);
  ZLeaves const z_leaves3 = order(dim, leaves3);
  EXPECT_EQ(get_min_level(dim, leaves1), 0);
  EXPECT_EQ(get_min_level(dim, leaves2), 2);
  EXPECT_EQ(get_min_level(dim, leaves3), 5);
  EXPECT_EQ(get_min_level(dim, z_leaves1), 0);
  EXPECT_EQ(get_min_level(dim, z_leaves2), 2);
  EXPECT_EQ(get_min_level(dim, z_leaves3), 5);
}

TEST(tree, get_min_level_uniform_2D)
{
  int const dim = 2;
  Leaves const leaves1 = create(dim, {1,1,0});
  Leaves const leaves2 = create(dim, {3,2,0});
  Leaves const leaves3 = create(dim, {22,2,0});
  ZLeaves const z_leaves1 = order(dim, leaves1);
  ZLeaves const z_leaves2 = order(dim, leaves2);
  ZLeaves const z_leaves3 = order(dim, leaves3);
  EXPECT_EQ(get_min_level(dim, leaves1), 0);
  EXPECT_EQ(get_min_level(dim, leaves2), 2);
  EXPECT_EQ(get_min_level(dim, leaves3), 5);
  EXPECT_EQ(get_min_level(dim, z_leaves1), 0);
  EXPECT_EQ(get_min_level(dim, z_leaves2), 2);
  EXPECT_EQ(get_min_level(dim, z_leaves3), 5);
}

TEST(tree, get_min_level_uniform_3D)
{
  int const dim = 3;
  Leaves const leaves1 = create(dim, {1,1,1});
  Leaves const leaves2 = create(dim, {3,2,1});
  Leaves const leaves3 = create(dim, {22,2,1});
  ZLeaves const z_leaves1 = order(dim, leaves1);
  ZLeaves const z_leaves2 = order(dim, leaves2);
  ZLeaves const z_leaves3 = order(dim, leaves3);
  EXPECT_EQ(get_min_level(dim, leaves1), 0);
  EXPECT_EQ(get_min_level(dim, leaves2), 2);
  EXPECT_EQ(get_min_level(dim, leaves3), 5);
  EXPECT_EQ(get_min_level(dim, z_leaves1), 0);
  EXPECT_EQ(get_min_level(dim, z_leaves2), 2);
  EXPECT_EQ(get_min_level(dim, z_leaves3), 5);
}

TEST(tree, get_min_level_non_uniform_1D)
{
  int const dim = 1;
  Leaves const leaves = get_example_refined(dim);
  ZLeaves const z_leaves = order(dim, leaves);
  EXPECT_EQ(get_min_level(dim, leaves), 2);
  EXPECT_EQ(get_min_level(dim, z_leaves), 2);
}

TEST(tree, get_min_level_non_uniform_2D)
{
  int const dim = 2;
  Leaves const leaves = get_example_refined(dim);
  ZLeaves const z_leaves = order(dim, leaves);
  EXPECT_EQ(get_min_level(dim, leaves), 2);
  EXPECT_EQ(get_min_level(dim, z_leaves), 2);
}

TEST(tree, get_min_level_non_uniform_3D)
{
  int const dim = 3;
  Leaves const leaves = get_example_refined(dim);
  ZLeaves const z_leaves = order(dim, leaves);
  EXPECT_EQ(get_min_level(dim, leaves), 2);
  EXPECT_EQ(get_min_level(dim, z_leaves), 2);
}

TEST(tree, get_max_level_uniform_1D)
{
  int const dim = 1;
  Leaves const leaves1 = create(dim, {1,0,0});
  Leaves const leaves2 = create(dim, {3,0,0});
  Leaves const leaves3 = create(dim, {22,0,0});
  ZLeaves const z_leaves1 = order(dim, leaves1);
  ZLeaves const z_leaves2 = order(dim, leaves2);
  ZLeaves const z_leaves3 = order(dim, leaves3);
  EXPECT_EQ(get_max_level(dim, leaves1), 0);
  EXPECT_EQ(get_max_level(dim, leaves2), 2);
  EXPECT_EQ(get_max_level(dim, leaves3), 5);
  EXPECT_EQ(get_max_level(dim, z_leaves1), 0);
  EXPECT_EQ(get_max_level(dim, z_leaves2), 2);
  EXPECT_EQ(get_max_level(dim, z_leaves3), 5);
}

TEST(tree, get_max_level_uniform_2D)
{
  int const dim = 2;
  Leaves const leaves1 = create(dim, {1,1,0});
  Leaves const leaves2 = create(dim, {3,2,0});
  Leaves const leaves3 = create(dim, {22,2,0});
  ZLeaves const z_leaves1 = order(dim, leaves1);
  ZLeaves const z_leaves2 = order(dim, leaves2);
  ZLeaves const z_leaves3 = order(dim, leaves3);
  EXPECT_EQ(get_max_level(dim, leaves1), 0);
  EXPECT_EQ(get_max_level(dim, leaves2), 2);
  EXPECT_EQ(get_max_level(dim, leaves3), 5);
  EXPECT_EQ(get_max_level(dim, z_leaves1), 0);
  EXPECT_EQ(get_max_level(dim, z_leaves2), 2);
  EXPECT_EQ(get_max_level(dim, z_leaves3), 5);
}

TEST(tree, get_max_level_uniform_3D)
{
  int const dim = 3;
  Leaves const leaves1 = create(dim, {1,1,1});
  Leaves const leaves2 = create(dim, {3,2,1});
  Leaves const leaves3 = create(dim, {22,2,1});
  ZLeaves const z_leaves1 = order(dim, leaves1);
  ZLeaves const z_leaves2 = order(dim, leaves2);
  ZLeaves const z_leaves3 = order(dim, leaves3);
  EXPECT_EQ(get_max_level(dim, leaves1), 0);
  EXPECT_EQ(get_max_level(dim, leaves2), 2);
  EXPECT_EQ(get_max_level(dim, leaves3), 5);
  EXPECT_EQ(get_max_level(dim, z_leaves1), 0);
  EXPECT_EQ(get_max_level(dim, z_leaves2), 2);
  EXPECT_EQ(get_max_level(dim, z_leaves3), 5);
}

TEST(tree, get_max_level_non_uniform_1D)
{
  int const dim = 1;
  Leaves const leaves = get_example_refined(dim);
  ZLeaves const z_leaves = order(dim, leaves);
  EXPECT_EQ(get_max_level(dim, leaves), 3);
  EXPECT_EQ(get_max_level(dim, z_leaves), 3);
}

TEST(tree, get_max_level_non_uniform_2D)
{
  int const dim = 2;
  Leaves const leaves = get_example_refined(dim);
  ZLeaves const z_leaves = order(dim, leaves);
  EXPECT_EQ(get_max_level(dim, leaves), 3);
  EXPECT_EQ(get_max_level(dim, z_leaves), 3);
}

TEST(tree, get_max_level_non_uniform_3D)
{
  int const dim = 3;
  Leaves const leaves = get_example_refined(dim);
  ZLeaves const z_leaves = order(dim, leaves);
  EXPECT_EQ(get_max_level(dim, leaves), 3);
  EXPECT_EQ(get_max_level(dim, z_leaves), 3);
}

TEST(tree, get_base_point_uniform_1D)
{
  int const dim = 1;
  Leaves const leaves1 = create(dim, {1,0,0});
  Leaves const leaves2 = create(dim, {3,0,0});
  Leaves const leaves3 = create(dim, {22,0,0});
  ZLeaves const z_leaves1 = order(dim, leaves1);
  ZLeaves const z_leaves2 = order(dim, leaves2);
  ZLeaves const z_leaves3 = order(dim, leaves3);
  EXPECT_EQ(get_base_point(dim, leaves1), Point(0, {1,0,0}));
  EXPECT_EQ(get_base_point(dim, leaves2), Point(2, {3,0,0}));
  EXPECT_EQ(get_base_point(dim, leaves3), Point(5, {22,0,0}));
  EXPECT_EQ(get_base_point(dim, z_leaves1), Point(0, {1,0,0}));
  EXPECT_EQ(get_base_point(dim, z_leaves2), Point(2, {3,0,0}));
  EXPECT_EQ(get_base_point(dim, z_leaves3), Point(5, {22,0,0}));
}

TEST(tree, get_base_point_uniform_2D)
{
  int const dim = 2;
  Leaves const leaves1 = create(dim, {1,1,0});
  Leaves const leaves2 = create(dim, {3,2,0});
  Leaves const leaves3 = create(dim, {22,2,0});
  ZLeaves const z_leaves1 = order(dim, leaves1);
  ZLeaves const z_leaves2 = order(dim, leaves2);
  ZLeaves const z_leaves3 = order(dim, leaves3);
  EXPECT_EQ(get_base_point(dim, leaves1), Point(0, {1,1,0}));
  EXPECT_EQ(get_base_point(dim, leaves2), Point(2, {3,2,0}));
  EXPECT_EQ(get_base_point(dim, leaves3), Point(5, {22,2,0}));
  EXPECT_EQ(get_base_point(dim, z_leaves1), Point(0, {1,1,0}));
  EXPECT_EQ(get_base_point(dim, z_leaves2), Point(2, {3,2,0}));
  EXPECT_EQ(get_base_point(dim, z_leaves3), Point(5, {22,2,0}));
}

TEST(tree, get_base_point_uniform_3D)
{
  int const dim = 3;
  Leaves const leaves1 = create(dim, {1,1,1});
  Leaves const leaves2 = create(dim, {3,2,1});
  Leaves const leaves3 = create(dim, {22,2,1});
  ZLeaves const z_leaves1 = order(dim, leaves1);
  ZLeaves const z_leaves2 = order(dim, leaves2);
  ZLeaves const z_leaves3 = order(dim, leaves3);
  EXPECT_EQ(get_base_point(dim, leaves1), Point(0, {1,1,1}));
  EXPECT_EQ(get_base_point(dim, leaves2), Point(2, {3,2,1}));
  EXPECT_EQ(get_base_point(dim, leaves3), Point(5, {22,2,1}));
  EXPECT_EQ(get_base_point(dim, z_leaves1), Point(0, {1,1,1}));
  EXPECT_EQ(get_base_point(dim, z_leaves2), Point(2, {3,2,1}));
  EXPECT_EQ(get_base_point(dim, z_leaves3), Point(5, {22,2,1}));
}

TEST(tree, get_base_point_non_uniform_1D)
{
  int const dim = 1;
  Leaves const leaves = get_example_refined(dim);
  ZLeaves const z_leaves = order(dim, leaves);
  EXPECT_EQ(get_base_point(dim, leaves), Point(2, {3,0,0}));
  EXPECT_EQ(get_base_point(dim, z_leaves), Point(2, {3,0,0}));
}

TEST(tree, get_base_point_non_uniform_2D)
{
  int const dim = 2;
  Leaves const leaves = get_example_refined(dim);
  ZLeaves const z_leaves = order(dim, leaves);
  EXPECT_EQ(get_base_point(dim, leaves), Point(2, {3,2,0}));
  EXPECT_EQ(get_base_point(dim, z_leaves), Point(2, {3,2,0}));
}

TEST(tree, get_base_point_non_uniform_3D)
{
  int const dim = 3;
  Leaves const leaves = get_example_refined(dim);
  ZLeaves const z_leaves = order(dim, leaves);
  EXPECT_EQ(get_base_point(dim, leaves), Point(2, {3,2,1}));
  EXPECT_EQ(get_base_point(dim, z_leaves), Point(2, {3,2,1}));
}

TEST(tree, get_domain_1D)
{
  int const dim = 1;
  ID const global_id = 1;
  Point const base_pt(0, {1,0,0});
  Box3<real> const domain({0.,0.,0.}, {1.,0.,0.});
  EXPECT_EQ(get_domain(dim, global_id, base_pt, domain).lower(), Vec3<real>(0.,0.,0.));
  EXPECT_EQ(get_domain(dim, global_id, base_pt, domain).upper(), Vec3<real>(0.5,0.,0.));
}

TEST(tree, get_domain_2D)
{
  int const dim = 2;
  ID const global_id = 1;
  Point const base_pt(0, {1,1,0});
  Box3<real> const domain({0.,0.,0.}, {1.,1.,0.});
  EXPECT_EQ(get_domain(dim, global_id, base_pt, domain).lower(), Vec3<real>(0.,0.,0.));
  EXPECT_EQ(get_domain(dim, global_id, base_pt, domain).upper(), Vec3<real>(0.5,0.5,0.));
}

TEST(tree, get_domain_3D)
{
  int const dim = 3;
  ID const global_id = 1;
  Point const base_pt(0, {1,1,1});
  Box3<real> const domain({0.,0.,0.}, {1.,1.,1.});
  EXPECT_EQ(get_domain(dim, global_id, base_pt, domain).lower(), Vec3<real>(0.,0.,0.));
  EXPECT_EQ(get_domain(dim, global_id, base_pt, domain).upper(), Vec3<real>(0.5,0.5,0.5));
}

TEST(tree, write_vtu_failure)
{
  int const dim = 3;
  Leaves const leaves = create(dim, {4,4,4});
  ZLeaves const z_leaves = order(dim, leaves);
  Box3<real> const domain({0.,0.,0.}, {1.,1.,1.});
  EXPECT_THROW(write_vtu(dim, "cant/do/that", z_leaves, domain), std::runtime_error);
}

TEST(tree, write_uniform_1D)
{
  int const dim = 1;
  Leaves const leaves = create(dim, {4,0,0});
  ZLeaves const z_leaves = order(dim, leaves);
  Box3<real> const domain({0.,0.,0.}, {1.,0.,0.});
  write_vtu(dim, "out_tree_uniform_1D", z_leaves, domain);
}

TEST(tree, write_uniform_2D)
{
  int const dim = 2;
  Leaves const leaves = create(dim, {4,4,0});
  ZLeaves const z_leaves = order(dim, leaves);
  Box3<real> const domain({0.,0.,0.}, {1.,1.,0.});
  write_vtu(dim, "out_tree_uniform_2D", z_leaves, domain);
}

TEST(tree, write_uniform_3D)
{
  int const dim = 3;
  Leaves const leaves = create(dim, {4,4,4});
  ZLeaves const z_leaves = order(dim, leaves);
  Box3<real> const domain({0.,0.,0.}, {1.,1.,1.});
  write_vtu(dim, "out_tree_uniform_2D", z_leaves, domain);
}

TEST(vtu, write_non_uniform_1D)
{
  int const dim = 1;
  Leaves const leaves = get_example_refined(dim);
  ZLeaves const z_leaves = order(dim, leaves);
  Box3<real> const domain({0.,0.,0.}, {1.,0.,0.});
  write_vtu(dim, "unit_non_uniform_1D", z_leaves, domain);
}

TEST(vtu, write_non_uniform_2D)
{
  int const dim = 2;
  Leaves const leaves = get_example_refined(dim);
  ZLeaves const z_leaves = order(dim, leaves);
  Box3<real> const domain({0.,0.,0.}, {1.,1.,0.});
  write_vtu(dim, "unit_non_uniform_2D", z_leaves, domain);
}

TEST(vtu, write_non_uniform_3D)
{
  int const dim = 3;
  Leaves const leaves = get_example_refined(dim);
  ZLeaves const z_leaves = order(dim, leaves);
  Box3<real> const domain({0.,0.,0.}, {1.,1.,1.});
  write_vtu(dim, "unit_non_uniform_3D", z_leaves, domain);
}
