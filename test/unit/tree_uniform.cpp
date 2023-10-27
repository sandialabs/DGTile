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

TEST(tree_uniform, get_min_level_1D)
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

TEST(tree_uniform, get_min_level_2D)
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

TEST(tree_uniform, get_min_level_3D)
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

TEST(tree_uniform, get_max_level_1D)
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

TEST(tree_uniform, get_max_level_2D)
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

TEST(tree_uniform, get_max_level_3D)
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

TEST(tree_uniform, get_base_point_1D)
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

TEST(tree_uniform, get_base_point_2D)
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

TEST(tree_uniform, get_base_point_3D)
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

TEST(tree_uniform, get_domain_1D)
{
  int const dim = 1;
  ID const global_id = 1;
  Point const base_pt(0, {1,0,0});
  Box3<real> const domain({0.,0.,0.}, {1.,0.,0.});
  EXPECT_EQ(get_domain(dim, global_id, base_pt, domain).lower(), Vec3<real>(0.,0.,0.));
  EXPECT_EQ(get_domain(dim, global_id, base_pt, domain).upper(), Vec3<real>(0.5,0.,0.));
}

TEST(tree_uniform, get_domain_2D)
{
  int const dim = 2;
  ID const global_id = 1;
  Point const base_pt(0, {1,1,0});
  Box3<real> const domain({0.,0.,0.}, {1.,1.,0.});
  EXPECT_EQ(get_domain(dim, global_id, base_pt, domain).lower(), Vec3<real>(0.,0.,0.));
  EXPECT_EQ(get_domain(dim, global_id, base_pt, domain).upper(), Vec3<real>(0.5,0.5,0.));
}

TEST(tree_uniform, get_domain_3D)
{
  int const dim = 3;
  ID const global_id = 1;
  Point const base_pt(0, {1,1,1});
  Box3<real> const domain({0.,0.,0.}, {1.,1.,1.});
  EXPECT_EQ(get_domain(dim, global_id, base_pt, domain).lower(), Vec3<real>(0.,0.,0.));
  EXPECT_EQ(get_domain(dim, global_id, base_pt, domain).upper(), Vec3<real>(0.5,0.5,0.5));
}

TEST(tree_uniform, write_vtu_failure)
{
  int const dim = 3;
  Leaves const leaves = create(dim, {4,4,4});
  ZLeaves const z_leaves = order(dim, leaves);
  Box3<real> const domain({0.,0.,0.}, {1.,1.,1.});
  EXPECT_THROW(write_vtu(dim, "cant/do/that", z_leaves, domain), std::runtime_error);
}

TEST(tree_uniform, write_uniform_1D)
{
  int const dim = 1;
  Leaves const leaves = create(dim, {4,0,0});
  ZLeaves const z_leaves = order(dim, leaves);
  Box3<real> const domain({0.,0.,0.}, {1.,0.,0.});
  write_vtu(dim, "out_tree_uniform_1D", z_leaves, domain);
}

TEST(tree_uniform, write_uniform_2D)
{
  int const dim = 2;
  Leaves const leaves = create(dim, {4,4,0});
  ZLeaves const z_leaves = order(dim, leaves);
  Box3<real> const domain({0.,0.,0.}, {1.,1.,0.});
  write_vtu(dim, "out_tree_uniform_2D", z_leaves, domain);
}

TEST(tree_uniform, write_uniform_3D)
{
  int const dim = 3;
  Leaves const leaves = create(dim, {4,4,4});
  ZLeaves const z_leaves = order(dim, leaves);
  Box3<real> const domain({0.,0.,0.}, {1.,1.,1.});
  write_vtu(dim, "out_tree_uniform_2D", z_leaves, domain);
}
