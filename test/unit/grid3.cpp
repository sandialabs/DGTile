#include <dgt_grid3.hpp>

#include <gtest/gtest.h>

#include <dgt_print.hpp> // debug

using namespace dgt;

TEST(grid3, construct_with_values)
{
  Grid3 a(1,2,3);
  EXPECT_EQ(a.extents().x(), 1);
  EXPECT_EQ(a.extents().y(), 2);
  EXPECT_EQ(a.extents().z(), 3);
}

TEST(grid3, construct_with_vec3)
{
  Grid3 a(Vec3<int>(1,2,3));
  EXPECT_EQ(a.extents().x(), 1);
  EXPECT_EQ(a.extents().y(), 2);
  EXPECT_EQ(a.extents().z(), 3);
}

TEST(grid3, extents)
{
  Grid3 a(Vec3<int>(1,2,3));
  EXPECT_EQ(a.extents(), Vec3<int>(1,2,3));
}

TEST(grid3, size)
{
  EXPECT_EQ(Grid3(1,2,3).size(), 6);
  EXPECT_EQ(Grid3(5,2,0).size(), 0);
  EXPECT_EQ(Grid3(27,0,0).size(), 0);
}

TEST(grid3, index)
{
  Grid3 const a(2,0,0);
  Grid3 const b(2,2,0);
  Grid3 const c(2,2,2);
  EXPECT_EQ(a.index({0,0,0}), 0);
  EXPECT_EQ(a.index({1,0,0}), 1);
  EXPECT_EQ(b.index({0,0,0}), 0);
  EXPECT_EQ(b.index({1,0,0}), 1);
  EXPECT_EQ(b.index({0,1,0}), 2);
  EXPECT_EQ(b.index({1,1,0}), 3);
  EXPECT_EQ(c.index({0,0,0}), 0);
  EXPECT_EQ(c.index({1,0,0}), 1);
  EXPECT_EQ(c.index({0,1,0}), 2);
  EXPECT_EQ(c.index({1,1,0}), 3);
  EXPECT_EQ(c.index({0,0,1}), 4);
  EXPECT_EQ(c.index({1,0,1}), 5);
  EXPECT_EQ(c.index({0,1,1}), 6);
  EXPECT_EQ(c.index({1,1,1}), 7);
}

TEST(grid3, ijk)
{
  Grid3 const a = generalize(1, Grid3(2,0,0));
  Grid3 const b = generalize(2, Grid3(2,2,0));
  Grid3 const c = generalize(3, Grid3(2,2,2));
  EXPECT_EQ(a.ijk(0), Vec3<int>(0,0,0));
  EXPECT_EQ(a.ijk(1), Vec3<int>(1,0,0));
  EXPECT_EQ(b.ijk(0), Vec3<int>(0,0,0));
  EXPECT_EQ(b.ijk(1), Vec3<int>(1,0,0));
  EXPECT_EQ(b.ijk(2), Vec3<int>(0,1,0));
  EXPECT_EQ(b.ijk(3), Vec3<int>(1,1,0));
  EXPECT_EQ(c.ijk(0), Vec3<int>(0,0,0));
  EXPECT_EQ(c.ijk(1), Vec3<int>(1,0,0));
  EXPECT_EQ(c.ijk(2), Vec3<int>(0,1,0));
  EXPECT_EQ(c.ijk(3), Vec3<int>(1,1,0));
  EXPECT_EQ(c.ijk(4), Vec3<int>(0,0,1));
  EXPECT_EQ(c.ijk(5), Vec3<int>(1,0,1));
  EXPECT_EQ(c.ijk(6), Vec3<int>(0,1,1));
  EXPECT_EQ(c.ijk(7), Vec3<int>(1,1,1));
}

TEST(grid3, infer_dimension)
{
  EXPECT_EQ(infer_dimension(Grid3(3,0,0)), 1);
  EXPECT_EQ(infer_dimension(Grid3(3,4,0)), 2);
  EXPECT_EQ(infer_dimension(Grid3(3,4,8)), 3);
  EXPECT_EQ(infer_dimension(Grid3(0,4,8)), -1);
  EXPECT_EQ(infer_dimension(Grid3(4,0,8)), -1);
  EXPECT_EQ(infer_dimension(Grid3(-1,0,0)), -1);
  EXPECT_EQ(infer_dimension(Grid3(-1,-2,0)), -1);
  EXPECT_EQ(infer_dimension(Grid3(-1,-2,-3)), -1);
}

TEST(grid3, dimensionalize)
{
  Grid3 const a(8,8,8);
  EXPECT_EQ(dimensionalize(1, a), Grid3(8,0,0));
  EXPECT_EQ(dimensionalize(2, a), Grid3(8,8,0));
  EXPECT_EQ(dimensionalize(3, a), Grid3(8,8,8));
}

TEST(grid3, generalize)
{
  Grid3 const a(8,0,0);
  Grid3 const b(8,8,0);
  Grid3 const c(8,8,8);
  EXPECT_EQ(generalize(1, a), Grid3(8,1,1));
  EXPECT_EQ(generalize(2, b), Grid3(8,8,1));
  EXPECT_EQ(generalize(3, c), Grid3(8,8,8));
}
