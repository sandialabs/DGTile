#include <dgt_grid3.hpp>

#include <gtest/gtest.h>

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

TEST(grid3, dimension)
{
  EXPECT_EQ(Grid3(1,2,3).dimension(), 3);
  EXPECT_EQ(Grid3(5,2,0).dimension(), 2);
  EXPECT_EQ(Grid3(27,0,0).dimension(), 1);
  EXPECT_EQ(Grid3(0,1,8).dimension(), -1);
  EXPECT_EQ(Grid3(1,0,3).dimension(), -1);
}

TEST(grid3, size)
{
  EXPECT_EQ(Grid3(1,2,3).size(), 6);
  EXPECT_EQ(Grid3(5,2,0).size(), 10);
  EXPECT_EQ(Grid3(27,0,0).size(), 27);
  EXPECT_EQ(Grid3(0,1,8).size(), -1);
  EXPECT_EQ(Grid3(1,0,3).size(), -1);
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
  Grid3 const a(2,0,0);
  Grid3 const b(2,2,0);
  Grid3 const c(2,2,2);
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
