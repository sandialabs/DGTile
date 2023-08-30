#include <dgt_subgrid3.hpp>

#include <gtest/gtest.h>

using namespace dgt;

TEST(subgrid3, construct_with_lower_upper)
{
  Vec3<int> const lower(1,1,1);
  Vec3<int> const upper(4,4,4);
  Subgrid3 const s(lower, upper);
  EXPECT_EQ(s.lower(), lower);
  EXPECT_EQ(s.upper(), upper);
}

TEST(subgrid3, construct_with_grid)
{
  Grid3 const grid(5,5,5);
  Subgrid3 const s(grid);
  EXPECT_EQ(s.lower(), Vec3<int>(0,0,0));
  EXPECT_EQ(s.upper(), Vec3<int>(5,5,5));
}

TEST(subgrid3, modify_lower)
{
  Subgrid3 s({5,5,5});
  s.lower().x() = 1;
  s.lower().y() = 1;
  s.lower().z() = 1;
  EXPECT_EQ(s.lower(), Vec3<int>(1,1,1));
}

TEST(subgrid3, modify_upper)
{
  Subgrid3 s({5,5,5});
  s.upper().x() = 4;
  s.upper().y() = 4;
  s.upper().z() = 4;
  EXPECT_EQ(s.upper(), Vec3<int>(4,4,4));
}

TEST(subgrid3, extents)
{
  Vec3<int> const lower(1,1,1);
  Vec3<int> const upper(5,5,5);
  Subgrid3 s(lower, upper);
  EXPECT_EQ(s.extents(), Vec3<int>(4,4,4));
}

TEST(subgrid3, size)
{
  EXPECT_EQ(Subgrid3({1,1,1}, {5,5,5}).size(), 64);
  EXPECT_EQ(Subgrid3({0,0,0}, {5,5,5}).size(), 125);
  EXPECT_EQ(Subgrid3({1,1,0}, {5,5,0}).size(), 0);
  EXPECT_EQ(generalize(2, Subgrid3({1,1,0}, {5,5,0})).size(), 16);
  EXPECT_EQ(generalize(2, Subgrid3({0,0,0}, {5,5,0})).size(), 25);
  EXPECT_EQ(Subgrid3({1,0,0}, {5,0,0}).size(), 0);
  EXPECT_EQ(generalize(1, Subgrid3({1,0,0}, {5,0,0})).size(), 4);
  EXPECT_EQ(generalize(1, Subgrid3({0,0,0}, {5,0,0})).size(), 5);
}

TEST(subgrid3, dimensionalize)
{
  EXPECT_EQ(dimensionalize(3, Subgrid3({1,1,1}, {8,8,8})), Subgrid3({1,1,1}, {8,8,8}));
  EXPECT_EQ(dimensionalize(2, Subgrid3({1,1,1}, {8,8,8})), Subgrid3({1,1,0}, {8,8,0}));
  EXPECT_EQ(dimensionalize(1, Subgrid3({1,1,1}, {8,8,8})), Subgrid3({1,0,0}, {8,0,0}));
  EXPECT_EQ(dimensionalize(3, Subgrid3({-1,-1,-1}, {8,8,8})), Subgrid3({-1,-1,-1}, {8,8,8}));
  EXPECT_EQ(dimensionalize(2, Subgrid3({-1,-1,-1}, {8,8,8})), Subgrid3({-1,-1,0}, {8,8,0}));
  EXPECT_EQ(dimensionalize(1, Subgrid3({-1,-1,-1}, {8,8,8})), Subgrid3({-1,0,0}, {8,0,0}));
}

TEST(subgrid3, generalize)
{
  EXPECT_EQ(generalize(3, Subgrid3({1,1,1}, {8,8,8})), Subgrid3({1,1,1}, {8,8,8}));
  EXPECT_EQ(generalize(2, Subgrid3({1,1,0}, {8,8,0})), Subgrid3({1,1,0}, {8,8,1}));
  EXPECT_EQ(generalize(1, Subgrid3({1,0,0}, {8,0,0})), Subgrid3({1,0,0}, {8,1,1}));
}
