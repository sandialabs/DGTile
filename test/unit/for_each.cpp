#include <dgt_for_each.hpp>

#include <gtest/gtest.h>

using namespace dgt;

TEST(for_each, sequenced_invalid_grid)
{
  int sum =  0;
  auto functor = [&] (Vec3<int> const& ijk) {
    sum += ijk.x() + ijk.y() + ijk.z();
  };
  seq_for_each(Grid3(2,0,2), functor);
  EXPECT_EQ(sum, 0);
}

TEST(for_each, sequenced_1D_grid)
{
  int sum =  0;
  auto functor = [&] (Vec3<int> const& ijk) {
    EXPECT_EQ(ijk.y(), 0);
    EXPECT_EQ(ijk.z(), 0);
    sum += ijk.x() + ijk.y() + ijk.z();
  };
  seq_for_each(Grid3(2,0,0), functor);
  EXPECT_EQ(sum, 1);
}

TEST(for_each, sequenced_2D_grid)
{
  int sum = 0;
  auto functor = [&] (Vec3<int> const& ijk) {
    EXPECT_EQ(ijk.z(), 0);
    sum += ijk.x() + ijk.y() + ijk.z();
  };
  seq_for_each(Grid3(2,2,0), functor);
  EXPECT_EQ(sum, 4);
}

TEST(for_each, sequenced_3D_grid)
{
  int sum = 0;
  auto functor = [&] (Vec3<int> const& ijk) {
    sum += ijk.x() + ijk.y() + ijk.z();
  };
  seq_for_each(Grid3(2,2,2), functor);
  EXPECT_EQ(sum, 12);
}

static Subgrid3 offset_grid = {{-1,-1,-1}, {2,2,2}};

TEST(for_each, sequenced_offset_grid_1D)
{
  int sum = 0;
  auto functor = [&] (Vec3<int> const& ijk) {
    (void)ijk;
    sum += 1;
  };
  seq_for_each(dimensionalize(1, offset_grid), functor);
  EXPECT_EQ(sum, 3);
}

TEST(for_each, sequenced_offset_grid_2D)
{
  int sum = 0;
  auto functor = [&] (Vec3<int> const& ijk) {
    (void)ijk;
    sum += 1;
  };
  seq_for_each(dimensionalize(2, offset_grid), functor);
  EXPECT_EQ(sum, 9);
}

TEST(for_each, sequenced_offset_grid_3D)
{
  int sum = 0;
  auto functor = [&] (Vec3<int> const& ijk) {
    (void)ijk;
    sum += 1;
  };
  seq_for_each(dimensionalize(3, offset_grid), functor);
  EXPECT_EQ(sum, 27);
}

//TODO: unit test 4D for each
