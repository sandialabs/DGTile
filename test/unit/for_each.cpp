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
  int sum =  0;
  auto functor = [&] (Vec3<int> const& ijk) {
    sum += ijk.x() + ijk.y() + ijk.z();
  };
  seq_for_each(Grid3(2,2,2), functor);
  EXPECT_EQ(sum, 12);
}
