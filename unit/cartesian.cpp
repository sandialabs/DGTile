#include <dgt_cartesian.hpp>

#include <gtest/gtest.h>

using namespace dgt;

TEST(cartesian, get_dir_sign)
{
  EXPECT_EQ(get_dir_sign(LEFT), -1);
  EXPECT_EQ(get_dir_sign(RIGHT), 1);
}

TEST(cartesian, invert_dir)
{
  EXPECT_EQ(invert_dir(LEFT), RIGHT);
  EXPECT_EQ(invert_dir(RIGHT), LEFT);
}

TEST(cartesian, get_cells_adj_face)
{
  EXPECT_EQ(get_cells_adj_face({2,2,2}, X, LEFT), vec3<int>(2,2,2));
  EXPECT_EQ(get_cells_adj_face({2,2,2}, X, RIGHT), vec3<int>(3,2,2));
  EXPECT_EQ(get_cells_adj_face({2,2,2}, Y, LEFT), vec3<int>(2,2,2));
  EXPECT_EQ(get_cells_adj_face({2,2,2}, Y, RIGHT), vec3<int>(2,3,2));
  EXPECT_EQ(get_cells_adj_face({2,2,2}, Z, LEFT), vec3<int>(2,2,2));
  EXPECT_EQ(get_cells_adj_face({2,2,2}, Z, RIGHT), vec3<int>(2,2,3));
}
