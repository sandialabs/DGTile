#include <dgt_defines.hpp>

#include <gtest/gtest.h>

using namespace dgt;

TEST(defines, spatial)
{
  EXPECT_EQ(X, 0);
  EXPECT_EQ(Y, 1);
  EXPECT_EQ(Z, 2);
  EXPECT_EQ(DIMENSIONS, 3);
}

TEST(defines, directions)
{
  EXPECT_EQ(MIN, 0);
  EXPECT_EQ(MAX, 1);
  EXPECT_EQ(LEFT, 0);
  EXPECT_EQ(RIGHT, 1);
  EXPECT_EQ(DIRECTIONS, 2);
}

TEST(defines, communication)
{
  EXPECT_EQ(SEND, 0);
  EXPECT_EQ(RECV, 1);
  EXPECT_EQ(OWNED, 0);
  EXPECT_EQ(GHOST, 1);
}

TEST(defines, errors)
{
  EXPECT_EQ(errors, 0);
}
