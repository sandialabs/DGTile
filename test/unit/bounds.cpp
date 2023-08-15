#include <dgt_bounds.hpp>

#include <gtest/gtest.h>

using namespace dgt;

static Bounds<2> get_test_bounds()
{
  Bounds<2> b;
  b(0, MIN) = -1.;
  b(1, MIN) = -2.;
  b(0, MAX) =  1.;
  b(1, MAX) =  2.;
  return b;
}

static void test_const_accessor(Bounds<2> const& b)
{
  EXPECT_EQ(b(0, MIN), -1.);
  EXPECT_EQ(b(1, MIN), -2.);
  EXPECT_EQ(b(0, MAX),  1.);
  EXPECT_EQ(b(1, MAX),  2.);
}

TEST(bounds, initialize)
{
  get_test_bounds();
}

TEST(bounds, const_accessor)
{
  Bounds<2> const b = get_test_bounds();
  test_const_accessor(b);
}

TEST(bounds, clamp_within_range)
{
  Bounds<2> const b = get_test_bounds();
  EXPECT_EQ(b.clamp(0, -0.5), -0.5);
  EXPECT_EQ(b.clamp(0, -0.12), -0.12);
  EXPECT_EQ(b.clamp(0, 0.12), 0.12);
  EXPECT_EQ(b.clamp(0, 0.5), 0.5);
  EXPECT_EQ(b.clamp(1, -1.2), -1.2);
  EXPECT_EQ(b.clamp(1, -0.7), -0.7);
  EXPECT_EQ(b.clamp(1, 0.7), 0.7);
  EXPECT_EQ(b.clamp(1, 1.2), 1.2);
}

TEST(bounds, clamp_beyond_min)
{
  Bounds<2> const b = get_test_bounds();
  EXPECT_EQ(b.clamp(0, -2.2), -1.);
  EXPECT_EQ(b.clamp(0, -1.2), -1.);
  EXPECT_EQ(b.clamp(1, -2.7), -2.);
  EXPECT_EQ(b.clamp(1, -3.1), -2.);
}

TEST(bounds, clamp_beyond_max)
{
  Bounds<2> const b = get_test_bounds();
  EXPECT_EQ(b.clamp(0, 2.2), 1.);
  EXPECT_EQ(b.clamp(0, 1.2), 1.);
  EXPECT_EQ(b.clamp(1, 2.7), 2.);
  EXPECT_EQ(b.clamp(1, 3.1), 2.);
}
