#include <dgt_vec.hpp>

#include <gtest/gtest.h>

using namespace dgt;

TEST(vec, construct)
{
  Vec<real, 4> a;
  (void)a;
}

TEST(vec, mutable_accessor)
{
  Vec<real, 1> a = decltype(a)::zero();
  a[0] = 1.;
  EXPECT_EQ(a[0], 1.);
}

TEST(vec, immutable_accessor)
{
  Vec<real, 1> const a = decltype(a)::ones();
  EXPECT_EQ(a[0], 1.);
}

TEST(vec, zero)
{
  Vec<real, 2> const a = decltype(a)::zero();
  EXPECT_EQ(a[0], 0.);
  EXPECT_EQ(a[1], 0.);
}

TEST(vec, ones)
{
  Vec<real, 2> const a = decltype(a)::ones();
  EXPECT_EQ(a[0], 1.);
  EXPECT_EQ(a[1], 1.);
}

TEST(vec, unary_adition)
{
  Vec<real, 2> a = decltype(a)::zero();
  Vec<real, 2> const b = decltype(b)::ones();
  a += b;
  EXPECT_EQ(a[0], 1.);
  EXPECT_EQ(a[1], 1.);
}

TEST(vec, unary_subtraction)
{
  Vec<real, 2> a = decltype(a)::zero();
  Vec<real, 2> const b = decltype(b)::ones();
  a -= b;
  EXPECT_EQ(a[0], -1.);
  EXPECT_EQ(a[1], -1.);
}

TEST(vec, unary_multiplication)
{
  Vec<real, 2> a = decltype(a)::ones();
  real const b = 2.;
  a *= b;
  EXPECT_EQ(a[0], 2.);
  EXPECT_EQ(a[1], 2.);
}

TEST(vec, unary_division)
{
  Vec<real, 2> a = decltype(a)::ones();
  real const b = 2.;
  a /= b;
  EXPECT_EQ(a[0], 0.5);
  EXPECT_EQ(a[1], 0.5);
}

TEST(vec, binary_addition)
{
  Vec<real, 2> const a = decltype(a)::ones();
  Vec<real, 2> const b = decltype(b)::ones();
  Vec<real, 2> const c = a + b;
  EXPECT_EQ(c[0], 2.);
  EXPECT_EQ(c[1], 2.);
}

TEST(vec, binary_subtraction)
{
  Vec<real, 2> const a = decltype(a)::ones();
  Vec<real, 2> const b = decltype(b)::ones();
  Vec<real, 2> const c = a - b;
  EXPECT_EQ(c[0], 0.);
  EXPECT_EQ(c[1], 0.);
}

TEST(vec, binary_multiplication)
{
  Vec<real, 2> const a = decltype(a)::ones();
  real const b = 2.;
  Vec<real, 2> const c = a * b;
  EXPECT_EQ(c[0], 2.);
  EXPECT_EQ(c[1], 2.);
}

TEST(vec, binary_division)
{
  Vec<real, 2> const a = decltype(a)::ones();
  real const b = 2.;
  Vec<real, 2> const c = a / b;
  EXPECT_EQ(c[0], 0.5);
  EXPECT_EQ(c[1], 0.5);
}

TEST(vec, negation)
{
  Vec<real, 2> const a = decltype(a)::ones();
  Vec<real, 2> const b = -a;
  EXPECT_EQ(b[0], -1.);
  EXPECT_EQ(b[1], -1.);
}

TEST(vec, left_multiply)
{
  Vec<real, 2> const a = decltype(a)::ones();
  real const b = 2.;
  Vec<real, 2>  c = b * a;
  EXPECT_EQ(c[0], 2.);
  EXPECT_EQ(c[1], 2.);
}
