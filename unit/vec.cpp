#include <dgt_vec.hpp>

#include <gtest/gtest.h>

using namespace dgt;

TEST(vec, construct)
{
  Vec<real, 4> a = decltype(a)::zero();
  EXPECT_EQ(a[0], 0.);
}
