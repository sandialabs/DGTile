#include <dgt_defines.hpp>
#include <dgt_view.hpp>

#include <gtest/gtest.h>

using namespace dgt;

template <class T>
void test_sum()
{
  View<real*> v("v", 10);
  Kokkos::deep_copy(v, T(1));
  T const result = sum(v);
  EXPECT_EQ(result, T(10));
}

TEST(view, sum) {
  test_sum<int>();
  test_sum<real>();
}
