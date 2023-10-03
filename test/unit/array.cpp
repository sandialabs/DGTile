#include <dgt_array.hpp>

#include <Kokkos_Core.hpp>

#include <gtest/gtest.h>

using namespace dgt;

void test_on_device()
{
  Array<int, 2> a;
  a[0] = 1.;
  a[1] = 2.;
  auto functor = [=] DGT_DEVICE (int const i) {
    Array<int, 2> b;
    b[0] = a[0];
    b[1] = a[1];
    (void)i;
  };
  Kokkos::parallel_for("unit_array", 1, functor);
}

TEST(array, construct)
{
  Array<real, 1> a;
  a[0] = 0.;
}

TEST(array, access)
{
  Array<real, 2> a;
  a[0] = 1.;
  a[1] = 7.;
  EXPECT_EQ(a[0], 1.);
  EXPECT_EQ(a[1], 7.);
}

TEST(array, on_device)
{
  test_on_device();
}
