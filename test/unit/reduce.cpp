#include <dgt_defines.hpp>
#include <dgt_reduce.hpp>
#include <dgt_subgrid3.hpp>

#include <gtest/gtest.h>

using namespace dgt;

static void test_4D_reduce()
{
  int const nblocks = 4;
  Grid3 const cell_grid(3,3,3);
  real result = 0.;
  auto functor = [=] DGT_HOST_DEVICE (
      int const block,
      Vec3<int> const& cell_ijk,
      real& result) {
    (void)block;
    (void)cell_ijk;
    result += 1.;
  };
  reduce_for_each("4D_reduce", nblocks, cell_grid, functor, result);
  EXPECT_EQ(result, 108);
}

TEST(reduce, whatever)
{
  test_4D_reduce();
}
