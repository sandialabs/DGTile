#include "gtest/gtest.h"

#include "dgt_basis.hpp"

TEST(basis, init_tensor) {
  for (int dim = 1; dim <= 3; ++dim) {
    for (int p = 0; p <= 2; ++p) {
      dgt::Basis basis;
      basis.init(dim, p, true);
    }
  }
}

TEST(basis, init_non_tensor) {
  for (int dim = 1; dim <= 3; ++dim) {
    for (int p = 0; p <= 2; ++p) {
      dgt::Basis basis;
      basis.init(dim, p, false);
    }
  }
}

template <class ViewT>
void test_sum(ViewT v, double exact, double tol = 1.e-14) {
  double const val = dgt::sum(p3a::execution::par, v);
  EXPECT_NEAR(exact, val, tol);
}

TEST(basis, pt_sums_tensor) {
  for (int dim = 1; dim <= 3; ++dim) {
    for (int p = 0; p <= 2; ++p) {
      dgt::Basis basis;
      basis.init(dim, p, true);
      test_sum(basis.wt_intr, dgt::ipow(2, dim));
      test_sum(basis.wt_side, dgt::ipow(2, dim-1));
      test_sum(basis.wt_fine, dgt::ipow(2, dim));
      test_sum(basis.pt_intr, 0.);
      test_sum(basis.pt_side, 0.);
      test_sum(basis.pt_child_intr, 0.);
      test_sum(basis.pt_fine, 0.);
      test_sum(basis.pt_viz, 0.);
      test_sum(basis.pt_corner, 0.);
    }
  }
}

TEST(basis, pt_sums_non_tensor) {
  for (int dim = 1; dim <= 3; ++dim) {
    for (int p = 0; p <= 2; ++p) {
      dgt::Basis basis;
      basis.init(dim, p, false);
      test_sum(basis.wt_intr, dgt::ipow(2, dim));
      test_sum(basis.wt_side, dgt::ipow(2, dim-1));
      test_sum(basis.wt_fine, dgt::ipow(2, dim));
      test_sum(basis.pt_intr, 0.);
      test_sum(basis.pt_side, 0.);
      test_sum(basis.pt_child_intr, 0.);
      test_sum(basis.pt_fine, 0.);
      test_sum(basis.pt_viz, 0.);
      test_sum(basis.pt_corner, 0.);
    }
  }
}

TEST(basis, copying) {
  dgt::Basis a;
  dgt::Basis b;
  a.init(3,1,true);
  a = b;
  ASSERT_EQ(a.wt_intr.data(), b.wt_intr.data());
}
