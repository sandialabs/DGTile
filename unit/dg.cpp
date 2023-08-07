#include <dgt_dg.hpp>

#include <gtest/gtest.h>

using namespace dgt;

TEST(basis, tabulated_locations)
{
  EXPECT_EQ(basis::INTERIOR, 0);
  EXPECT_EQ(basis::VERTICES, 1);
  EXPECT_EQ(basis::XMIN_FACE, 2);
  EXPECT_EQ(basis::XMAX_FACE, 3);
  EXPECT_EQ(basis::YMIN_FACE, 4);
  EXPECT_EQ(basis::YMAX_FACE, 5);
  EXPECT_EQ(basis::ZMIN_FACE, 6);
  EXPECT_EQ(basis::ZMAX_FACE, 7);
  EXPECT_EQ(basis::EVALUATION, 8);
  EXPECT_EQ(basis::LOCATIONS, 9);
}

TEST(basis, ipow)
{
  EXPECT_EQ(ipow(1, 1), 1);
  EXPECT_EQ(ipow(1, 2), 1);
  EXPECT_EQ(ipow(1, 3), 1);
  EXPECT_EQ(ipow(2, 1), 2);
  EXPECT_EQ(ipow(2, 2), 4);
  EXPECT_EQ(ipow(2, 3), 8);
  EXPECT_EQ(ipow(3, 1), 3);
  EXPECT_EQ(ipow(3, 2), 9);
  EXPECT_EQ(ipow(3, 3), 27);
  EXPECT_EQ(ipow(4, 1), 4);
  EXPECT_EQ(ipow(4, 2), 16);
  EXPECT_EQ(ipow(4, 3), 64);
}

TEST(basis, num_tensor_modes)
{
  EXPECT_EQ(num_tensor_modes(1, 0), 1);
  EXPECT_EQ(num_tensor_modes(1, 1), 2);
  EXPECT_EQ(num_tensor_modes(1, 2), 3);
  EXPECT_EQ(num_tensor_modes(1, 3), 4);
  EXPECT_EQ(num_tensor_modes(2, 0), 1);
  EXPECT_EQ(num_tensor_modes(2, 1), 4);
  EXPECT_EQ(num_tensor_modes(2, 2), 9);
  EXPECT_EQ(num_tensor_modes(2, 3), 16);
  EXPECT_EQ(num_tensor_modes(3, 0), 1);
  EXPECT_EQ(num_tensor_modes(3, 1), 8);
  EXPECT_EQ(num_tensor_modes(3, 2), 27);
  EXPECT_EQ(num_tensor_modes(3, 3), 64);
}

TEST(basis, num_non_tensor_modes)
{
  EXPECT_EQ(num_non_tensor_modes(1, 0), 1);
  EXPECT_EQ(num_non_tensor_modes(1, 1), 2);
  EXPECT_EQ(num_non_tensor_modes(1, 2), 3);
  EXPECT_EQ(num_non_tensor_modes(1, 3), 4);
  EXPECT_EQ(num_non_tensor_modes(1, 4), 5);
  EXPECT_EQ(num_non_tensor_modes(1, 5), 6);
  EXPECT_EQ(num_non_tensor_modes(2, 0), 1);
  EXPECT_EQ(num_non_tensor_modes(2, 1), 3);
  EXPECT_EQ(num_non_tensor_modes(2, 2), 6);
  EXPECT_EQ(num_non_tensor_modes(2, 3), 10);
  EXPECT_EQ(num_non_tensor_modes(2, 4), 15);
  EXPECT_EQ(num_non_tensor_modes(2, 5), 21);
  EXPECT_EQ(num_non_tensor_modes(3, 0), 1);
  EXPECT_EQ(num_non_tensor_modes(3, 1), 4);
  EXPECT_EQ(num_non_tensor_modes(3, 2), 10);
  EXPECT_EQ(num_non_tensor_modes(3, 3), 20);
  EXPECT_EQ(num_non_tensor_modes(3, 4), 35);
  EXPECT_EQ(num_non_tensor_modes(3, 5), 56);
}

TEST(basis, num_modes)
{
  EXPECT_EQ(num_modes(2, 0, true), num_tensor_modes(2, 0));
  EXPECT_EQ(num_modes(2, 1, true), num_tensor_modes(2, 1));
  EXPECT_EQ(num_modes(2, 2, true), num_tensor_modes(2, 2));
  EXPECT_EQ(num_modes(3, 0, true), num_tensor_modes(3, 0));
  EXPECT_EQ(num_modes(3, 1, true), num_tensor_modes(3, 1));
  EXPECT_EQ(num_modes(3, 2, true), num_tensor_modes(3, 2));
  EXPECT_EQ(num_modes(2, 0, false), num_non_tensor_modes(2, 0));
  EXPECT_EQ(num_modes(2, 1, false), num_non_tensor_modes(2, 1));
  EXPECT_EQ(num_modes(2, 2, false), num_non_tensor_modes(2, 2));
  EXPECT_EQ(num_modes(3, 0, false), num_non_tensor_modes(3, 0));
  EXPECT_EQ(num_modes(3, 1, false), num_non_tensor_modes(3, 1));
  EXPECT_EQ(num_modes(3, 2, false), num_non_tensor_modes(3, 2));
}

TEST(basis, num_gauss_points)
{
  EXPECT_EQ(num_gauss_points(1, 1), 1);
  EXPECT_EQ(num_gauss_points(1, 2), 2);
  EXPECT_EQ(num_gauss_points(1, 3), 3);
  EXPECT_EQ(num_gauss_points(1, 4), 4);
  EXPECT_EQ(num_gauss_points(2, 1), 1);
  EXPECT_EQ(num_gauss_points(2, 2), 4);
  EXPECT_EQ(num_gauss_points(2, 3), 9);
  EXPECT_EQ(num_gauss_points(2, 4), 16);
  EXPECT_EQ(num_gauss_points(3, 1), 1);
  EXPECT_EQ(num_gauss_points(3, 2), 8);
  EXPECT_EQ(num_gauss_points(3, 3), 27);
  EXPECT_EQ(num_gauss_points(3, 4), 64);
}

TEST(basis, num_evaluation_points)
{
  EXPECT_EQ(num_evaluation_points(1, 1), 3);
  EXPECT_EQ(num_evaluation_points(1, 2), 4);
  EXPECT_EQ(num_evaluation_points(1, 3), 5);
  EXPECT_EQ(num_evaluation_points(1, 4), 6);
  EXPECT_EQ(num_evaluation_points(2, 1), 5);
  EXPECT_EQ(num_evaluation_points(2, 2), 12);
  EXPECT_EQ(num_evaluation_points(2, 3), 21);
  EXPECT_EQ(num_evaluation_points(2, 4), 32);
  EXPECT_EQ(num_evaluation_points(3, 1), 7);
  EXPECT_EQ(num_evaluation_points(3, 2), 32);
  EXPECT_EQ(num_evaluation_points(3, 3), 81);
  EXPECT_EQ(num_evaluation_points(3, 4), 160);
}

TEST(basis, num_vertices)
{
  EXPECT_EQ(num_vertices(1), 2);
  EXPECT_EQ(num_vertices(2), 4);
  EXPECT_EQ(num_vertices(3), 8);
}

static real Pn(int const n, real const x)
{
  if (n == 0) return 1.;
  if (n == 1) return x;
  if (n == 2) return 0.5*(3.*x*x-1.);
  if (n == 3) return 0.5*(5.*x*x*x-3.*x);
  if (n == 4) return 0.125*(35.*x*x*x*x - 30.*x*x + 3.);
  if (n == 5) return 0.125*(63.*x*x*x*x*x -70.*x*x*x + 15.*x);
  return 0.;
}

static real dPn(int const n, real const x)
{
  if (n == 0) return 0.;
  if (n == 1) return 1.;
  if (n == 2) return 3.*x;
  if (n == 3) return 1.5*(5.*x*x - 1.);
  if (n == 4) return 2.5*(7.*x*x*x - 3.*x);
  if (n == 5) return 15./8.*(21.*x*x*x*x - 14.*x*x + 1.);
  return 0.;
}

static void test_legendre(int const n, real const tol)
{
  EXPECT_NEAR(get_legendre(n, -1.),   Pn(n, -1.),   tol);
  EXPECT_NEAR(get_legendre(n, -0.75), Pn(n, -0.75), tol);
  EXPECT_NEAR(get_legendre(n, -0.5),  Pn(n, -0.5),  tol);
  EXPECT_NEAR(get_legendre(n, -0.25), Pn(n, -0.25), tol);
  EXPECT_NEAR(get_legendre(n,  0.),   Pn(n,  0.),   tol);
  EXPECT_NEAR(get_legendre(n,  0.25), Pn(n,  0.25), tol);
  EXPECT_NEAR(get_legendre(n,  0.5),  Pn(n,  0.5),  tol);
  EXPECT_NEAR(get_legendre(n,  0.75), Pn(n,  0.75), tol);
  EXPECT_NEAR(get_legendre(n,  1.),   Pn(n,  1.),   tol);
}

static void test_dlegendre(int const n, real const tol)
{
  EXPECT_NEAR(get_dlegendre(n, -1.),   dPn(n, -1.),   tol);
  EXPECT_NEAR(get_dlegendre(n, -0.75), dPn(n, -0.75), tol);
  EXPECT_NEAR(get_dlegendre(n, -0.5),  dPn(n, -0.5),  tol);
  EXPECT_NEAR(get_dlegendre(n, -0.25), dPn(n, -0.25), tol);
  EXPECT_NEAR(get_dlegendre(n,  0.),   dPn(n,  0.),   tol);
  EXPECT_NEAR(get_dlegendre(n,  0.25), dPn(n,  0.25), tol);
  EXPECT_NEAR(get_dlegendre(n,  0.5),  dPn(n,  0.5),  tol);
  EXPECT_NEAR(get_dlegendre(n,  0.75), dPn(n,  0.75), tol);
  EXPECT_NEAR(get_dlegendre(n,  1.),   dPn(n,  1.),   tol);
}

TEST(basis, get_legendre_n0)
{
  test_legendre(0, 0.);
}

TEST(basis, get_legendre_n1)
{
  test_legendre(1, 0.);
}

TEST(basis, get_legendre_n2)
{
  test_legendre(2, 1.e-15);
}

TEST(basis, get_legendre_n3)
{
  test_legendre(3, 1.e-15);
}

TEST(basis, get_legendre_n4)
{
  test_legendre(4, 1.e-15);
}

TEST(basis, get_legendre_n5)
{
  test_legendre(5, 1.e-15);
}

TEST(basis, get_dlegendre_n0)
{
  test_dlegendre(0, 0.);
}

TEST(basis, get_dlegendre_n1)
{
  test_dlegendre(1, 0.);
}

TEST(basis, get_dlegendre_n2)
{
  test_dlegendre(2, 1.e-15);
}

TEST(basis, get_dlegendre_n3)
{
  test_dlegendre(3, 1.e-15);
}

TEST(basis, get_dlegendre_n4)
{
  test_dlegendre(4, 1.e-15);
}

TEST(basis, get_dlegendre_n5)
{
  test_dlegendre(5, 1.e-15);
}
