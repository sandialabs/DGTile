#include <dgt_dg.hpp>

#include <gtest/gtest.h>

using namespace dgt;

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
