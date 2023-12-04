#include <dgt_interpolate.hpp>

#include <gtest/gtest.h>

using namespace dgt;

struct Table
{
  std::vector<double> x;
  std::vector<double> y;
  std::vector<double> f;
};

static double f_example(real const x, real const y)
{
  return 1. + 2.*x + 3.*y + 4.*x*y;
}

static Table setup_example_table()
{
  Table t;
  int const nx = 4;
  int const ny = 4;
  t.x.resize(nx);
  t.y.resize(ny);
  t.f.resize(nx*ny);
  t.x[0] = 0.;
  t.x[1] = 1.;
  t.x[2] = 2.;
  t.x[3] = 3.;
  t.y[0] = 0.;
  t.y[1] = 1.;
  t.y[2] = 2.;
  t.y[3] = 3.;
  for (int i = 0; i < nx; ++i) {
    for (int j = 0; j < ny; ++j) {
      real const x = t.x[i];
      real const y = t.y[j];
      t.f[j*nx + i] = f_example(x, y);
    }
  }
  return t;
}

TEST(interpolate, binary_search)
{
  std::vector<double> a = {1., 2.12, 4.87, 7.5, 7.6, 278.1, 401.2};
  EXPECT_EQ(binary_search(a, 0.), 0);
  EXPECT_EQ(binary_search(a, 1.5), 0);
  EXPECT_EQ(binary_search(a, 3.2111), 1);
  EXPECT_EQ(binary_search(a, 7.1), 2);
  EXPECT_EQ(binary_search(a, 7.5), 3);
  EXPECT_EQ(binary_search(a, 7.555), 3);
  EXPECT_EQ(binary_search(a, 101.555), 4);
  EXPECT_EQ(binary_search(a, 281.1), 5);
  EXPECT_EQ(binary_search(a, 450.), 5);
}


TEST(interpolate, bilinear_2D_base)
{
  Table const t = setup_example_table();
  EXPECT_NEAR(
      bilinear::interpolate(t.f, t.x, t.y, 0.5, 0.5, 0, 0),
      f_example(0.5, 0.5),
      1.e-14);
  EXPECT_NEAR(
      bilinear::interpolate(t.f, t.x, t.y, 1.5, 1.5, 1, 1),
      f_example(1.5, 1.5),
      1.e-14);
  EXPECT_NEAR(
      bilinear::interpolate(t.f, t.x, t.y, 2.5, 2.5, 2, 2),
      f_example(2.5, 2.5),
      1.e-14);
  EXPECT_NEAR(
      bilinear::interpolate(t.f, t.x, t.y, 0.75, 1.75, 0, 1),
      f_example(0.75, 1.75),
      1.e-14);
}
