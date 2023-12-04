#include "dg.hpp"

TEST(dg, maxes)
{
  EXPECT_EQ(max_polynomial_order, 3);
  EXPECT_EQ(max_1D_quadrature_points, 5);
}

TEST(dg, basis_locations)
{
  EXPECT_EQ(basis_locations::CELL, 0);
  EXPECT_EQ(basis_locations::VERTICES, 1);
  EXPECT_EQ(basis_locations::XMIN_FACE, 2);
  EXPECT_EQ(basis_locations::XMAX_FACE, 3);
  EXPECT_EQ(basis_locations::YMIN_FACE, 4);
  EXPECT_EQ(basis_locations::YMAX_FACE, 5);
  EXPECT_EQ(basis_locations::ZMIN_FACE, 6);
  EXPECT_EQ(basis_locations::ZMAX_FACE, 7);
  EXPECT_EQ(basis_locations::EVALUATION, 8);
  EXPECT_EQ(basis_locations::NUM, 9);
}

TEST(dg, face_basis_locations)
{
  EXPECT_EQ(basis_locations::face(X, LEFT), basis_locations::XMIN_FACE);
  EXPECT_EQ(basis_locations::face(X, RIGHT), basis_locations::XMAX_FACE);
  EXPECT_EQ(basis_locations::face(Y, LEFT), basis_locations::YMIN_FACE);
  EXPECT_EQ(basis_locations::face(Y, RIGHT), basis_locations::YMAX_FACE);
  EXPECT_EQ(basis_locations::face(Z, LEFT), basis_locations::ZMIN_FACE);
  EXPECT_EQ(basis_locations::face(Z, RIGHT), basis_locations::ZMAX_FACE);
}

TEST(dg, ipow)
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

TEST(dg, num_tensor_modes)
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

TEST(dg, num_non_tensor_modes)
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

TEST(dg, num_modes)
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

TEST(dg, num_gauss_points)
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

TEST(dg, num_evaluation_points)
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

TEST(dg, num_vertices)
{
  EXPECT_EQ(num_vertices(1), 2);
  EXPECT_EQ(num_vertices(2), 4);
  EXPECT_EQ(num_vertices(3), 8);
}

TEST(dg, num_edges)
{
  EXPECT_EQ(num_edges(1), 2);
  EXPECT_EQ(num_edges(2), 4);
  EXPECT_EQ(num_edges(3), 12);
}

TEST(dg, num_faces)
{
  EXPECT_EQ(num_faces(1), 2);
  EXPECT_EQ(num_faces(2), 4);
  EXPECT_EQ(num_faces(3), 6);
}

TEST(dg, get_guass_weight_sums)
{
  for (int q = 1; q <= max_1D_quadrature_points; ++q) {
    real sum = 0.;
    for (int pt = 0; pt < max_1D_quadrature_points; ++pt) {
      sum += get_gauss_weight(q, pt);
    }
    EXPECT_NEAR(sum, 2., 1.e-15);
  }
}

TEST(dg, get_gauss_point_sums)
{
  for (int q = 1; q <= max_1D_quadrature_points; ++q) {
    real sum = 0.;
    for (int pt = 0; pt < max_1D_quadrature_points; ++pt) {
      sum += get_gauss_point(q, pt);
    }
    EXPECT_NEAR(sum, 0., 1.e-15);
  }
}

TEST(dg, tensor_bounds)
{
  EXPECT_EQ(tensor_bounds(1, 0), Vec3<int>(1,0,0));
  EXPECT_EQ(tensor_bounds(1, 1), Vec3<int>(2,0,0));
  EXPECT_EQ(tensor_bounds(1, 2), Vec3<int>(3,0,0));
  EXPECT_EQ(tensor_bounds(1, 3), Vec3<int>(4,0,0));
  EXPECT_EQ(tensor_bounds(1, 4), Vec3<int>(5,0,0));
  EXPECT_EQ(tensor_bounds(2, 0), Vec3<int>(1,1,0));
  EXPECT_EQ(tensor_bounds(2, 1), Vec3<int>(2,2,0));
  EXPECT_EQ(tensor_bounds(2, 2), Vec3<int>(3,3,0));
  EXPECT_EQ(tensor_bounds(2, 3), Vec3<int>(4,4,0));
  EXPECT_EQ(tensor_bounds(2, 4), Vec3<int>(5,5,0));
  EXPECT_EQ(tensor_bounds(3, 0), Vec3<int>(1,1,1));
  EXPECT_EQ(tensor_bounds(3, 1), Vec3<int>(2,2,2));
  EXPECT_EQ(tensor_bounds(3, 2), Vec3<int>(3,3,3));
  EXPECT_EQ(tensor_bounds(3, 3), Vec3<int>(4,4,4));
  EXPECT_EQ(tensor_bounds(3, 4), Vec3<int>(5,5,5));
}

static real function(int p, Vec3<real> const& x)
{
  real val = 1.;
  if (p > 0) {
    val +=
      std::pow(x.x(), p) +
      std::pow(x.y(), p) +
      std::pow(x.z(), p);
  }
  return val;
}

static void test_integrates_polynomial(Basis<HostView> const& B)
{
  real integral = 0.;
  for (int pt = 0; pt < B.num_cell_pts; ++pt) {
    Vec3<real> const xi = get_point(B, basis_locations::CELL, pt);
    real const f = function(B.p, xi);
    real const wt = B.cell_weights(pt);
    integral += f * wt;
  }
  if (B.p != 2) {
    EXPECT_NEAR(integral, ipow(2, B.dim), 1.e-14);
  }
}

void test_basis(
    int const dim,
    int const p,
    int const q,
    bool const tensor)
{
  auto B = build_basis<HostView>(dim, p, q, tensor);
  test_integers(B, dim, p, q, tensor);
  test_view_sizes(B);
  test_integrates_polynomial(B);
}
