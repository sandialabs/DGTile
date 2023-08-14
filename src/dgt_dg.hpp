#pragma once

#include "dgt_grid3.hpp"

namespace dgt {

static constexpr int max_polynomial_order = 3;
static constexpr int max_1D_quadrature_points = 5;

namespace basis_locations {
static constexpr int CELL = 0;
static constexpr int VERTICES = 1;
static constexpr int XMIN_FACE = 2;
static constexpr int XMAX_FACE = 3;
static constexpr int YMIN_FACE = 4;
static constexpr int YMAX_FACE = 5;
static constexpr int ZMIN_FACE = 6;
static constexpr int ZMAX_FACE = 7;
static constexpr int EVALUATION = 8;
static constexpr int NUM = 9;
}

template <template <class> class ViewT>
struct TabulatedBasis
{
  ViewT<real**> points;
  ViewT<real**> phis;
  ViewT<real***> grad_phis;
};

template <template <class> class ViewT>
struct Basis
{
  public:
    int dim = -1;
    int p = -1;
    int q = -1;
    bool tensor = true;
  public:
    int num_modes = -1;
    int num_cell_pts = -1;
    int num_face_pts = -1;
  public:
    ViewT<real*> cell_weights;
    ViewT<real*> face_weights;
    TabulatedBasis<ViewT> modes[basis_locations::NUM];
};

template <template <class> class ViewT>
Basis<ViewT> build_basis(
    int const dim,
    int const p,
    int const q,
    bool const tensor);

DGT_METHOD constexpr int ipow(int const a, int const b)
{
  int result = 1;
  for (int mult = b; mult > 0; mult--) {
    (void)mult;
    result *= a;
  }
  return result;
}

DGT_METHOD constexpr int num_gauss_points(int const dim, int const q)
{
  return ipow(q, dim);
}

DGT_METHOD constexpr int num_evaluation_points(int const dim, int const q)
{
  return
    num_gauss_points(dim, q) +
    num_gauss_points(dim-1, q) * dim * DIRECTIONS;
}

DGT_METHOD constexpr int num_vertices(int const dim)
{
  return ipow(2, dim);
}

DGT_METHOD constexpr int num_tensor_modes(int const dim, int const p)
{
  return ipow(p+1, dim);
}

DGT_METHOD constexpr int num_non_tensor_modes(int const dim, int const p)
{
  if (dim == 1) return p+1;
  if (dim == 2) return (p+1)*(p+2)/2;
  if (dim == 3) return (p+1)*(p+2)*(p+3)/6;
  return -1;
}

DGT_METHOD constexpr int num_modes(int const dim, int const p, bool const tensor)
{
  if (tensor) return num_tensor_modes(dim, p);
  else return num_non_tensor_modes(dim, p);
}

DGT_METHOD constexpr real get_gauss_weight(int const q, int const pt)
{
  real table[max_1D_quadrature_points][max_1D_quadrature_points] = {
    {2.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000},
    {1.0000000000000000, 1.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000},
    {0.5555555555555556, 0.8888888888888888, 0.5555555555555556, 0.0000000000000000, 0.0000000000000000},
    {0.3478548451374538, 0.6521451548625461, 0.6521451548625461, 0.3478548451374538, 0.0000000000000000},
    {0.2369268850561891, 0.4786286704993665, 0.5688888888888889, 0.4786286704993665, 0.2369268850561891}
  };
  return table[q-1][pt];
}

DGT_METHOD constexpr real get_gauss_point(int const q, int const pt)
{
  real table[max_1D_quadrature_points][max_1D_quadrature_points] = {
    { 0.0000000000000000,  0.0000000000000000,  0.0000000000000000,  0.0000000000000000,  0.0000000000000000},
    {-0.5773502691896257,  0.5773502691896257,  0.0000000000000000,  0.0000000000000000,  0.0000000000000000},
    {-0.7745966692414834,  0.0000000000000000,  0.7745966692414834,  0.0000000000000000,  0.0000000000000000},
    {-0.8611363115940526, -0.3399810435848563,  0.3399810435848563,  0.8611363115940526,  0.0000000000000000},
    {-0.9061798459386640, -0.5384693101056831,  0.0000000000000000,  0.5384693101056831,  0.9061798459386640}
  };
  return table[q-1][pt];
}

DGT_METHOD constexpr real get_legendre(int const p, int const deriv, real const x)
{
  real table[max_polynomial_order+1][max_polynomial_order+1] = {
    {1.,                  0.,               0.,    0.},
    {x,                   1.,               0.,    0.},
    {0.5*(3.*x*x-1.),     3.*x,             3.,    0.},
    {0.5*(5.*x*x*x-3.*x), 1.5*(5.*x*x-1.),  15.*x, 15.}
  };
  return table[p][deriv];
}

DGT_METHOD constexpr Grid3 tensor_bounds(int const dim, int const p)
{
  Vec3<int> b = Vec3<int>::zero();
  if (dim > 0) b.x() = p+1;
  if (dim > 1) b.y() = p+1;
  if (dim > 2) b.z() = p+1;
  return Grid3(b);
}

template <class ModalT, class BasisT>
DGT_METHOD inline real eval(
    ModalT const U,
    int const cell,
    int const eq,
    BasisT const B,
    int const location,
    int const pt) {
  real val = U(cell, eq, 0) * B.modes[location](pt, 0);
  for (int mode = 1; mode < B.num_modes; ++mode) {
    val += U(cell, eq, mode) * B.modes[location](pt, mode);
  }
  return val;
}

}
