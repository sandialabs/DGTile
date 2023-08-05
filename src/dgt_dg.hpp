#pragma once

#include "dgt_defines.hpp"
#include "dgt_macros.hpp"

namespace dgt {

namespace basis {
static constexpr int INTERIOR = 0;
static constexpr int VERTICES = 1;
static constexpr int XMIN_FACE = 2;
static constexpr int XMAX_FACE = 3;
static constexpr int YMIN_FACE = 4;
static constexpr int YMAX_FACE = 5;
static constexpr int ZMIN_FACE = 6;
static constexpr int ZMAX_FACE = 7;
static constexpr int EVALUATION = 8;
static constexpr int LOCATIONS = 9;
}

template <template <class> class ViewT>
struct GaussianQuadrature
{
  ViewT<real*> wt;
  ViewT<real**> pts;
};

template <template <class> class ViewT>
struct TabulatedBasis
{
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
    int num_intr_pts = -1;
    int num_face_pts = -1;
  public:
    GaussianQuadrature<ViewT> intr_quadrature;
    GaussianQuadrature<ViewT> face_quadrature;
    TabulatedBasis<ViewT> modes[basis::LOCATIONS];
};

template <template <class> class ViewT>
Basis<ViewT> build_basis();

DGT_METHOD constexpr int ipow(int a, int b)
{
  int result = 1;
  for (int mult = b; mult > 0; mult--) {
    (void)mult;
    result *= a;
  }
  return result;
}

DGT_METHOD real get_legendre(
    int const n,
    double const x)
{
  if (n == 0) return 1;
  if (n == 1) return x;
  real const term1 = (2.*real(n)-1.) * x * get_legendre(n-1, x);
  real const term2 = (real(n)-1.) * get_legendre(n-2, x);
  return (term1 - term2) / real(n);
}

DGT_METHOD real get_dlegendre(
    int const n,
    double const x)
{
  if (n == 0) return 0.;
  if (n == 1) return 1.;
  real const term1 = (2.*real(n)-1.) * get_legendre(n-1, x);
  real const term2 = (2.*real(n)-1.) * x * get_dlegendre(n-1, x);
  real const term3 = (real(n)-1.) * get_dlegendre(n-2, x);
  return (term1 + term2 - term3) / real(n);
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
}

}
