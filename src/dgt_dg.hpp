#pragma once

#include "dgt_defines.hpp"
#include "dgt_macros.hpp"

namespace dgt {

enum {
  INTERIOR = 0,
  XMIN_FACE = 1,
  XMAX_FACE = 2,
  YMIN_FACE = 3,
  YMAX_FACE = 4,
  ZMIN_FACE = 5,
  ZMAX_FACE = 6,
  VERTICES = 7,
  BASIS_LOCATIONS = 8
};

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
    TabulatedBasis<ViewT> modes[BASIS_LOCATIONS];
};

template <template <class> class ViewT>
Basis<ViewT> build_basis();

DGT_METHOD constexpr real get_legendre(
    int const n,
    double const x)
{
  if (n == 0) return 1;
  if (n == 1) return x;
  real const term1 = (2.*real(n)-1.) * x * get_legendre(n-1, x);
  real const term2 = (real(n)-1.) * get_legendre(n-2, x);
  return (term1 - term2) / real(n);
}

DGT_METHOD constexpr real get_dlegendre(
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
