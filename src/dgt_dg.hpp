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
    val += U(cell, eq, mode) * phi(pt, mode);
  }
}

}