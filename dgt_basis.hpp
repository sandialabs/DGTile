#pragma once

#include "p3a_for_each.hpp"
#include "p3a_static_vector.hpp"

#include "dgt_defines.hpp"
#include "dgt_views.hpp"

namespace dgt {

static constexpr int max_q = 7;
static constexpr int max_p = 2;

struct Basis {
  int dim = -1;
  int p = -1;
  int nmodes = -1;
  bool tensor = true;
  View<double*>     wt_intr;          // (intr_pt)
  View<double*>     wt_side;          // (side_pt)
  View<double*>     wt_fine;          // (fine_pt)
  View<double**>    pt_intr;          // (intr_pt, dim)
  View<double****>  pt_side;          // (axis, dir, side_pt, dim)
  View<double***>   pt_child_intr;    // (which_child, intr_pt, dim)
  View<double*****> pt_child_side;    // (axis, dir, which_child, side_pt, dim)
  View<double**>    pt_fine;          // (fine_intr_pt, dim)
  View<double**>    pt_eval;          // (all_pt, dim)
  View<double**>    pt_corner;        // (corner_pt, dim)
  View<double**>    phi_intr;         // (intr_pt, mode)
  View<double***>   dphi_intr;        // (axis, intr_pt, mode)
  View<double****>  phi_side;         // (axis, dir, side_pt, mode)
  View<double***>   phi_child_intr;   // (which_child, intr_pt, mode)
  View<double*****> phi_child_side;   // (axis, dir, which_child, side_pt, mode)
  View<double**>    phi_fine;         // (fine_pt, mode)
  View<double**>    phi_eval;         // (eval_pt, mode)
  View<double**>    phi_corner;       // (corner_pt, mode)
  View<double*>     mass;             // (mode)
  void init(int dim, int p, bool tensor);
};

struct HostBasis {
  int dim = -1;
  int p = -1;
  int nmodes = -1;
  bool tensor = true;
  HView<double*>     wt_intr;         // (intr_pt)
  HView<double*>     wt_side;         // (side_pt)
  HView<double*>     wt_fine;         // (fine_pt)
  HView<double**>    pt_intr;         // (intr_pt, dim)
  HView<double****>  pt_side;         // (axis, dir, side_pt, dim)
  HView<double***>   pt_child_intr;   // (which_child, intr_pt, dim)
  HView<double*****> pt_child_side;   // (axis, dir, which_child, side_pt, dim)
  HView<double**>    pt_fine;         // (fine_intr_pt, dim)
  HView<double**>    pt_eval;         // (all_pt, dim)
  HView<double**>    pt_corner;       // (corner_pt, dim)
  HView<double**>    phi_intr;        // (intr_pt, mode)
  HView<double***>   dphi_intr;       // (axis, intr_pt, mode)
  HView<double****>  phi_side;        // (axis, dir, side_pt, mode)
  HView<double***>   phi_child_intr;  // (which_child, intr_pt, mode)
  HView<double*****> phi_child_side;  // (axis, dir, which_child, side_pt, mode)
  HView<double**>    phi_fine;        // (fine_pt, mode)
  HView<double**>    phi_eval;        // (eval_pt, mode)
  HView<double**>    phi_corner;      // (corner_pt, mode)
  HView<double*>     mass;            // (mode)
  void init(int dim, int p, bool tensor);
};

[[nodiscard]] P3A_ALWAYS_INLINE P3A_HOST_DEVICE constexpr
int ipow(int a, int b) {
  int result = 1;
  for (int mult = b; mult > 0; mult--) {
    (void)mult;
    result *= a;
  }
  return result;
}

[[nodiscard]] P3A_ALWAYS_INLINE P3A_HOST_DEVICE constexpr
int num_tensor_modes(int dim, int p) {
  return ipow(p+1, dim);
}

[[nodiscard]] P3A_ALWAYS_INLINE P3A_HOST_DEVICE constexpr
int num_non_tensor_modes(int dim, int p) {
  if (dim == 1) return p+1;
  if (dim == 2) return (p+1)*(p+2)/2;
  if (dim == 3) return (p+1)*(p+2)*(p+3)/6;
  return -1;
}

[[nodiscard]] P3A_ALWAYS_INLINE P3A_HOST_DEVICE constexpr
int num_modes(int dim, int p, bool tensor) {
  if (tensor) return num_tensor_modes(dim, p);
  else return num_non_tensor_modes(dim, p);
}

[[nodiscard]] P3A_ALWAYS_INLINE P3A_HOST_DEVICE constexpr
int num_pts(int dim, int p) {
  return ipow(p+1, dim);
}

[[nodiscard]] P3A_ALWAYS_INLINE P3A_HOST_DEVICE constexpr
int num_corner_pts(int dim) {
  return ipow(2, dim);
}

[[nodiscard]] P3A_ALWAYS_INLINE P3A_HOST_DEVICE constexpr
int num_eval_pts(int dim, int p) {
  return num_pts(dim, p) + dim*ndirs*num_pts(dim-1, p);
}

[[nodiscard]] P3A_ALWAYS_INLINE P3A_HOST_DEVICE constexpr
int num_child(int dim) {
  return ipow(2, dim);
}

[[nodiscard]] P3A_ALWAYS_INLINE P3A_HOST_DEVICE constexpr
int index(p3a::vector3<int> const& ijk, p3a::vector3<int> const& b) {
  return (ijk.z()*b.y() + ijk.y())*b.x() + ijk.x();
}

[[nodiscard]] P3A_ALWAYS_INLINE P3A_HOST_DEVICE constexpr
double legendre(int p, int deriv, double x) {
  double table[max_p+1][max_p+1] = {
    {1.,              0.,   0.},
    {x,               1.,   0.},
    {0.5*(3.*x*x-1.), 3.*x, 3.}
  };
  return table[p][deriv];
}

[[nodiscard]] P3A_ALWAYS_INLINE P3A_HOST_DEVICE constexpr
p3a::vector3<int> tensor_bounds(int dim, int p) {
  p3a::vector3<int> b(1,1,1);
  if (dim > 0) b.x() = p+1;
  if (dim > 1) b.y() = p+1;
  if (dim > 2) b.z() = p+1;
  return b;
}

using ModeVector = p3a::static_vector<double, num_tensor_modes(DIMS, max_p)>;

[[nodiscard]] P3A_ALWAYS_INLINE P3A_HOST_DEVICE inline
ModeVector modes(
    int dim,
    int p,
    bool tensor,
    p3a::vector3<double> const& xi) {
  int m = 0;
  ModeVector phi;
  p3a::vector3<int> const bounds = tensor_bounds(dim, p);
  for (int block = 0; block < p+1; ++block) {
    for (int deg = 0; deg < dim*p + 1; ++deg) {
      for (int k = 0; k < bounds.z(); ++k) {
        for (int j = 0; j < bounds.y(); ++j) {
          for (int i = 0; i < bounds.x(); ++i) {
            int const sum = i+j+k;
            int const idx = p3a::max(i, p3a::max(j, k));
            if ((!tensor) && (sum > p)) continue;
            if ((idx == block) && (sum == deg)) {
              phi[m]               = legendre(i, 0, xi.x());
              if (dim > 1) phi[m] *= legendre(j, 0, xi.y());
              if (dim > 2) phi[m] *= legendre(k, 0, xi.z());
              m++;
            }
          }
        }
      }
    }
  }
  return phi;
}

[[nodiscard]] P3A_ALWAYS_INLINE P3A_HOST_DEVICE inline
ModeVector dmodes(
    int dim,
    int p,
    bool tensor,
    p3a::vector3<int> const& d,
    p3a::vector3<double> const& xi) {
  int m = 0;
  ModeVector dphi;
  p3a::vector3<int> const bounds = tensor_bounds(dim, p);
  for (int block = 0; block < p+1; ++block) {
    for (int deg = 0; deg < dim*p + 1; ++deg) {
      for (int k = 0; k < bounds.z(); ++k) {
        for (int j = 0; j < bounds.y(); ++j) {
          for (int i = 0; i < bounds.x(); ++i) {
            int const sum = i+j+k;
            int const idx = p3a::max(i, p3a::max(j, k));
            if ((!tensor) && (sum > p)) continue;
            if ((idx == block) && (sum == deg)) {
              dphi[m]               = legendre(i, d.x(), xi.x());
              if (dim > 1) dphi[m] *= legendre(j, d.y(), xi.y());
              if (dim > 2) dphi[m] *= legendre(k, d.z(), xi.z());
              m++;
            }
          }
        }
      }
    }
  }
  return dphi;
}

template <class BasisT>
[[nodiscard]] P3A_ALWAYS_INLINE P3A_HOST_DEVICE inline
p3a::vector3<double> get_intr_pt(BasisT const& b, int pt) {
  p3a::vector3<double> xi(0,0,0);
  for (int d = 0; d < b.dim; ++d) {
    xi[d] = b.pt_intr(pt, d);
  }
  return xi;
}

template <class BasisT>
[[nodiscard]] P3A_ALWAYS_INLINE P3A_HOST_DEVICE inline
p3a::vector3<double> get_fine_pt(BasisT const& b, int pt) {
  p3a::vector3<double> xi(0,0,0);
  for (int d = 0; d < b.dim; ++d) {
    xi[d] = b.pt_fine(pt, d);
  }
  return xi;
}

template <class BasisT>
[[nodiscard]] P3A_ALWAYS_INLINE P3A_HOST_DEVICE inline
p3a::vector3<double> get_side_pt(BasisT const& b, int axis, int dir, int pt) {
  p3a::vector3<double> xi(0,0,0);
  for (int d = 0; d < b.dim; ++d) {
    xi[d] = b.pt_side(axis, dir, pt, d);
  }
  return xi;
}

void print_modal_ordering(int dim, int p, bool tensor);

}
