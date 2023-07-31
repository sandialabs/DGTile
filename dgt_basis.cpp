#include <cmath>
#include <limits>

#include "caliper/cali.h"

#include "dgt_amr.hpp"
#include "dgt_basis.hpp"
#include "dgt_spatial.hpp"

namespace dgt {

static double nan = std::numeric_limits<double>::quiet_NaN();

static void verify_dim(int dim) {
  if ((dim <= 0) || (dim > DIMS)) {
    throw std::runtime_error("Basis - invalid dim");
  }
}

static void verify_p(int p) {
  if ((p < 0) || (p > max_p)) {
    throw std::runtime_error("Basis - invalid p");
  }
}

static double gauss_wt(int q, int pt) {
  static double table[max_q][max_q] = {
    {2.,                 nan,                nan,                nan,                nan,                nan,                nan},
    {1.0000000000000000, 1.0000000000000000, nan,                nan,                nan,                nan,                nan},
    {0.5555555555555556, 0.8888888888888888, 0.5555555555555556, nan,                nan,                nan,                nan},
    {0.3478548451374538, 0.6521451548625461, 0.6521451548625461, 0.3478548451374538, nan,                nan,                nan},
    {0.2369268850561891, 0.4786286704993665, 0.5688888888888889, 0.4786286704993665, 0.2369268850561891, nan,                nan},
    {0.1713244923791704, 0.3607615730481386, 0.4679139345726910, 0.4679139345726910, 0.3607615730481386, 0.1713244923791704, nan},
    {0.1294849661688697, 0.2797053914892766, 0.3818300505051189, 0.4179591836734694, 0.3818300505051189, 0.2797053914892766, 0.1294849661688697}
  };
  return table[q][pt];
}

static double gauss_pt(int q, int pt) {
  static double table[max_q][max_q] = {
    { 0.0000000000000000,  nan,                 nan,                 nan,                 nan,                 nan,                 nan},
    {-0.5773502691896257,  0.5773502691896257,  nan,                 nan,                 nan,                 nan,                 nan},
    {-0.7745966692414834,  0.0000000000000000,  0.7745966692414834,  nan,                 nan,                 nan,                 nan},
    {-0.8611363115940526, -0.3399810435848563,  0.3399810435848563,  0.8611363115940526,  nan,                 nan,                 nan},
    {-0.9061798459386640, -0.5384693101056831,  0.0000000000000000,  0.5384693101056831,  0.9061798459386640,  nan,                 nan},
    {-0.9324695142031521, -0.6612093864662645, -0.2386191860831969,  0.2386191860831969,  0.6612093864662645,  0.9324695142031521,  nan},
    {-0.9491079123427585, -0.7415311855993945, -0.4058451513773972,  0.0000000000000000,  0.4058451513773972,  0.7415311855993945,  0.9491079123427585}
  };
  return table[q][pt];
}

static double corner_pt(int pt) {
  static constexpr int nc = 2;
  static double table[nc] = {-1.0, 1.0};
  return table[pt];
}

static double shift(double xi, double x0, double x1) {
  return 0.5*x0*(1.-xi) + 0.5*x1*(1.+xi);
}

static p3a::vector3<double> get_child_xi(
    p3a::vector3<double> const& xi, int which_child) {
  p3a::vector3<double> child_xi;
  p3a::vector3<int> const local = get_local(which_child);
  double const x[3] = {-1., 0., 1.};
  for (int axis = 0; axis < DIMS; ++axis) {
    int const idx = local[axis];
    child_xi[axis] = shift(xi[axis], x[idx], x[idx + 1]);
  }
  return child_xi;
}

static HView<double*> get_wt_intr(int dim, int p) {
  HView<double*> wts("", num_pts(dim, p));
  p3a::vector3<int> const bounds = tensor_bounds(dim, p);
  auto f = [&] (p3a::vector3<int> const& ijk) {
    int const pt = index(ijk, bounds);
    wts(pt)               = gauss_wt(p, ijk.x());
    if (dim > 1) wts(pt) *= gauss_wt(p, ijk.y());
    if (dim > 2) wts(pt) *= gauss_wt(p, ijk.z());
  };
  p3a::for_each(p3a::execution::seq, bounds, f);
  return wts;
}

static HView<double*> get_wt_side(int dim, int p) {
  HView<double*> wts("", num_pts(dim-1, p));
  p3a::vector3<int> const bounds = tensor_bounds(dim-1, p);
  auto f = [&] (p3a::vector3<int> const& ijk) {
    int const pt = index(ijk, bounds);
    wts(pt)               = 1.;
    if (dim > 1) wts(pt) *= gauss_wt(p, ijk.x());
    if (dim > 2) wts(pt) *= gauss_wt(p, ijk.y());
  };
  p3a::for_each(p3a::execution::seq, bounds, f);
  return wts;
}

static HView<double*> get_wt_fine(int dim, int p) {
  (void)p;
  return get_wt_intr(dim, 3);
}

static HView<double**> get_pt_intr(int dim, int p) {
  HView<double**> pts("", num_pts(dim, p), dim);
  p3a::vector3<int> const bounds = tensor_bounds(dim, p);
  auto f = [&] (p3a::vector3<int> const& ijk) {
    int const pt = index(ijk, bounds);
    for (int axis = 0; axis < dim; ++axis) {
      pts(pt, axis) = gauss_pt(p, ijk[axis]);
    }
  };
  p3a::for_each(p3a::execution::seq, bounds, f);
  return pts;
}

static HView<double**> get_pt_corner(int dim) {
  int constexpr p = 1;
  HView<double**> pts("", num_corner_pts(dim), dim);
  p3a::vector3<int> const bounds = tensor_bounds(dim, p);
  auto f = [&] (p3a::vector3<int> const& ijk) {
    int const pt = index(ijk, bounds);
    for (int axis = 0; axis < dim; ++axis) {
      pts(pt, axis) = corner_pt(ijk[axis]);
    }
  };
  p3a::for_each(p3a::execution::seq, bounds, f);
  return pts;
}

static HView<double****> get_pt_side(int dim, int p) {
  HView<double****> pts("", dim, ndirs, num_pts(dim-1, p), dim);
  p3a::vector3<int> const bounds = tensor_bounds(dim-1, p);
  for (int axis = 0; axis < dim; ++axis) {
    for (int dir = 0; dir < ndirs; ++dir) {
      auto f = [&] (p3a::vector3<int> const& ijk) {
        int const pt = index(ijk, bounds);
        int const ii = (axis == X) ? Y : X;
        int const jj = (axis == Z) ? Y : Z;
        pts(axis, dir, pt, axis) = get_dir_sign(dir);
        if (dim > 1) pts(axis, dir, pt, ii) = gauss_pt(p, ijk.x());
        if (dim > 2) pts(axis, dir, pt, jj) = gauss_pt(p, ijk.y());
      };
      p3a::for_each(p3a::execution::seq, bounds, f);
    }
  }
  return pts;
}

static HView<double***> get_pt_child_intr(int dim, int p) {
  HView<double***> pts("", num_child(dim), num_pts(dim, p), dim);
  p3a::vector3<double> xi(0,0,0);
  p3a::vector3<int> const bounds = tensor_bounds(dim, p);
  HView<double**> const intr_pts = get_pt_intr(dim, p);
  for (int which_child = 0; which_child < num_child(dim); ++which_child) {
    auto f = [&] (p3a::vector3<int> const& ijk) {
      int const pt = index(ijk, bounds);
      xi.x()              = intr_pts(pt, X);
      if (dim > 1) xi.y() = intr_pts(pt, Y);
      if (dim > 2) xi.z() = intr_pts(pt, Z);
      p3a::vector3<double> const child_xi = get_child_xi(xi, which_child);
      for (int axis = 0; axis < dim; ++axis) {
        pts(which_child, pt, axis) = child_xi[axis];
      }
    };
    p3a::for_each(p3a::execution::seq, bounds, f);
  }
  return pts;
}

static HView<double*****> get_pt_child_side(int dim, int p) {
  HView<double*****> pts("", dim, ndirs, num_child(dim-1), num_pts(dim-1, p), dim);
  p3a::vector3<double> xi(0,0,0);
  p3a::vector3<int> const bounds = tensor_bounds(dim-1, p);
  HView<double****> const side_pts = get_pt_side(dim, p);
  for (int axis = 0; axis < dim; ++axis) {
    for (int dir = 0; dir < ndirs; ++dir) {
      for (int which_child = 0; which_child < num_child(dim-1); ++which_child) {
        auto f = [&] (p3a::vector3<int> const& ijk) {
          int const pt = index(ijk, bounds);
          int const ii = (axis == X) ? Y : X;
          int const jj = (axis == Z) ? Y : Z;
          if (dim > 1) xi.x() = side_pts(axis, dir, pt, ii);
          if (dim > 2) xi.y() = side_pts(axis, dir, pt, jj);
          p3a::vector3<double> const child_xi = get_child_xi(xi, which_child);
          pts(axis, dir, which_child, pt, axis) = side_pts(axis, dir, pt, axis); 
          if (dim > 1) pts(axis, dir, which_child, pt, ii) = child_xi.x();
          if (dim > 2) pts(axis, dir, which_child, pt, jj) = child_xi.y();
        };
        p3a::for_each(p3a::execution::seq, bounds, f);
      }
    }
  }
  return pts;
}

static HView<double**> get_pt_fine(int dim, int p) {
  (void)p;
  return get_pt_intr(dim, 3);
}

static HView<double**> get_pt_eval(int dim, int p) {
  int const neval_pts = num_eval_pts(dim, p);
  HView<double**> const intr_pts = get_pt_intr(dim, p);
  HView<double****> const side_pts = get_pt_side(dim, p);
  HView<double**> pts("", neval_pts, dim);
  int eval_pt = 0;
  for (int pt = 0; pt < num_pts(dim, p); ++pt) {
    for (int d = 0; d < dim; ++d) {
      pts(eval_pt, d) = intr_pts(pt, d);
    }
    eval_pt++;
  }
  for (int axis = 0; axis < dim; ++axis) {
    for (int dir = 0; dir < ndirs; ++dir) {
      for (int pt = 0; pt < num_pts(dim-1, p); ++pt) {
        for (int d = 0; d < dim; ++d) {
          pts(eval_pt, d) = side_pts(axis, dir, pt, d);
        }
        eval_pt++;
      }
    }
  }
  return pts;
}

static HView<double**> get_phi(int dim, int p, bool tensor, HView<double**> pts) {
  int const npts = pts.extent(0);
  int const nmodes = num_modes(dim, p, tensor);
  p3a::vector3<double> xi(0,0,0);
  HView<double**> phi("", npts, nmodes);
  for (int pt = 0; pt < npts; ++pt) {
    for (int d = 0; d < dim; ++d) {
      xi[d] = pts(pt, d);
    }
    auto const phi_xi = modes(dim, p, tensor, xi);
    for (int m = 0; m < nmodes; ++m) {
      phi(pt, m) = phi_xi[m];
    }
  }
  return phi;
}

static HView<double***> get_dphi(int dim, int p, bool tensor, HView<double**> pts) {
  int const npts = pts.extent(0);
  int const nmodes = num_modes(dim, p, tensor);
  p3a::vector3<double> xi(0,0,0);
  HView<double***> dphi("", dim, npts, nmodes);
  for (int axis = 0; axis < dim; ++axis) {
    for (int pt = 0; pt < npts; ++pt) {
      for (int d = 0; d < dim; ++d) {
        xi[d] = pts(pt, d);
      }
      p3a::vector3<int> const deriv = p3a::vector3<int>::axis(axis);
      auto const dphi_xi = dmodes(dim, p, tensor, deriv, xi);
      for (int m = 0; m < nmodes; ++m) {
        dphi(axis, pt, m) = dphi_xi[m];
      }
    }
  }
  return dphi;
}

static HView<double**> get_phi_intr(int dim, int p, bool tensor) {
  HView<double**> const pts = get_pt_intr(dim, p);
  return get_phi(dim, p, tensor, pts);
}

static HView<double***> get_dphi_intr(int dim, int p, bool tensor) {
  HView<double**> const pts = get_pt_intr(dim, p);
  return get_dphi(dim, p, tensor, pts);
}

static HView<double**> get_phi_corner(int dim, int p, bool tensor) {
  HView<double**> const pts = get_pt_corner(dim);
  return get_phi(dim, p, tensor, pts);
}

static HView<double****> get_phi_side(int dim, int p, bool tensor) {
  HView<double****> const pts = get_pt_side(dim, p);
  int const npts = num_pts(dim-1, p);
  int const nmodes = num_modes(dim, p, tensor);
  p3a::vector3<double> xi(0,0,0);
  HView<double****> phi("", dim, ndirs, npts, nmodes);
  for (int axis = 0; axis < dim; ++axis) {
    for (int dir = 0; dir < ndirs; ++dir) {
      for (int pt = 0; pt < npts; ++pt) {
        for (int d = 0; d < dim; ++d) {
          xi[d] = pts(axis, dir, pt, d);
        }
        auto const phi_xi = modes(dim, p, tensor, xi);
        for (int m = 0; m < nmodes; ++m) {
          phi(axis, dir, pt, m) = phi_xi[m];
        }
      }
    }
  }
  return phi;
}

static HView<double***> get_phi_child_intr(int dim, int p, bool tensor) {
  HView<double***> const pts = get_pt_child_intr(dim, p);
  int const nchild = num_child(dim);
  int const npts = num_pts(dim, p);
  int const nmodes = num_modes(dim, p, tensor);
  p3a::vector3<double> xi(0,0,0);
  HView<double***> phi("", nchild, npts, nmodes);
  for (int which_child = 0; which_child < nchild; ++which_child) {
    for (int pt = 0; pt < npts; ++pt) {
      for (int d = 0; d < dim; ++d) {
        xi[d] = pts(which_child, pt, d);
      }
      auto const phi_xi = modes(dim, p, tensor, xi);
      for (int m = 0; m < nmodes; ++m) {
        phi(which_child, pt, m) = phi_xi[m];
      }
    }
  }
  return phi;
}

static HView<double*****> get_phi_child_side(int dim, int p, bool tensor) {
  HView<double*****> const pts = get_pt_child_side(dim, p);
  int const nchild = num_child(dim-1);
  int const npts = num_pts(dim-1, p);
  int const nmodes = num_modes(dim, p, tensor);
  p3a::vector3<double> xi(0,0,0);
  HView<double*****> phi("", dim, ndirs, nchild, npts, nmodes);
  for (int axis = 0; axis < dim; ++axis) {
    for (int dir = 0; dir < ndirs; ++dir) {
      for (int which_child = 0; which_child < nchild; ++which_child) {
        for (int pt = 0; pt < npts; ++pt) {
          for (int d = 0; d < dim; ++d) {
            xi[d] = pts(axis, dir, which_child, pt, d);
          }
          auto const phi_xi = modes(dim, p, tensor, xi);
          for (int m = 0; m < nmodes; ++m) {
            phi(axis, dir, which_child, pt, m) = phi_xi[m];
          }
        }
      }
    }
  }
  return phi;
}

static HView<double**> get_phi_fine(int dim, int p, bool tensor) {
  HView<double**> const pts = get_pt_fine(dim, p);
  return get_phi(dim, p, tensor, pts);
}

static HView<double**> get_phi_eval(int dim, int p, bool tensor) {
  HView<double**> const pts = get_pt_eval(dim, p);
  return get_phi(dim, p, tensor, pts);
}

static HView<double*> get_mass(int dim, int p, bool tensor) {
  HView<double*> const wts = get_wt_intr(dim, p);
  HView<double**> const phi = get_phi_intr(dim, p, tensor);
  int const npts = num_pts(dim, p);
  int const nmodes = num_modes(dim, p, tensor);
  HView<double*> mass("", nmodes);
  for (int pt = 0; pt < npts; ++pt) {
    for (int m = 0; m < nmodes; ++m) {
      mass(m) += phi(pt, m) * phi(pt, m) * wts(pt);
    }
  }
  return mass;
}

void Basis::init(int in_dim, int in_p, bool tensor_in) {
  CALI_CXX_MARK_FUNCTION;
  verify_dim(in_dim);
  verify_p(in_p);
  dim = in_dim;
  p = in_p;
  tensor = tensor_in;
  if (tensor) nmodes = num_tensor_modes(dim, p);
  else nmodes = num_non_tensor_modes(dim, p);
  copy(get_wt_intr(dim, p), wt_intr);
  copy(get_wt_side(dim, p), wt_side);
  copy(get_wt_fine(dim, p), wt_fine);
  copy(get_pt_intr(dim, p), pt_intr);
  copy(get_pt_side(dim, p), pt_side);
  copy(get_pt_child_intr(dim, p), pt_child_intr);
  copy(get_pt_child_side(dim, p), pt_child_side);
  copy(get_pt_fine(dim, p), pt_fine);
  copy(get_pt_eval(dim, p), pt_eval);
  copy(get_pt_corner(dim), pt_corner);
  copy(get_phi_intr(dim, p, tensor), phi_intr);
  copy(get_dphi_intr(dim, p, tensor), dphi_intr);
  copy(get_phi_side(dim, p, tensor), phi_side);
  copy(get_phi_child_intr(dim, p, tensor), phi_child_intr);
  copy(get_phi_child_side(dim, p, tensor), phi_child_side);
  copy(get_phi_fine(dim, p, tensor), phi_fine);
  copy(get_phi_eval(dim, p, tensor), phi_eval);
  copy(get_phi_corner(dim, p, tensor), phi_corner);
  copy(get_mass(dim, p, tensor), mass);
}

void HostBasis::init(int in_dim, int in_p, bool tensor_in) {
  CALI_CXX_MARK_FUNCTION;
  verify_dim(in_dim);
  verify_p(in_p);
  dim = in_dim;
  p = in_p;
  tensor = tensor_in;
  if (tensor) nmodes = num_tensor_modes(dim, p);
  else nmodes = num_non_tensor_modes(dim, p);
  wt_intr = get_wt_intr(dim, p);
  wt_side = get_wt_side(dim, p);
  wt_fine = get_wt_fine(dim, p);
  pt_intr = get_pt_intr(dim, p);
  pt_side = get_pt_side(dim, p);
  pt_child_intr = get_pt_child_intr(dim, p);
  pt_child_side = get_pt_child_side(dim, p);
  pt_fine = get_pt_fine(dim, p);
  pt_eval = get_pt_eval(dim, p);
  pt_corner = get_pt_corner(dim);
  phi_intr = get_phi_intr(dim, p, tensor);
  dphi_intr = get_dphi_intr(dim, p, tensor);
  phi_side = get_phi_side(dim, p, tensor);
  phi_child_intr = get_phi_child_intr(dim, p, tensor);
  phi_child_side = get_phi_child_side(dim, p, tensor);
  phi_fine = get_phi_fine(dim, p, tensor);
  phi_eval = get_phi_eval(dim, p, tensor);
  phi_corner = get_phi_corner(dim, p, tensor);
  mass = get_mass(dim, p, tensor);
}

static std::string mode_comp_name(int axis, int p) {
  if (p == 0) return "1";
  std::string c = "xi_" + std::to_string(axis);
  if (p > 1) c += "^" + std::to_string(p);
  return c;
}

void print_modal_ordering(int dim, int p, bool tensor) {
  p3a::vector3<int> const bounds = tensor_bounds(dim, p);
  for (int block = 0; block < p + 1; ++block) {
    for (int deg = 0; deg < dim * p + 1; ++deg) {
      for (int k = 0; k < bounds.z(); ++k) {
        for (int j = 0; j < bounds.y(); ++j) {
          for (int i = 0; i < bounds.x(); ++i) {
            int const sum = i + j + k;
            int const idx = p3a::max(i, p3a::max(j, k));
            if ((!tensor) && (sum > p)) continue;
            if ((idx == block) && (sum == deg)) {
              if (dim > 0) std::cout        << mode_comp_name(X, i);
              if (dim > 1) std::cout << "*" << mode_comp_name(Y, j);
              if (dim > 2) std::cout << "*" << mode_comp_name(Z, k);
              std::cout << "\n";
            }
          }
        }
      }
    }
  }
}

}
