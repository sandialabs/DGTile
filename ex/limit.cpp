#include <stdexcept>

#include "caliper/cali.h"

#include "dgt_amr.hpp"
#include "dgt_grid.hpp"
#include "dgt_interp.hpp"
#include "dgt_interp_simd.hpp"
#include "dgt_spatial.hpp"
#include "dgt_views.hpp"

#include "hydro.hpp"

namespace hydro {

static constexpr int center = 2;

static constexpr double rho_floor_val = 1.e-12;
static constexpr double P_floor_val = 1.e-12;

[[nodiscard]] P3A_ALWAYS_INLINE P3A_HOST_DEVICE inline
p3a::static_vector<double, NEQ> gather_dg_slopes(
    View<double***> U,
    int cell,
    int mode) {
  p3a::static_vector<double, NEQ> dc;
  for (int eq = 0; eq < NEQ; ++eq) {
    dc[eq] = U(cell, eq, mode);
  }
  return dc;
}

[[nodiscard]] P3A_ALWAYS_INLINE P3A_HOST_DEVICE inline
p3a::static_vector<double, NEQ> gather_border_avg(
    View<double**> soln,
    int const bside) {
  p3a::static_vector<double, NEQ> U_avg;
  for (int eq = 0; eq < NEQ; ++eq) {
    U_avg[eq] = soln(bside, eq);
  }
  return U_avg;
}

[[nodiscard]] P3A_ALWAYS_INLINE P3A_HOST_DEVICE inline
p3a::static_vector<double, NEQ> gather_amr_border_avg(
    View<double***> soln,
    int const nchild,
    int const bside) {
  p3a::static_vector<double, NEQ> U_avg = decltype(U_avg)::zero();
  for (int eq = 0; eq < NEQ; ++eq) {
    for (int child = 0; child < nchild; ++child) {
      U_avg[eq] += soln(bside, child, eq);
    }
    U_avg[eq] /= nchild;
  }
  return U_avg;
}

[[nodiscard]] P3A_ALWAYS_INLINE P3A_HOST_DEVICE inline
p3a::static_vector<double, NEQ> convert_to_primitive(
    p3a::static_vector<double, NEQ> const& c,
    double gamma) {
  p3a::static_vector<double, NEQ> p;
  p[RH] = c[RH];
  p[VX] = c[MX]/c[RH];
  p[VY] = c[MY]/c[RH];
  p[VZ] = c[MZ]/c[RH];
  p[PR] = get_pressure(c, gamma);
  return p;
}

[[nodiscard]] P3A_ALWAYS_INLINE P3A_HOST_DEVICE inline
p3a::static_vector<double, NEQ> convert_to_dprimitive(
    p3a::static_vector<double, NEQ> const& c,
    p3a::static_vector<double, NEQ> const& dc,
    double gamma) {
  p3a::static_vector<double, NEQ> dp;
  double const rho = c[RH];
  p3a::vector3<double> const v = get_vec3(c, MM) / rho;
  double const v2 = dot_product(v, v);
  dp[RH] = dc[RH];
  dp[VX] = (dc[MX] - dp[RH]*v.x())/rho;
  dp[VY] = (dc[MY] - dp[RH]*v.y())/rho;
  dp[VZ] = (dc[MZ] - dp[RH]*v.z())/rho;
  dp[PR] = (gamma-1.)*(
      dc[EN] -
      0.5*v2*dp[RH] -
      rho*v.x()*dp[VX] -
      rho*v.y()*dp[VY] -
      rho*v.z()*dp[VZ]); // ideal gas
  return dp;
}

[[nodiscard]] P3A_ALWAYS_INLINE P3A_HOST_DEVICE inline
p3a::static_vector<double, NEQ> convert_to_dconservative(
    p3a::static_vector<double, NEQ> const c,
    p3a::static_vector<double, NEQ> const& dp,
    double gamma) {
  p3a::static_vector<double, NEQ> dc;
  double const rho = c[RH];
  p3a::vector3<double> const v = get_vec3(c, MM) / rho;
  double const v2 = dot_product(v, v);
  dc[RH] = dp[RH];
  dc[MX] = dp[RH]*v.x() + rho*dp[VX];
  dc[MY] = dp[RH]*v.y() + rho*dp[VY];
  dc[MZ] = dp[RH]*v.z() + rho*dp[VZ];
  dc[EN] = 0.5*v2*dp[RH] +
    rho*v.x()*dp[VX] +
    rho*v.y()*dp[VY] +
    rho*v.z()*dp[VZ] +
    (1./(gamma-1.))*dp[PR]; // ideal gas
  return dc;
}

[[nodiscard]] P3A_ALWAYS_INLINE P3A_HOST_DEVICE inline
double minmod(double a, double b) {
  if (p3a::sign(a) != p3a::sign(b)) return 0.;
  else return p3a::sign(a)*p3a::min(p3a::abs(a), p3a::abs(b));
}

[[nodiscard]] P3A_ALWAYS_INLINE P3A_HOST_DEVICE inline
double minmod(double a, double b, double c) {
  return minmod(a,minmod(b,c));
}

[[nodiscard]] P3A_ALWAYS_INLINE P3A_HOST_DEVICE inline
double minmodB(
    double beta,
    double M,
    double dx,
    double dp,
    double p,
    double p_left,
    double p_right) {
  double const dp_left = beta*(p-p_left);
  double const dp_right = beta*(p_right-p);
  if (std::abs(dp) < M*dx*dx) return dp;
  else return minmod(dp, dp_left, dp_right);
}

struct BorderData {
  public:
    p3a::static_array<p3a::static_array<int, ndirs>, DIMS> border_type;
    p3a::static_array<p3a::static_array<p3a::grid3, ndirs>, DIMS> bside_grid;
    p3a::static_array<p3a::static_array<View<double**>, ndirs>, DIMS> avg_U;
    p3a::static_array<p3a::static_array<View<double***>, ndirs>, DIMS> amr_avg_U;
  public:
    BorderData(Block& block) {
      int const dim = block.dim();
      p3a::grid3 const g = block.cell_grid();
      for (int axis = 0; axis < dim; ++axis) {
        for (int dir = 0; dir < ndirs; ++dir) {
          Border& border = block.border(axis, dir);
          border_type[axis][dir] = border.type();
          p3a::subgrid3 const adj_sides = dgt::get_adj_sides(g, axis, dir);
          bside_grid[axis][dir] = dgt::generalize(adj_sides.extents());
          if (border_type[axis][dir] == dgt::COARSE_TO_FINE) {
            amr_avg_U[axis][dir] = border.amr(recv).avg_soln;
          } else {
            avg_U[axis][dir] = border.avg_soln(recv).val;
          }
        }
      }
    }
};

[[nodiscard]] P3A_ALWAYS_INLINE P3A_HOST_DEVICE inline
bool needs_border(
    p3a::grid3 const& cell_grid,
    p3a::vector3<int> const& cell_ijk,
    int axis,
    int dir) {
  int const start = 0;
  int const end = cell_grid.extents()[axis] - 1;
  if ((cell_ijk[axis] == start) && (dir == left)) return true;
  if ((cell_ijk[axis] == end) && (dir == right)) return true;
  return false;
}

[[nodiscard]] P3A_ALWAYS_INLINE P3A_HOST_DEVICE inline
p3a::static_vector<double, NEQ> get_adj_avg(
    View<double***> soln,
    BorderData const& dat,
    p3a::grid3 const& cell_grid,
    p3a::vector3<int> const& cell_ijk,
    int axis,
    int dir,
    int nchild) {
  p3a::static_vector<double, NEQ> c;
  if (!needs_border(cell_grid, cell_ijk, axis, dir)) {
    p3a::vector3<int> const adj_cell_ijk = dgt::get_cells_adj_cell(cell_ijk, axis, dir);
    int const adj_cell = cell_grid.index(adj_cell_ijk);
    c = dgt::gather_avg<NEQ>(soln, adj_cell);
  } else {
    p3a::vector3<int> const side_ijk = dgt::get_cells_adj_side(cell_ijk, axis, dir);
    p3a::vector3<int> const bside_ijk = dgt::get_border_ijk(side_ijk, axis);
    int const bside = dat.bside_grid[axis][dir].index(bside_ijk);
    if (dat.border_type[axis][dir] != dgt::COARSE_TO_FINE) {
      c = gather_border_avg(dat.avg_U[axis][dir], bside);
    } else {
      c = gather_amr_border_avg(dat.amr_avg_U[axis][dir], nchild, bside);
    }
  }
  return c;
}

void limit(State& state, Block& block, int soln_idx, View<double***> soln_lim) {
  CALI_CXX_MARK_FUNCTION;
  if (block.basis().p == 0) return;
  int const dim = block.dim();
  p3a::grid3 const g = block.cell_grid();
  p3a::grid3 const cell_grid = dgt::generalize(g);
  p3a::vector3<double> const dx = block.dx();
  Basis const b = block.basis();
  int const nchild = dgt::num_child(dim-1);
  BorderData dat(block);
  double const gamma = state.in.gamma;
  double const M = state.in.M;
  double const beta = state.in.beta;
  View<double***> soln = block.soln(soln_idx);
  Kokkos::deep_copy(soln_lim, soln);
  auto f = [=] P3A_DEVICE (p3a::vector3<int> const& cell_ijk) {
    p3a::static_vector<double, NEQ> dc, dp, dc_lim, dp_lim;
    p3a::static_vector<double, NEQ> c[ndirs+1], p[ndirs+1];
    bool should_zero = false;
    int const cell = cell_grid.index(cell_ijk);
    c[center] = dgt::gather_avg<NEQ>(soln, cell);
    p[center] = convert_to_primitive(c[center], gamma);
    for (int axis = 0; axis < dim; ++axis) {
      int const dg_mode = 1 + axis;
      dc = gather_dg_slopes(soln, cell, dg_mode);
      dp = convert_to_dprimitive(c[center], dc, gamma);
      for (int dir = 0; dir < ndirs; ++dir) {
        c[dir] = get_adj_avg(soln, dat, cell_grid, cell_ijk, axis, dir, nchild);
        p[dir] = convert_to_primitive(c[dir], gamma);
      }
      for (int eq = 0; eq < NEQ; ++eq) {
        dp_lim[eq] = minmodB(beta, M, dx[axis], dp[eq], p[center][eq], p[left][eq], p[right][eq]);
        if (dp_lim[eq] != dp[eq]) {
          should_zero = true;
        }
      }
      dc_lim = convert_to_dconservative(c[center], dp_lim, gamma);
      for (int eq = 0; eq < NEQ; ++eq) {
        soln_lim(cell, eq, dg_mode) = dc_lim[eq];
      }
    }
    if (should_zero) {
      for (int eq = 0; eq < NEQ; ++eq) {
        for (int m = dim+1; m < b.nmodes; ++m) {
          soln_lim(cell, eq, m) = 0.;
        }
      }
    }
  };
  p3a::for_each(p3a::execution::par, cell_grid, f);
  Kokkos::deep_copy(soln, soln_lim);
}

[[nodiscard]] P3A_ALWAYS_INLINE P3A_HOST_DEVICE inline
double get_min_eval(View<double***> U, Basis const& b, int cell, int eq) {
  double min_val = p3a::maximum_value<double>();
  int const neval_pts = dgt::num_eval_pts(b.dim, b.p);
  for (int pt = 0; pt < neval_pts; ++pt) {
    double const val = dgt::interp_scalar_eval(U, b, cell, pt, eq);
    min_val = p3a::min(min_val, val);
  }
  return min_val;
}

[[nodiscard]] P3A_ALWAYS_INLINE P3A_HOST_DEVICE inline
double get_min_amr(
    View<double***> U,
    Basis const& b,
    int cell,
    int axis,
    int dir,
    int eq) {
  int const nchild_sides = dgt::num_child(b.dim-1);
  int const nside_pts = dgt::num_pts(b.dim-1, b.p);
  double min_val = p3a::maximum_value<double>();
  for (int child = 0; child < nchild_sides; ++child) {
    for (int pt = 0; pt < nside_pts; ++pt) {
      double const val = dgt::interp_scalar_child_side(
          U, b, cell, axis, dir, child, pt, eq);
      min_val = p3a::min(min_val, val);
    }
  }
  return min_val;
}

[[nodiscard]] P3A_ALWAYS_INLINE P3A_HOST_DEVICE inline
double bound_theta(double theta) {
  double btheta = theta;
  btheta = p3a::min(theta, 1.);
  btheta = p3a::max(theta, 0.);
  return btheta;
}

void preserve_bounds(State& state, Block& block, int soln_idx) {
  CALI_CXX_MARK_FUNCTION;
  if (block.basis().p == 0) return;
  p3a::grid3 const g = block.cell_grid();
  p3a::grid3 const cell_grid = dgt::generalize(g);
  Basis const b = block.basis();
  int const neval_pts = dgt::num_eval_pts(b.dim, b.p);
  View<double***> U = block.soln(soln_idx);
  double const gamma = state.in.gamma;
  double const rho_floor = rho_floor_val;
  double const P_floor = P_floor_val;
  double const E_floor = P_floor/(gamma-1.); // ideal gas
  auto f = [=] P3A_DEVICE (p3a::vector3<int> const& cell_ijk) {
    int const cell = cell_grid.index(cell_ijk);
    p3a::static_vector<double, NEQ> U_pt;
    { // small densities
      if (U(cell, RH, 0) < rho_floor) {
        U(cell, RH, 0) = rho_floor;
        U(cell, MX, 0) = 0.;
        U(cell, MY, 0) = 0.;
        U(cell, MZ, 0) = 0.;
        U(cell, EN, 0) = E_floor;
        for (int eq = 0; eq < NEQ; ++eq) {
          for (int m = 1; m < b.nmodes; ++m) {
            U(cell, eq, m) = 0.;
          }
        }
      } else {
        double const rho_avg = U(cell, RH, 0);
        double const rho_min = get_min_eval(U, b, cell, RH);
        if (rho_min < rho_floor) {
          double const theta = (rho_avg - rho_floor) / (rho_avg - rho_min);
          double const bounded_theta = bound_theta(theta);
          for (int m = 1; m < b.nmodes; ++m) {
            U(cell, RH, m) *= bounded_theta;
          }
        }
      }
    }
    { // small pressures
      p3a::static_vector<double, NEQ> const U_avg = dgt::gather_avg<NEQ>(U, cell);
      double const P_avg = get_pressure(U_avg, gamma);
      if (P_avg < P_floor) {
        for (int eq = 0; eq < NEQ; ++eq) {
          for (int m = 1; m < b.nmodes; ++m) {
            U(cell, eq, m) = 0.;
          }
        }
      } else {
        double theta_min = 1.;
        for (int pt = 0; pt < neval_pts; ++pt) {
          U_pt = dgt::interp_vec_eval<NEQ>(U, b, cell, pt);
          double const P_pt = get_pressure(U_pt, gamma);
          if (P_pt < P_floor) {
            double const theta = (P_avg-P_floor)/(P_avg-P_pt);
            theta_min = p3a::min(theta_min, theta);
          }
          double const bounded_theta_min = bound_theta(theta_min);
          for (int eq = 0; eq < NEQ; ++eq) {
            for (int m = 1; m < b.nmodes; ++m) {
              U(cell, eq, m) *= bounded_theta_min;
            }
          }
        }
      }
    }
  };
  p3a::for_each(p3a::execution::par, cell_grid, f);
}

void preserve_bounds_amr(
    State& state,
    Block& block,
    int axis,
    int dir,
    int soln_idx) {
  CALI_CXX_MARK_FUNCTION;
  if (block.basis().p == 0) return;
  p3a::grid3 const g = block.cell_grid();
  p3a::grid3 const cell_grid = dgt::generalize(g);
  Basis const b = block.basis();
  View<double***> U = block.soln(soln_idx);
  p3a::subgrid3 const border_cells = dgt::generalize(dgt::get_adj_cells(g, axis, dir));
  double const gamma = state.in.gamma;
  double const rho_floor = rho_floor_val;
  double const P_floor = P_floor_val;
  int const nchild_sides = dgt::num_child(b.dim-1);
  int const nside_pts = dgt::num_pts(b.dim-1, b.p);
  auto f = [=] P3A_DEVICE (p3a::vector3<int> const& cell_ijk) {
    p3a::static_vector<double, NEQ> U_pt;
    int const cell = cell_grid.index(cell_ijk);
    { // small densities
      double const rho_avg = U(cell, RH, 0);
      double const rho_min = get_min_amr(U, b, cell, axis, dir, RH);
      if (rho_min < rho_floor) {
        if (rho_avg != rho_min) {
          double const theta = (rho_avg - rho_floor) / (rho_avg - rho_min);
          double const bounded_theta = bound_theta(theta);
          for (int m = 1; m < b.nmodes; ++m) {
            U(cell, RH, m) *= bounded_theta;
          }
        }
      }
    }
    { // small pressures
      double theta_min = 1.;
      p3a::static_vector<double, NEQ> const U_avg = dgt::gather_avg<NEQ>(U, cell);
      double const P_avg = get_pressure(U_avg, gamma);
      for (int child = 0; child < nchild_sides; ++child) {
        for (int pt = 0; pt < nside_pts; ++pt) {
          U_pt = dgt::interp_vec_child_side<NEQ>(U, b, cell, axis, dir, child, pt);
          double const P_pt = get_pressure(U_pt, gamma);
          if (P_pt < P_floor) {
            double const theta = (P_avg-P_floor)/(P_avg-P_pt);
            theta_min = p3a::min(theta_min, theta);
          }
        }
      }
      double const bounded_theta_min = bound_theta(theta_min);
      for (int eq = 0; eq < NEQ; ++eq) {
        for (int m = 1; m < b.nmodes; ++m) {
          U(cell, eq, m) *= bounded_theta_min;
        }
      }
    }
  };
  p3a::for_each(p3a::execution::par, border_cells, f);
}

}
