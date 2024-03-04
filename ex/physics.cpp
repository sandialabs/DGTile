#include <stdexcept>

#include "caliper/cali.h"

#include "p3a_for_each.hpp"

#include "dgt_amr.hpp"
#include "dgt_grid.hpp"
#include "dgt_interp.hpp"
#include "dgt_interp_simd.hpp"
#include "dgt_spatial.hpp"
#include "dgt_views.hpp"

#include "hydro.hpp"

namespace hydro {

static void handle_error_code(
    State const& state,
    std::string const& function) {
  CALI_CXX_MARK_FUNCTION;
  std::int8_t host_error;
  p3a::copy(p3a::execution::par,
      state.error_code.begin(), state.error_code.end(), &host_error);
  if (host_error == 1) {
    throw std::runtime_error(
        "error_code - invalid pressure, " + function);
  }
  if (host_error == 2) {
    throw std::runtime_error(
        "error code - invalid wave speed, " + function);
  }
}

double compute_stable_time_step(State& state, Block const& block) {
  CALI_CXX_MARK_FUNCTION;
  int const dim = block.dim();
  int const p = block.basis().p;
  double const factor = 2*p+1;
  double const gamma = state.in.gamma;
  p3a::vector3<double> const dx = block.dx();
  p3a::grid3 const g = block.cell_grid();
  p3a::grid3 const cell_grid = dgt::generalize(g);
  p3a::simd_view<double***> U = block.simd_soln(0);
  volatile std::int8_t* error_ptr = state.error_code.begin();
  auto f = [=] P3A_DEVICE (
      p3a::vector3<int> const& cell_ijk,
      p3a::device_simd_mask<double> const& mask) {
    p3a::device_simd<double> c;
    int const cell = cell_grid.index(cell_ijk);
    p3a::static_vector<p3a::device_simd<double>, NEQ> const U_avg =
      dgt::gather_avg<NEQ>(U, cell, mask);
    p3a::vector3<p3a::device_simd<double>> const v_avg =
      get_vec3(U_avg, MM) / U_avg[RH];
    c = get_wave_speed(U_avg, gamma);
    if (any_of((c != c) && mask)) { *error_ptr = 2; }
    p3a::device_simd<double> dvdx = 0.;
    for (int axis = 0; axis < dim; ++axis) {
      dvdx += (abs(v_avg[axis]) + c) / dx[axis];
    }
    p3a::device_simd<double> const dt = 1./(factor*dvdx);
    return dt;
  };
  double constexpr identity_value = p3a::maximum_value<double>();
  auto constexpr binary_op = p3a::minimizer<double>();
  double const result = p3a::simd_transform_reduce(
      p3a::execution::par, cell_grid, identity_value, binary_op, f);
  handle_error_code(state, "compute_stable_time_step");
  return result;
}

void compute_intr_fluxes(
    State& state,
    Block& block,
    int axis,
    int soln_idx) {
  CALI_CXX_MARK_FUNCTION;
  int const dim = block.dim();
  p3a::grid3 const g = block.cell_grid();
  p3a::grid3 const cell_grid = dgt::generalize(g);
  p3a::grid3 const side_grid = dgt::generalize(dgt::get_side_grid(g, axis));
  p3a::subgrid3 const intr_sides = dgt::generalize(dgt::get_intr_sides(g, axis));
  Basis const b = block.basis();
  int const nside_pts = dgt::num_pts(dim-1, b.p);
  double const gamma = state.in.gamma;
  p3a::simd_view<double***> soln = block.simd_soln(soln_idx);
  p3a::simd_view<double***> fluxes = block.simd_flux(axis);
  volatile std::int8_t* error_ptr = state.error_code.begin();
  auto f = [=] P3A_DEVICE (
      p3a::vector3<int> const& side_ijk,
      p3a::device_simd_mask<double> const& mask) {
    p3a::device_simd<double> P[ndirs], c[ndirs];
    p3a::static_vector<p3a::device_simd<double>, NEQ> U[ndirs], F[ndirs], F_hllc;
    int const side = side_grid.index(side_ijk);
    for (int pt = 0; pt < nside_pts; ++pt) {
      for (int lr = 0; lr < ndirs; ++lr) {
        int const ilr = dgt::invert_dir(lr);
        p3a::vector3<int> const cell_ijk = dgt::get_sides_adj_cell(side_ijk, axis, lr);
        int const cell = cell_grid.index(cell_ijk);
        U[lr] = dgt::interp_vec_side<NEQ>(soln, b, cell, axis, ilr, pt, mask);
        P[lr] = get_pressure(U[lr], gamma);
        c[lr] = get_wave_speed(U[lr], gamma);
        F[lr] = get_flux(U[lr], P[lr], axis);
        if (any_of((P[lr] != P[lr]) && mask)) { *error_ptr = 1; }
        if (any_of((c[lr] != c[lr]) && mask)) { *error_ptr = 2; }
      }
      F_hllc = get_hllc_flux(U, F, P, c, axis);
      for (int eq = 0; eq < NEQ; ++eq) {
        fluxes.store(F_hllc[eq], side, pt, eq, mask);
      }
    }
  };
  p3a::simd_for_each<double>(p3a::execution::par, intr_sides, f);
  handle_error_code(state, "compute_intr_fluxes");
}

[[nodiscard]] P3A_ALWAYS_INLINE P3A_HOST_DEVICE inline
p3a::static_vector<double, NEQ> gather_side(View<double***> soln, int bside, int pt) {
  p3a::static_vector<double, NEQ> U;
  for (int eq = 0; eq < NEQ; ++eq) {
    U[eq] = soln(bside, pt, eq);
  }
  return U;
}

void compute_border_fluxes(
    State& state,
    Block& block,
    int axis,
    int dir) {
  CALI_CXX_MARK_FUNCTION;
  int const dim = block.dim();
  p3a::grid3 const g = block.cell_grid();
  p3a::grid3 const side_grid = dgt::generalize(dgt::get_side_grid(g, axis));
  p3a::subgrid3 const border_sides = dgt::generalize(dgt::get_adj_sides(g, axis, dir));
  p3a::grid3 const bside_grid(border_sides.extents());
  Basis const b = block.basis();
  int const nside_pts = dgt::num_pts(dim-1, b.p);
  double const gamma = state.in.gamma;
  Border const& border = block.border(axis, dir);
  p3a::static_array<View<double***>, ndirs> bsoln = border.soln();
  View<double***> fluxes = block.flux(axis);
  volatile std::int8_t* error_ptr = state.error_code.begin();
  auto f = [=] P3A_DEVICE (p3a::vector3<int> const& side_ijk) {
    double P[ndirs], c[ndirs];
    p3a::static_vector<double, NEQ> U[ndirs], F[ndirs], F_hllc;
    p3a::vector3<int> const bside_ijk = dgt::get_border_ijk(side_ijk, axis);
    int const side = side_grid.index(side_ijk);
    int const bside = bside_grid.index(bside_ijk);
    for (int pt = 0; pt < nside_pts; ++pt) {
      for (int lr = 0; lr < ndirs; ++lr) {
        U[lr] = gather_side(bsoln[lr], bside, pt);
        P[lr] = get_pressure(U[lr], gamma);
        c[lr] = get_wave_speed(U[lr], gamma);
        F[lr] = get_flux(U[lr], P[lr], axis);
        if (P[lr] != P[lr]) { *error_ptr = 1; }
        if (c[lr] != c[lr]) { *error_ptr = 2; }
      }
      F_hllc = get_hllc_flux(U, F, P, c, axis);
      for (int eq = 0; eq < NEQ; ++eq) {
        fluxes(side, pt, eq) = F_hllc[eq];
      }
    }
  };
  p3a::for_each(p3a::execution::par, border_sides, f);
  handle_error_code(state, "compute_border_fluxes");
}

[[nodiscard]] P3A_ALWAYS_INLINE P3A_HOST_DEVICE inline
p3a::static_vector<double, NEQ> gather_amr_side(
    View<double****> soln, int bside, int child, int pt) {
  p3a::static_vector<double, NEQ> U;
  for (int eq = 0; eq < NEQ; ++eq) {
    U[eq] = soln(bside, child, pt, eq);
  }
  return U;
}

void compute_amr_border_fluxes(
    State& state,
    Block& block,
    int axis,
    int dir) {
  CALI_CXX_MARK_FUNCTION;
  int const dim = block.dim();
  p3a::grid3 const g = block.cell_grid();
  p3a::subgrid3 const border_sides = dgt::generalize(dgt::get_adj_sides(g, axis, dir));
  p3a::grid3 const bside_grid(border_sides.extents());
  Basis const b = block.basis();
  int const nside_pts = dgt::num_pts(dim-1, b.p);
  int const nchild_sides = dgt::num_child(dim-1);
  double const gamma = state.in.gamma;
  Border const& border = block.border(axis, dir);
  p3a::static_array<View<double****>, ndirs> bsoln = border.amr_soln();
  View<double****> fluxes = border.amr_flux();
  volatile std::int8_t* error_ptr = state.error_code.begin();
  auto f = [=] P3A_DEVICE (p3a::vector3<int> const& side_ijk) {
    double P[ndirs], c[ndirs];
    p3a::static_vector<double, NEQ> U[ndirs], F[ndirs], F_hllc;
    p3a::vector3<int> const bside_ijk = dgt::get_border_ijk(side_ijk, axis);
    int const bside = bside_grid.index(bside_ijk);
    for (int child = 0; child < nchild_sides; ++child) {
      for (int pt = 0; pt < nside_pts; ++pt) {
        for (int lr = 0; lr < ndirs; ++lr) {
          U[lr] = gather_amr_side(bsoln[lr], bside, child, pt);
          P[lr] = get_pressure(U[lr], gamma);
          c[lr] = get_wave_speed(U[lr], gamma);
          F[lr] = get_flux(U[lr], P[lr], axis);
          if (P[lr] != P[lr]) { *error_ptr = 1; }
          if (c[lr] != c[lr]) { *error_ptr = 2; }
        }
        F_hllc = get_hllc_flux(U, F, P, c, axis);
        for (int eq = 0; eq < NEQ; ++eq) {
          fluxes(bside, child, pt, eq) = F_hllc[eq];
        }
      }
    }
  };
  p3a::for_each(p3a::execution::par, border_sides, f);
  handle_error_code(state, "compute_amr_border_fluxes");
}

void compute_vol_integral(State& state, Block& block, int soln_idx) {
  CALI_CXX_MARK_FUNCTION;
  int const dim = block.dim();
  p3a::grid3 const g = block.cell_grid();
  p3a::grid3 const cell_grid = dgt::generalize(g);
  p3a::vector3<double> const dx = block.dx();
  double const detJ = block.cell_detJ();
  Basis const b = block.basis();
  double const gamma = state.in.gamma;
  int const nintr_pts = dgt::num_pts(dim, b.p);
  p3a::simd_view<double***> R = block.simd_resid();
  p3a::simd_view<double***> soln = block.simd_soln(soln_idx);
  volatile std::int8_t* error_ptr = state.error_code.begin();
  auto f = [=] P3A_DEVICE (
      p3a::vector3<int> const& cell_ijk,
      p3a::device_simd_mask<double> const& mask) {
    p3a::device_simd<double> P, val;
    p3a::static_vector<p3a::device_simd<double>, NEQ> U, F;
    int const cell = cell_grid.index(cell_ijk);
    for (int pt = 0; pt < nintr_pts; ++pt) {
      double const wt = b.wt_intr(pt);
      U = dgt::interp_vec_intr<NEQ>(soln, b, cell, pt, mask);
      //for (int axis = 0; axis < b.dim; ++axis) {
      //   dUdx[axis] = dgt::d_interp_vec_intr<NEQ>(soln, b, dx, axis, cell, pt, mask);
      //}
      P = get_pressure(U, gamma);
      if (any_of((P != P) && mask)) { *error_ptr = 1; }
      for (int axis = 0; axis < dim; ++axis) {
        F = get_flux(U, P, axis);
        for (int eq = 0; eq < NEQ; ++eq) {
          for (int m = 0; m < b.nmodes; ++m) {
            double const dphi_dx = b.dphi_intr(axis, pt, m) * (2./dx[axis]);
            val = F[eq] * dphi_dx * detJ * wt;
            R.sum_store(val, cell, eq, m, mask);
          }
        }
      }
    }
  };
  p3a::simd_for_each<double>(p3a::execution::par, cell_grid, f);
  handle_error_code(state, "compute_volume_integral");
}

void compute_side_integral(Block& block, int axis) {
  CALI_CXX_MARK_FUNCTION;
  p3a::grid3 const g = block.cell_grid();
  p3a::grid3 const cell_grid = dgt::generalize(g);
  p3a::grid3 const side_grid = dgt::generalize(dgt::get_side_grid(g, axis));
  double const detJ = block.side_detJ(axis);
  Basis const b = block.basis();
  int const nside_pts = dgt::num_pts(block.dim()-1, b.p);
  p3a::simd_view<double***> F = block.simd_flux(axis);
  p3a::simd_view<double***> R = block.simd_resid();
  auto f = [=] P3A_DEVICE (
      p3a::vector3<int> const& cell_ijk,
      p3a::device_simd_mask<double> const& mask) {
    p3a::device_simd<double> F_eq, val;
    int const cell = cell_grid.index(cell_ijk);
    for (int dir = 0; dir < ndirs; ++dir) {
      p3a::vector3<int> const side_ijk = dgt::get_cells_adj_side(cell_ijk, axis, dir);
      int const side = side_grid.index(side_ijk);
      double const sgn = dgt::get_dir_sign(dir);
      for (int pt = 0; pt < nside_pts; ++pt) {
        double const wt = b.wt_side(pt);
        for (int m = 0; m < b.nmodes; ++m) {
          double const phi = b.phi_side(axis, dir, pt, m);
          for (int eq = 0; eq < NEQ; ++eq) {
            F_eq = F.load(side, pt, eq, mask);
            val = -sgn * F_eq * phi * detJ * wt;
            R.sum_store(val, cell, eq, m, mask);
          }
        }
      }
    }
  };
  p3a::simd_for_each<double>(p3a::execution::par, cell_grid, f);
}

void compute_amr_side_integral(Block& block, int axis, int dir) {
  CALI_CXX_MARK_FUNCTION;
  int const dim = block.dim();
  p3a::grid3 const g = block.cell_grid();
  p3a::grid3 const cell_grid = dgt::generalize(g);
  p3a::subgrid3 const ac = dgt::get_adj_cells(g, axis, dir);
  p3a::subgrid3 const as = dgt::get_adj_sides(g, axis, dir);
  p3a::subgrid3 const border_cells = dgt::generalize(ac);
  p3a::subgrid3 const border_sides = dgt::generalize(as);
  p3a::grid3 const bside_grid(border_sides.extents());
  double const detJ = block.amr_side_detJ(axis);
  Basis const b = block.basis();
  int const nside_pts = dgt::num_pts(dim-1, b.p);
  int const nchild_sides = dgt::num_child(dim-1);
  View<double****> F = block.border(axis, dir).amr_flux();
  View<double***> R = block.resid();
  auto f = [=] P3A_DEVICE (p3a::vector3<int> const& cell_ijk) {
    int const cell = cell_grid.index(cell_ijk);
    p3a::vector3<int> const side_ijk = dgt::get_cells_adj_side(cell_ijk, axis, dir);
    double const sgn = dgt::get_dir_sign(dir);
    p3a::vector3<int> const bside_ijk = dgt::get_border_ijk(side_ijk, axis);
    int const bside = bside_grid.index(bside_ijk);
    for (int child = 0; child < nchild_sides; ++child) {
      for (int pt = 0; pt < nside_pts; ++pt) {
        double const wt = b.wt_side(pt);
        for (int m = 0; m < b.nmodes; ++m) {
          double const phi = b.phi_child_side(axis, dir, child, pt, m);
          for (int eq = 0; eq < NEQ; ++eq) {
            R(cell, eq, m) -= sgn * F(bside, child, pt, eq) * phi * detJ * wt;
          }
        }
      }
    }
  };
  p3a::for_each(p3a::execution::par, border_cells, f);
}

void compute_gravity_source(Block& block, int soln_idx, double g, int axis) {
  CALI_CXX_MARK_FUNCTION;
  Basis const b = block.basis();
  int const dim = block.dim();
  int const nintr_pts = dgt::num_pts(dim, b.p);
  p3a::grid3 const cell_grid = dgt::generalize(block.cell_grid());
  double const detJ = block.cell_detJ();
  View<double***> R = block.resid();
  View<double***> U = block.soln(soln_idx);
  auto f = [=] P3A_DEVICE (p3a::vector3<int> const& cell_ijk) {
    int const cell = cell_grid.index(cell_ijk);
    for (int pt = 0; pt < nintr_pts; ++pt) {
      double const wt = b.wt_intr(pt);
      double const rho = dgt::interp_scalar_intr(U, b, cell, pt, RH);
      double const m_a = dgt::interp_scalar_intr(U, b, cell, pt, MM + axis);
      for (int m = 0; m < b.nmodes; ++m) {
        double const phi = b.phi_intr(pt, m);
        R(cell, MM + axis, m) += (rho * g) * phi * detJ * wt;
        R(cell, EN, m) += (m_a * g) * phi * detJ * wt;
      }
    }
  };
  p3a::for_each(p3a::execution::par, cell_grid, f);
}

void advance_explicitly(Block& block, int from_idx, int to_idx, double dt) {
  CALI_CXX_MARK_FUNCTION;
  p3a::grid3 const g = block.cell_grid();
  p3a::grid3 const cell_grid = dgt::generalize(g);
  double const detJ = block.cell_detJ();
  Basis const b = block.basis();
  p3a::simd_view<double***> from = block.simd_soln(from_idx);
  p3a::simd_view<double***> to = block.simd_soln(to_idx);
  p3a::simd_view<double***> R = block.simd_resid();
  auto f = [=] P3A_DEVICE (
      p3a::vector3<int> const& cell_ijk,
      p3a::device_simd_mask<double> const& mask) {
    p3a::device_simd<double> from_eq, R_eq, val;
    int const cell = cell_grid.index(cell_ijk);
    for (int m = 0; m < b.nmodes; ++m) {
      double const mass = detJ * b.mass(m);
      for (int eq = 0; eq < NEQ; ++eq) {
        from_eq = from.load(cell, eq, m, mask);
        R_eq = R.load(cell, eq, m, mask);
        val = from_eq + (dt/mass) * R_eq;
        to.store(val, cell, eq, m, mask);
      }
    }
  };
  p3a::simd_for_each<double>(p3a::execution::par, cell_grid, f);
}

void reflect_boundary(Border& border) {
  CALI_CXX_MARK_FUNCTION;
  Block const block = border.node()->block;
  Basis const b = block.basis();
  int const axis = border.axis();
  int const dir = border.dir();
  int const dim = block.dim();
  p3a::grid3 const g = block.cell_grid();
  p3a::subgrid3 const bsides = dgt::generalize(dgt::get_adj_sides(g, axis, dir));
  p3a::grid3 const bside_grid(bsides.extents());
  int const nside_pts = dgt::num_pts(dim-1, b.p);
  View<double***> bsoln = border.soln(recv).val;
  auto f = [=] P3A_DEVICE (p3a::vector3<int> const& bside_ijk) {
    int const bside = bside_grid.index(bside_ijk);
    for (int pt = 0; pt < nside_pts; ++pt) {
      bsoln(bside, pt, MM + axis) *= -1.0;
    }
  };
  p3a::for_each(p3a::execution::par, bside_grid, f);
}

double compute_tally(Block& block, int eq) {
  CALI_CXX_MARK_FUNCTION;
  p3a::grid3 const cell_grid = dgt::generalize(block.cell_grid());
  double const detJ = block.cell_detJ();
  Basis const b = block.basis();
  View<double***> U = block.soln(0);
  int const nintr_pts = dgt::num_pts(b.dim, b.p);
  auto f = [=] P3A_DEVICE (p3a::vector3<int> const& cell_ijk) P3A_NEVER_INLINE {
    double tally = 0.;
    int const cell = cell_grid.index(cell_ijk);
    for (int pt = 0; pt < nintr_pts; ++pt) {
      double const wt = b.wt_intr(pt);
      double const U_eq = dgt::interp_scalar_intr(U, b, cell, pt, eq);
      tally += U_eq * wt * detJ;
    }
    return tally;
  };
  double constexpr identity_value = p3a::zero_value<double>();
  auto constexpr binary_op = p3a::adder<double>();
  double const result = p3a::transform_reduce(
      p3a::execution::par, cell_grid, identity_value, binary_op, f);
  return result;
}

double compute_L1_error(
    Block& block,
    View<double***> U_ex,
    int eq) {
  CALI_CXX_MARK_FUNCTION;
  p3a::grid3 const g = block.cell_grid();
  p3a::grid3 const cell_grid = dgt::generalize(g);
  Basis const b = block.basis();
  double const detJ = block.cell_detJ();
  int const nfine_intr_pts = b.phi_fine.extent(0);
  View<double***> Uh = block.soln(0);
  auto f = [=] P3A_DEVICE (p3a::vector3<int> const& cell_ijk) {
    double error = 0.;
    int const cell = cell_grid.index(cell_ijk);
    for (int pt = 0; pt < nfine_intr_pts; ++pt) {
      double const wt = b.wt_fine(pt);
      double const u = U_ex(cell, pt, eq);
      double const uh = dgt::interp_scalar_fine(Uh, b, cell, pt, eq);
      error += std::abs(u-uh) * wt * detJ;
    }
    return error;
  };
  double constexpr identity_value = p3a::zero_value<double>();
  auto constexpr binary_op = p3a::adder<double>();
  double const result = p3a::transform_reduce(
      p3a::execution::par, cell_grid, identity_value, binary_op, f);
  return result;
}

double compute_L2_error(
    Block& block,
    View<double***> U_ex,
    int eq) {
  CALI_CXX_MARK_FUNCTION;
  p3a::grid3 const g = block.cell_grid();
  p3a::grid3 const cell_grid = dgt::generalize(g);
  Basis const b = block.basis();
  double const detJ = block.cell_detJ();
  int const nfine_intr_pts = b.phi_fine.extent(0);
  View<double***> Uh = block.soln(0);
  auto f = [=] P3A_DEVICE (p3a::vector3<int> const& cell_ijk) {
    double error = 0.;
    int const cell = cell_grid.index(cell_ijk);
    for (int pt = 0; pt < nfine_intr_pts; ++pt) {
      double const wt = b.wt_fine(pt);
      double const u = U_ex(cell, pt, eq);
      double const uh = dgt::interp_scalar_fine(Uh, b, cell, pt, eq);
      double const e = (u-uh)*(u-uh);
      error += e * wt * detJ;
    }
    return error;
  };
  double constexpr identity_value = p3a::zero_value<double>();
  auto constexpr binary_op = p3a::adder<double>();
  double const result = p3a::transform_reduce(
      p3a::execution::par, cell_grid, identity_value, binary_op, f);
  return result;
}

}
