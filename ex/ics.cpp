#include <stdexcept>

#include "caliper/cali.h"

#include "dgt_grid.hpp"
#include "dgt_spatial.hpp"

#include "hydro.hpp"

namespace hydro {

#define GET_SHARED_DATA \
  Basis const b = block.basis(); \
  int const nintr_pts = dgt::num_pts(b.dim, b.p); \
  p3a::grid3 const g = block.cell_grid(); \
  p3a::grid3 const cell_grid = dgt::generalize(g); \
  p3a::vector3<double> const dx = block.dx(); \
  p3a::vector3<double> const origin = block.domain().lower(); \
  View<double***> U = block.soln(0);

static void verify_advect(Block const& block) {
  int const dim = block.dim();
  for (int axis = 0; axis < dim; ++axis) {
    if (!block.mesh()->periodic()[axis]) {
      throw std::runtime_error("advect - invalid periodicity");
    }
    if ((block.mesh()->domain().lower()[axis] != 0) ||
        (block.mesh()->domain().upper()[axis] != 1)) {
      throw std::runtime_error("advect - invalid domain");
    }
  }
}

static void verify_isentropic_vortex(Block const& block) {
  int const dim = block.dim();
  if (dim != 2) {
    throw std::runtime_error("isentropic_vortex - invalid dim");
  }
  for (int axis = 0; axis < dim; ++axis) {
    if (!block.mesh()->periodic()[axis]) {
      throw std::runtime_error("isentropic_vortex - invalid periodicity");
    }
    if ((block.mesh()->domain().lower()[axis] != 0) ||
        (block.mesh()->domain().upper()[axis] != 10)) {
      throw std::runtime_error("isentropic_vortex - invalid domain");
    }
  }
}

static void verify_sod(Block const& block) {
  if ((block.mesh()->domain().lower().x() != 0) ||
      (block.mesh()->domain().upper().x() != 1)) {
    throw std::runtime_error("sod - invalid domain");
  }
}

static void verify_rt(Block const& block, int axis) {
  if ((block.dim() != 2) ||
      (block.mesh()->domain().lower().x() != 0.) ||
      (block.mesh()->domain().upper().x() != 0.5) ||
      (block.mesh()->domain().lower().y() != 0.) ||
      (block.mesh()->domain().upper().y() != 1.5)) {
    throw std::runtime_error("rt - invalid domain");
  }
  if ((!block.mesh()->periodic().x()) ||
      (block.mesh()->periodic().y()) ||
      (block.mesh()->periodic().z())) {
    throw std::runtime_error("rt - invalid periodicity");
  }
  if (axis != Y) {
    throw std::runtime_error("rt - invalid gravity axis");
  }
}

static void set_advect_ics(Block& block) {
  CALI_CXX_MARK_FUNCTION;
  GET_SHARED_DATA;
  verify_advect(block);
  p3a::vector3<double> v(0,0,0);
  if (b.dim > 0) v.x() = 1.;
  if (b.dim > 1) v.y() = 1.;
  if (b.dim > 2) v.z() = 1.;
  v /= std::sqrt(b.dim);
  double const pi = p3a::pi_value<double>();
  auto f = [=] P3A_DEVICE (p3a::vector3<int> const& cell_ijk) {
    int const cell = cell_grid.index(cell_ijk);
    for (int pt = 0; pt < nintr_pts; ++pt) {
      double const wt = b.wt_intr(pt);
      p3a::vector3<double> const xi = dgt::get_intr_pt(b, pt);
      p3a::vector3<double> const x = dgt::get_x(cell_ijk, origin, dx, xi);
      double const ksi = x.x() + x.y() + x.z();
      double const rho = 1. + 0.25 * std::sin(2. * pi * ksi);
      p3a::vector3<double> const mmtm = v * rho;
      double const rho_eint = 0.8;
      double const E = rho_eint + 0.5 * rho * dot_product(v,v);
      for (int mode = 0; mode < b.nmodes; ++mode) {
        double const phi = b.phi_intr(pt, mode);
        double const m = b.mass(mode);
        U(cell, RH, mode) += rho * phi * wt / m;
        U(cell, MX, mode) += mmtm.x() * phi * wt / m;
        U(cell, MY, mode) += mmtm.y() * phi * wt / m;
        U(cell, MZ, mode) += mmtm.z() * phi * wt / m;
        U(cell, EN, mode) += E * phi * wt / m;
      }
    }
  };
  p3a::for_each(p3a::execution::par, cell_grid, f);
}

static void set_advect_exact(State&, Block& block, View<double***> U_ex) {
  CALI_CXX_MARK_FUNCTION;
  GET_SHARED_DATA;
  verify_advect(block);
  p3a::vector3<double> v(0,0,0);
  if (b.dim > 0) v.x() = 1.;
  if (b.dim > 1) v.y() = 1.;
  if (b.dim > 2) v.z() = 1.;
  v /= std::sqrt(b.dim);
  double const pi = p3a::pi_value<double>();
  int const nfine_pts = b.wt_fine.extent(0);
  auto f = [=] P3A_DEVICE (p3a::vector3<int> const& cell_ijk) {
    int const cell = cell_grid.index(cell_ijk);
    for (int pt = 0; pt < nfine_pts; ++pt) {
      p3a::vector3<double> const xi = dgt::get_fine_pt(b, pt);
      p3a::vector3<double> const x = dgt::get_x(cell_ijk, origin, dx, xi);
      double const ksi = x.x() + x.y() + x.z();
      double const rho = 1. + 0.25 * std::sin(2. * pi * ksi);
      p3a::vector3<double> const mmtm = v * rho;
      double const rho_eint = 0.8;
      double const E = rho_eint + 0.5 * rho * dot_product(v,v);
      U_ex(cell, pt, RH) = rho;
      U_ex(cell, pt, MX) = mmtm.x();
      U_ex(cell, pt, MY) = mmtm.y();
      U_ex(cell, pt, MZ) = mmtm.z();
      U_ex(cell, pt, EN) = E;
    }
  };
  p3a::for_each(p3a::execution::par, cell_grid, f);
  (void)nintr_pts;
}

static void set_isentropic_vortex_ics(Block& block, double gamma) {
  CALI_CXX_MARK_FUNCTION;
  GET_SHARED_DATA;
  verify_isentropic_vortex(block);
  p3a::vector3<double> const c(5,5,0);
  p3a::vector3<double> const v0(1,1,0);
  double const beta = 5.;
  double const pi = p3a::pi_value<double>();
  auto f = [=] P3A_DEVICE (p3a::vector3<int> const& cell_ijk) {
    int const cell = cell_grid.index(cell_ijk);
    for (int pt = 0; pt < nintr_pts; ++pt) {
      double const wt = b.wt_intr(pt);
      p3a::vector3<double> const xi = dgt::get_intr_pt(b, pt);
      p3a::vector3<double> const x = dgt::get_x(cell_ijk, origin, dx, xi);
      double const r = magnitude(x - c);
      double const base =
        (1. - (gamma-1.)*beta*beta / (8.*gamma*pi*pi) * std::exp(1.-r*r));
      double const rho = std::pow(base, 1./(gamma-1.));
      double const P = std::pow(rho, gamma);
      p3a::vector3<double> v(0.,0.,0.);
      v.x() = -(x.y() - c.y()) * (0.5*beta/pi) * std::exp(0.5*(1.-r*r));
      v.y() =  (x.x() - c.x()) * (0.5*beta/pi) * std::exp(0.5*(1.-r*r));
      v += v0;
      p3a::vector3<double> const mmtm = rho * v;
      double const half_v2 = 0.5*dot_product(v,v);
      double const E = P/(gamma-1.) + rho*half_v2; // ideal gas
      for (int mode = 0; mode < b.nmodes; ++mode) {
        double const phi = b.phi_intr(pt, mode);
        double const m = b.mass(mode);
        U(cell, RH, mode) += rho * phi * wt / m;
        U(cell, MX, mode) += mmtm.x() * phi * wt / m;
        U(cell, MY, mode) += mmtm.y() * phi * wt / m;
        U(cell, MZ, mode) += mmtm.z() * phi * wt / m;
        U(cell, EN, mode) += E * phi * wt / m;
      }
    }
  };
  p3a::for_each(p3a::execution::par, cell_grid, f);
}

static void set_isentropic_vortex_exact(
    State& state,
    Block& block,
    View<double***> U_ex) {
  CALI_CXX_MARK_FUNCTION;
  GET_SHARED_DATA;
  verify_isentropic_vortex(block);
  p3a::vector3<double> const c(5,5,0);
  p3a::vector3<double> const v0(1,1,0);
  double const beta = 5.;
  double const pi = p3a::pi_value<double>();
  double const gamma = state.in.gamma;
  int const nfine_pts = b.wt_fine.extent(0);
  auto f = [=] P3A_DEVICE (p3a::vector3<int> const& cell_ijk) {
    int const cell = cell_grid.index(cell_ijk);
    for (int pt = 0; pt < nfine_pts; ++pt) {
      p3a::vector3<double> const xi = dgt::get_fine_pt(b, pt);
      p3a::vector3<double> const x = dgt::get_x(cell_ijk, origin, dx, xi);
      double const r = magnitude(x - c);
      double const base =
        (1. - (gamma-1.)*beta*beta / (8.*gamma*pi*pi) * std::exp(1.-r*r));
      double const rho = std::pow(base, 1./(gamma-1.));
      double const P = std::pow(rho, gamma);
      p3a::vector3<double> v(0.,0.,0.);
      v.x() = -(x.y() - c.y()) * (0.5*beta/pi) * std::exp(0.5*(1.-r*r));
      v.y() =  (x.x() - c.x()) * (0.5*beta/pi) * std::exp(0.5*(1.-r*r));
      v += v0;
      p3a::vector3<double> const mmtm = rho * v;
      double const half_v2 = 0.5*dot_product(v,v);
      double const E = P/(gamma-1.) + rho*half_v2; // ideal gas
      U_ex(cell, pt, RH) = rho;
      U_ex(cell, pt, MX) = mmtm.x();
      U_ex(cell, pt, MY) = mmtm.y();
      U_ex(cell, pt, MZ) = mmtm.z();
      U_ex(cell, pt, EN) = E;
    }
  };
  p3a::for_each(p3a::execution::par, cell_grid, f);
  (void)nintr_pts;
}




static void set_sod_ics(Block& block, double gamma) {
  CALI_CXX_MARK_FUNCTION;
  GET_SHARED_DATA;
  verify_sod(block);
  auto f = [=] P3A_DEVICE (p3a::vector3<int> const& cell_ijk) {
    int const cell = cell_grid.index(cell_ijk);
    for (int pt = 0; pt < nintr_pts; ++pt) {
      double const wt = b.wt_intr(pt);
      p3a::vector3<double> const xi = dgt::get_intr_pt(b, pt);
      p3a::vector3<double> const x = dgt::get_x(cell_ijk, origin, dx, xi);
      double rho = 0.;
      double P = 0.;
      if (x.x() < 0.5) {
        rho = 1.;
        P = 1.;
      } else {
        rho = 0.125;
        P = 0.1;
      }
      double const E = P/(gamma-1.); // ideal gas
      for (int mode = 0; mode < b.nmodes; ++mode) {
        double const phi = b.phi_intr(pt, mode);
        double const m = b.mass(mode);
        U(cell, RH, mode) += rho * phi * wt / m;
        U(cell, EN, mode) += E * phi * wt / m;
      }
    }
  };
  p3a::for_each(p3a::execution::par, cell_grid, f);
}

static void set_rt_ics(Block& block, double gamma, double gravity, int axis) {
  CALI_CXX_MARK_FUNCTION;
  GET_SHARED_DATA;
  verify_rt(block, axis);
  double const pi = p3a::pi_value<double>();
  auto f = [=] P3A_DEVICE (p3a::vector3<int> const& cell_ijk) {
    int const cell = cell_grid.index(cell_ijk);
    for (int pt = 0; pt < nintr_pts; ++pt) {
      double const wt = b.wt_intr(pt);
      p3a::vector3<double> const xi = dgt::get_intr_pt(b, pt);
      p3a::vector3<double> const x = dgt::get_x(cell_ijk, origin, dx, xi);
      double rho = (x.y() < 0.75) ? 1.0 : 2.0;
      double const P = 2.5 + gravity * rho * (x.y() - 0.75);
      double const v =
        0.0025 *
        (1.0 - std::cos(4.*pi*x.x())) *
        (1.0 - std::cos(4.*pi*x.y()));
      double const half_v2 = 0.5*v*v;
      double const E = P/(gamma-1.) + rho*half_v2; // ideal gas
      for (int mode = 0; mode < b.nmodes; ++mode) {
        double const phi = b.phi_intr(pt, mode);
        double const m = b.mass(mode);
        U(cell, RH, mode) += rho * phi * wt / m;
        U(cell, MY, mode) += rho * v * phi * wt / m;
        U(cell, EN, mode) += E * phi * wt / m;
      }
    }
  };
  p3a::for_each(p3a::execution::par, cell_grid, f);
}

void set_ics(State& state) {
  CALI_CXX_MARK_FUNCTION;
  Input const& in = state.in;
  Mesh& mesh = state.mesh;
  double const gamma = state.in.gamma;
  for (Node* leaf : mesh.owned_leaves()) {
    Block& block = leaf->block;
    if (in.ics == "advect") {
      set_advect_ics(block);
      state.in.exact_solution = set_advect_exact;
    } else if (in.ics == "isentropic_vortex") {
      set_isentropic_vortex_ics(block, gamma);
      state.in.exact_solution = set_isentropic_vortex_exact;
    } else if (in.ics == "sod") {
      set_sod_ics(block, gamma);
    } else if (in.ics == "rt") {
      int const axis = state.in.gravity_axis;
      double const g = state.in.gravity;
      set_rt_ics(block, gamma, g, axis);
    } else {
      throw std::runtime_error("invalid ics");
    }
  }
}

}
