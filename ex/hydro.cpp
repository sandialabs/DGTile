#include <iomanip>
#include <iostream>
#include <stdexcept>

#include "caliper/cali.h"

#include "dgt_amr.hpp"
#include "dgt_grid.hpp"
#include "dgt_interp.hpp"
#include "dgt_interp_simd.hpp"
#include "dgt_marks.hpp"
#include "dgt_views.hpp"

#include "hydro.hpp"

namespace hydro {

enum {L1=0,L2=1};

static std::vector<std::string> error_names = {
  "L1", "L2"
};

static std::vector<std::string> cons_var_names = {
  "rh", "mx", "my", "mz", "en"
};

static void allocate_scratch(State& state) {
  Mesh const& mesh = state.mesh;
  int const ncells = dgt::generalize(mesh.cell_grid()).size();
  int const nmodes = mesh.basis().nmodes;
  Kokkos::resize(state.scratch, ncells, NEQ, nmodes);
}

static void setup(State& state) {
  CALI_CXX_MARK_FUNCTION;
  Input const in = state.in;
  int const dim = dgt::get_dim(in.cell_grid);
  Mesh& mesh = state.mesh;
  mesh.set_comm(in.comm);
  mesh.set_domain(p3a::box3<double>(in.xmin, in.xmax));
  mesh.set_periodic(in.periodic);
  mesh.set_cell_grid(in.cell_grid);
  mesh.set_nsoln(p3a::min(in.p+1, 2));
  mesh.set_nmodal_eq(NEQ);
  mesh.set_nflux_eq(NEQ);
  mesh.add_field("test", dim-1, 1);
  mesh.init(in.block_grid, in.p, in.tensor);
  mesh.rebuild();
  do_initial_amr(state);
  mesh.verify();
  mesh.allocate();
  allocate_scratch(state);
  state.step = 0;
  state.t = 0;
  state.dt = 0;
  state.ssp_rk_stages = in.p+1;
  state.error_code.resize(1, 0);
}

static double compute_stable_time_step(State& state) {
  CALI_CXX_MARK_FUNCTION;
  double dt = p3a::maximum_value<double>();
  for (Node* leaf : state.mesh.owned_leaves()) {
    dt = p3a::min(dt, compute_stable_time_step(state, leaf->block));
  }
  mpicpp::comm* comm = state.mesh.comm();
  mpicpp::request req = comm->iallreduce(&dt, 1, mpicpp::op::min());
  req.wait();
  dt = p3a::min(state.in.tfinal - state.t, state.in.CFL * dt);
  return dt;
}

static void compute_fluxes(State& state, int soln_idx) {
  CALI_CXX_MARK_FUNCTION;
  int const dim = state.mesh.dim();
  for (Node* leaf : state.mesh.owned_leaves()) {
    Block& block = leaf->block;
    for (int axis = 0; axis < dim; ++axis) {
      compute_intr_fluxes(state, block, axis, soln_idx);
    }
    for (int axis = 0; axis < dim; ++axis) {
      for (int dir = 0; dir < ndirs; ++dir) {
        Border const& border = block.border(axis, dir);
        if (border.type() == dgt::COARSE_TO_FINE) {
          compute_amr_border_fluxes(state, block, axis, dir);
        } else {
          compute_border_fluxes(state, block, axis, dir);
        }
      }
    }
  }
}

static void zero_residual(State& state) {
  CALI_CXX_MARK_FUNCTION;
  for (Node* leaf : state.mesh.owned_leaves()) {
    Kokkos::deep_copy(leaf->block.resid(), 0.);
  }
}

static void compute_vol_integral(State& state, int soln_idx) {
  CALI_CXX_MARK_FUNCTION;
  for (Node* leaf : state.mesh.owned_leaves()) {
    compute_vol_integral(state, leaf->block, soln_idx);
  }
}

static void compute_side_integral(State& state) {
  CALI_CXX_MARK_FUNCTION;
  int const dim = state.mesh.dim();
  for (Node* leaf : state.mesh.owned_leaves()) {
    Block& block = leaf->block;
    for (int axis = 0; axis < dim; ++axis) {
      compute_side_integral(block, axis);
    }
    for (int axis = 0; axis < dim; ++axis) {
      for (int dir = 0; dir < ndirs; ++dir) {
        Border const& border = block.border(axis, dir);
        if (border.type() == dgt::COARSE_TO_FINE) {
          compute_amr_side_integral(block, axis, dir);
        }
      }
    }
  }
}

static void compute_gravity_source(State& state, int soln_idx) {
  CALI_CXX_MARK_FUNCTION;
  if (state.in.gravity == 0.) return;
  int const axis = state.in.gravity_axis;
  double const g = state.in.gravity;
  for (Node* leaf : state.mesh.owned_leaves()) {
    Block& block = leaf->block;
    compute_gravity_source(block, soln_idx, g, axis);
  }
}

static void advance_explicitly(State& state, int from_idx, int to_idx, double dt) {
  CALI_CXX_MARK_FUNCTION;
  for (Node* leaf : state.mesh.owned_leaves()) {
    advance_explicitly(leaf->block, from_idx, to_idx, dt);
  }
}

static void limit(State& state, int soln_idx) {
  CALI_CXX_MARK_FUNCTION;
  begin_border_transfer(state.mesh, soln_idx);
  end_border_transfer(state.mesh);
  for (Node* leaf : state.mesh.owned_leaves()) {
    limit(state, leaf->block, soln_idx, state.scratch);
  }
}

static void reflect_boundary(State& state) {
  CALI_CXX_MARK_FUNCTION;
  if (state.in.ics != "rt") return;
  int const dim = state.mesh.dim();
  for (Node* leaf : state.mesh.owned_leaves()) {
    Block& block = leaf->block;
    for (int axis = 0; axis < dim; ++axis) {
      for (int dir = 0; dir < ndirs; ++dir) {
        Border& border = block.border(axis, dir);
        if (border.type() == dgt::BOUNDARY) {
          reflect_boundary(border);
        }
      }
    }
  }
}

static void preserve_bounds(State& state, int soln_idx) {
  CALI_CXX_MARK_FUNCTION;
  int const dim = state.mesh.dim();
  for (Node* leaf : state.mesh.owned_leaves()) {
    Block& block = leaf->block;
    preserve_bounds(state, block, soln_idx);
    for (int axis = 0; axis < dim; ++axis) {
      for (int dir = 0; dir < ndirs; ++dir) {
        Border const& border = block.border(axis, dir);
        if (border.type() == dgt::COARSE_TO_FINE) {
          preserve_bounds_amr(state, block, axis, dir, soln_idx);
        }
      }
    }
  }
}

static void axpby(State& state, int r, double a, int x, double b, int y) {
  CALI_CXX_MARK_FUNCTION;
  for (Node* leaf : state.mesh.owned_leaves()) {
    Block& block = leaf->block;
    View<double***> Ur = block.soln(r);
    View<double***> Ux = block.soln(x);
    View<double***> Uy = block.soln(y);
    dgt::axpby(p3a::execution::par, Ur, a, Ux, b, Uy);
  }
}

static double compute_tally(State& state, int eq) {
  CALI_CXX_MARK_FUNCTION;
  double tally = 0.;
  for (Node* leaf : state.mesh.owned_leaves()) {
    Block& block = leaf->block;
    tally += compute_tally(block, eq);
  }
  mpicpp::comm* comm = state.mesh.comm();
  mpicpp::request req = comm->iallreduce(&tally, 1, mpicpp::op::sum());
  req.wait();
  return tally;
}

static void print_tallies(State& state) {
  p3a::static_vector<double, NEQ> tallies;
  for (int eq = 0; eq < NEQ; ++eq) {
    tallies[eq] = compute_tally(state, eq);
  }
  if (state.mesh.comm()->rank() != 0) return;
  std::cout << std::scientific << std::setprecision(16);
  std::cout << "--- tallies ---\n";
  for (int eq = 0; eq < NEQ; ++eq) {
    std::cout << "[" << cons_var_names[eq] << "]: " << tallies[eq] << "\n";
  }
  std::cout << "---\n";
}

static void print_step(
    mpicpp::comm* comm,
    int freq,
    int step,
    double t,
    double dt) {
  if (step % freq) return;
  if (comm->rank() == 0) {
    std::cout << std::scientific;
    std::cout << std::setprecision(7);
    std::cout << "[step] " << std::setw(7) << step;
    std::cout << " [t] " << t;
    std::cout << " [dt] " << dt;
    std::cout << "\n";
  }
}

static std::array<int, 2> ssp_idx(int nstages, int stage) {
  if (nstages == 1)     return std::array<int, 2>{0, 0};
  else if (stage == 0)  return std::array<int, 2>{0, 1};
  else                  return std::array<int, 2>{1, 1};
}

static void axpby(State& state, int nstages, int stage) {
  if (stage == 0) return;
  if (nstages == 2 && stage == 1) axpby(state, 0, 0.5, 0, 0.5, 1);
  if (nstages == 3 && stage == 1) axpby(state, 1, 0.75, 0, 0.25, 1);
  if (nstages == 3 && stage == 2) axpby(state, 0, 1./3., 0, 2./3., 1);
}


template <class F>
double compute_error(int type, F const& f, State& state, int eq) {
  CALI_CXX_MARK_FUNCTION;
  if (!state.in.exact_solution) return -1.;
  int const dim = state.mesh.dim();
  int const nfine_pts = state.mesh.basis().wt_fine.extent(0);
  int const ncells = dgt::generalize(state.mesh.cell_grid()).size();
  double const volume = dgt::get_volume(dim, state.mesh.domain().extents());
  Kokkos::resize(state.scratch, ncells, nfine_pts, NEQ);
  View<double***> U_ex = state.scratch;
  double error = 0.;
  for (Node* leaf : state.mesh.owned_leaves()) {
    Block& block = leaf->block;
    state.in.exact_solution(state, block, U_ex);
    error += f(block, U_ex, eq);
  }
  mpicpp::comm* comm = state.mesh.comm();
  mpicpp::request req = comm->iallreduce(&error, 1, mpicpp::op::sum());
  req.wait();
  if (type == L2) { error = std::sqrt(error / volume); }
  if (state.in.comm->rank() == 0) {
    std::cout << error_names[type] << " " << cons_var_names[eq]
      << " error: " << error << "\n";
  }
  return error;
}

static void check_error_regression(State& state, double computed) {
  if (!state.in.exact_solution) return;
  if (state.in.error_regression == 0.) return;
  double const expected = state.in.error_regression;
  double const absval = std::abs(expected - computed);
  if (absval > 1.e-5) {
    throw std::runtime_error("test failed!\n");
  }
}

static void run(mpicpp::comm* comm, std::string const& name) {
  CALI_CXX_MARK_FUNCTION;
  State state;
  parse_input(state.in, comm, name);
  print_input(state.in);
  verify_input(state.in);
  setup(state);
  set_ics(state);
  preserve_bounds(state, 0);
  print_tallies(state);
  while (true) {
    do_amr(state);
    write_out(state);
    if (state.t >= state.in.tfinal) break;
    state.dt = compute_stable_time_step(state);
    print_step(comm, state.in.step_frequency, state.step, state.t, state.dt);
    int const nstages = state.ssp_rk_stages;
    for (int stage = 0; stage < nstages; ++stage) {
      int const from = ssp_idx(nstages, stage)[0];
      int const to = ssp_idx(nstages, stage)[1];
      begin_border_transfer(state.mesh, from);
      end_border_transfer(state.mesh);
      reflect_boundary(state);
      zero_residual(state);
      compute_fluxes(state, from);
      compute_vol_integral(state, from);
      compute_side_integral(state);
      compute_gravity_source(state, from);
      advance_explicitly(state, from, to, state.dt);
      limit(state, to);
      preserve_bounds(state, to);
      axpby(state, nstages, stage);
    }
    state.t += state.dt;
    state.step++;
  }
  print_step(comm, 1, state.step, state.t, state.dt);
  print_tallies(state);
  write_pvd(state);
  double const L1_error = compute_error(L1, compute_L1_error, state, RH);
  double const L2_error = compute_error(L2, compute_L2_error, state, RH);
  check_error_regression(state, L2_error);
  (void)L1_error;
}

}

int main(int argc, char** argv) {
  dgt::Library library(argc, argv);
  if (argc != 2) {
    std::string exe = argv[0];
    throw std::runtime_error("usage: " + exe + " input");
  }
  mpicpp::comm comm = mpicpp::comm::world();
  hydro::run(&comm, argv[1]);
}
