#include <dgt_cartesian.hpp>
#include <dgt_for_each.hpp>
#include <dgt_reduce.hpp>

#include "app.hpp"
#include "hydro.hpp"

namespace app {

enum {DENS = 0, MMTM = 1, ENER = 4, NEQ = 5};

template <class FieldT>
DGT_METHOD inline Vec<real, NEQ> get_avg(
    FieldT const& U,
    int const block,
    int const cell)
{
  Vec<real, NEQ> U_avg;
  for (int eq = 0; eq < NEQ; ++eq) {
    U_avg[eq] = U[block](cell, eq, 0);
  }
  return U_avg;
}

DGT_METHOD inline Vec3<real> get_vec3(
    Vec<real, NEQ> const& U,
    int const eq)
{
  return Vec3<real>(U[eq + X], U[eq + Y], U[eq + Z]);
}

DGT_METHOD real e_from_rho_p(
    real const rho,
    real const p,
    real const gamma)
{
  return p/(rho*(gamma-1.));
}

DGT_METHOD real p_from_rho_e(
    real const rho,
    real const e,
    real const gamma)
{
  return rho*e*(gamma-1.);
}

DGT_METHOD real c_from_rho_p(
    real const rho,
    real const p,
    real const gamma)
{
  return std::sqrt(gamma*p/rho);
}

Hydro::Hydro(State* state_in, Input const* input_in)
{
  state = state_in;
  input = input_in;
  int const num_stored_solns = state->integrator->required_containers();
  state->mesh.add_modal({"hydro", num_stored_solns, NEQ, true});
}

void Hydro::apply_initial_conditions()
{
  static constexpr int CELL = basis_locations::CELL;
  Mesh& mesh = state->mesh;
  int const nblocks = mesh.num_owned_blocks();
  Grid3 const cell_grid = mesh.cell_grid();
  Subgrid3 const owned_cells = get_owned_cells(cell_grid);
  int const ncells = generalize(mesh.dim(), cell_grid).size();
  real const gamma = input->hydro.gamma;
  inputs::Hydro const& ics = input->hydro;
  Basis<HostView> const& B = mesh.basis_h();
  BlockInfo<HostView> const& info = mesh.block_info_h();
  HostView<real***> U_host("U_ic", ncells, NEQ, B.num_modes);
  Field<real***> U_field = mesh.get_solution("hydro", 0);
  for (int block = 0; block < nblocks; ++block) {
    Kokkos::deep_copy(U_host, 0.);
    Vec3<real> const origin = info.domains[block].lower();
    Vec3<real> const dx = info.cell_dxs[block];
    auto functor = [&] (Vec3<int> const& cell_ijk) {
      int const cell = cell_grid.index(cell_ijk);
      for (int pt = 0; pt < B.num_cell_pts; ++pt) {
        real const wt = B.cell_weights(pt);
        Vec3<real> const xi = get_point(B, CELL, pt);
        Vec3<real> const x = map_to_physical(cell_ijk, origin, dx, xi);
        real const rho = ics.density->operator()(x);
        real const p = ics.pressure->operator()(x);
        Vec3<real> const v = ics.velocity->operator()(x);
        Vec3<real> const mmtm = rho * v;
        real const e = e_from_rho_p(rho, p, gamma);
        real const half_v2 = 0.5 * dot(v,v);
        real const En = rho*e + rho * half_v2;
        for (int mode = 0; mode < B.num_modes; ++mode) {
          real const phi = B.modes[CELL].phis(pt, mode);
          real const M_inv = 1./B.mass(mode);
          U_host(cell, DENS, mode) += rho * phi * wt * M_inv;
          U_host(cell, MMTM + X, mode) += mmtm.x() * phi * wt * M_inv;
          U_host(cell, MMTM + Y, mode) += mmtm.y() * phi * wt * M_inv;
          U_host(cell, MMTM + Z, mode) += mmtm.z() * phi * wt * M_inv;
          U_host(cell, ENER, mode) += En * phi * wt * M_inv;
        }
      }
    };
    seq_for_each(owned_cells, functor);
    Kokkos::deep_copy(U_field.get_view(block), U_host);
  }
}

double Hydro::compute_time_step()
{
  Mesh const& mesh = state->mesh;
  int const dim = mesh.dim();
  int const nblocks = mesh.num_owned_blocks();
  Grid3 const cell_grid = mesh.cell_grid();
  Subgrid3 const owned_cells = get_owned_cells(cell_grid);
  BlockInfo<View> const blocks = mesh.block_info();
  real const factor = 2*mesh.basis().p + 1;
  real const gamma = input->hydro.gamma;
  auto const U_field = mesh.get_solution("hydro", 0).get();
  auto functor = [=] DGT_HOST_DEVICE (
      int const block,
      Vec3<int> const& cell_ijk,
      real& dt) DGT_ALWAYS_INLINE
  {
    int const cell = cell_grid.index(cell_ijk);
    Vec3<real> const dx = blocks.cell_dxs[block];
    Vec<real, NEQ> const U = get_avg(U_field, block, cell);
    Vec3<real> const v = get_vec3(U, MMTM) / U[DENS];
    real const rho = U[DENS];
    real const En = U[ENER];
    real const e = En/rho - 0.5 * dot(v,v);
    real const p = p_from_rho_e(rho, e, gamma);
    real const c = c_from_rho_p(rho, p, gamma);
    real dvdx = 0.;
    for (int axis = 0; axis < dim; ++axis) {
      dvdx += (std::abs(v[axis]) + c) / dx[axis];
    }
    real const cell_dt = 1./(factor*dvdx);
    dt = std::min(dt, cell_dt);
  };
  real dt = DBL_MAX;
  reduce_for_each<real>("dt",
      nblocks, owned_cells, functor, Kokkos::Min<real>(dt));
  real const cfl = input->time.cfl;
  real sim_dt = cfl*dt;
  return sim_dt;
}

void Hydro::handle_visualization()
{
}

void Hydro::handle_history()
{
}

void Hydro::compute_explicit_residual(
    int const from,
    int const,
    int const into,
    real const dt,
    real const)
{
  (void)from;
  (void)into;
  (void)dt;
}

void Hydro::advance_explicitly(
    int const from,
    int const,
    int const into,
    real const dt,
    real const)
{
  (void)from;
  (void)into;
  (void)dt;
}

void Hydro::axpby(
    int const r,
    real const a,
    int const x,
    real const b,
    int const y)
{
  (void)r;
  (void)a;
  (void)x;
  (void)b;
  (void)y;
}

}
