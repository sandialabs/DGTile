#include <dgt_array.hpp>
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

DGT_METHOD inline real get_eint(Vec<real, NEQ> const& U)
{
  real const E = U[ENER];
  real const rho = U[DENS];
  Vec3<real> const v = get_vec3(U, MMTM) / rho;
  real const half_v2 = 0.5 * dot(v,v);
  return E/rho - half_v2;
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

DGT_METHOD inline Vec<real, NEQ> get_physical_flux(
    Vec<real, NEQ> const& U,
    real const p,
    int const j)
{
  Vec<real, NEQ> F;
  real const rho = U[DENS];
  real const En = U[ENER];
  Vec3<real> const v = get_vec3(U, MMTM) / rho;
  F[DENS]     = v[j] * rho;
  F[MMTM + X] = v[j] * v[X] * rho;
  F[MMTM + Y] = v[j] * v[Y] * rho;
  F[MMTM + Z] = v[j] * v[Z] * rho;
  F[ENER]     = v[j] * (En + p);
  return F;
}

DGT_METHOD inline Vec<real, NEQ> get_hllc_flux(
    Vec<real, NEQ> const U[DIRECTIONS],
    Vec<real, NEQ> const F[DIRECTIONS],
    real const p[DIRECTIONS],
    real const c[DIRECTIONS],
    int const axis)
{
  real rho[DIRECTIONS], v[DIRECTIONS], s[DIRECTIONS];
  Vec<real, NEQ> U_star[DIRECTIONS], F_star[DIRECTIONS];
  for (int dir = 0; dir < DIRECTIONS; ++dir) {
    rho[dir] = U[dir][DENS];
    v[dir] = U[dir][MMTM + axis] / rho[dir];
  }
  s[LEFT] = std::min(v[LEFT] - c[LEFT], v[RIGHT] - c[RIGHT]);
  s[RIGHT] = std::max(v[LEFT] + c[LEFT], v[RIGHT] + c[RIGHT]);
  real const term1 = p[RIGHT] - p[LEFT];
  real const term2 = rho[LEFT] * (s[LEFT] - v[LEFT]);
  real const term3 = rho[RIGHT] * (s[RIGHT] - v[RIGHT]);
  real const num = term1 + v[LEFT] * term2 - v[RIGHT] * term3;
  real const den = term2 - term3;
  real s_star = 0.;
  if (den != 0.) s_star = num / den;
  else s_star = v[RIGHT] - v[LEFT];
  int const mpar = MMTM + permute(X, axis);
  int const mtr1 = MMTM + permute(Y, axis);
  int const mtr2 = MMTM + permute(Z, axis);
  for (int dir = 0; dir < DIRECTIONS; ++dir) {
    real const fac = (s[dir] - v[dir]) / (s[dir] - s_star);
    U_star[dir][DENS] = fac * U[dir][DENS];
    U_star[dir][mpar] = fac * U[dir][DENS] * s_star;
    U_star[dir][mtr1] = fac * U[dir][mtr1];
    U_star[dir][mtr2] = fac * U[dir][mtr2];
    U_star[dir][ENER] = fac * (U[dir][ENER] + (s_star - v[dir]) *
        (rho[dir] * s_star + p[dir]/(s[dir] - v[dir])));
    for (int eq = 0; eq < NEQ; ++eq) {
      F_star[dir][eq] = F[dir][eq] + s[dir] * (U_star[dir][eq] - U[dir][eq]);
    }
  }
  if (s[LEFT] > 0.) return F[LEFT];
  else if ((s[LEFT] <= 0) && (0. < s_star)) return F_star[LEFT];
  else if ((s_star <= 0.) && (0. <= s[RIGHT])) return F_star[RIGHT];
  else if (s[RIGHT] < 0.) return F[RIGHT];
  else return Vec<real, NEQ>::zero();
}

Hydro::Hydro(State* state_in, Input const* input_in)
{
  m_state = state_in;
  m_input = input_in;
  int const num_stored_solns = m_state->integrator->required_containers();
  m_state->mesh.add_modal({"hydro", num_stored_solns, NEQ, true});
}

static void apply_initial_conditions(Input const* input, State* state)
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

static double compute_time_step(Input const* input, State* state)
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

static void zero_residual(State* state)
{
  Mesh& mesh = state->mesh;
  int const num_blocks = mesh.num_owned_blocks();
  int const num_modes = mesh.basis().num_modes;
  Grid3 const cell_grid = mesh.cell_grid();
  auto R = mesh.get_residual("hydro").get();
  auto functor = [=] DGT_DEVICE (
      int const block,
      Vec3<int> const& cell_ijk) DGT_ALWAYS_INLINE
  {
    int const cell = cell_grid.index(cell_ijk);
    for (int eq = 0; eq < NEQ; ++eq) {
      for (int mode = 0; mode < num_modes; ++mode) {
        R[block](cell, eq, mode) = 0.;
      }
    }
  };
  for_each("zero_residual", num_blocks, cell_grid, functor);
}

static void compute_fluxes(
    State* state,
    Input const* input,
    int const soln_idx,
    int const axis)
{
  Mesh& mesh = state->mesh;
  int const num_blocks = mesh.num_owned_blocks();
  Grid3 const cell_grid = mesh.cell_grid();
  Grid3 const face_grid = get_face_grid(cell_grid, axis);
  Subgrid3 const owned_faces = get_owned_faces(cell_grid, axis);
  real const gamma = input->hydro.gamma;
  auto const B = mesh.basis();
  auto const soln = mesh.get_solution("hydro", soln_idx).get();
  auto const flux = mesh.get_fluxes("hydro", axis).get();
  int const min_loc = basis_locations::face(axis, LEFT);
  int const max_loc = basis_locations::face(axis, RIGHT);
  int const locs[DIRECTIONS] = {max_loc, min_loc};
  auto functor = [=] DGT_DEVICE (
      int const block,
      Vec3<int> const& face_ijk) DGT_ALWAYS_INLINE
  {
    real p[DIRECTIONS], c[DIRECTIONS];
    Vec<real, NEQ> U[DIRECTIONS], F[DIRECTIONS], F_hat;
    int const face = face_grid.index(face_ijk);
    for (int pt = 0; pt < B.num_face_pts; ++pt) {
      for (int dir = 0; dir < DIRECTIONS; ++dir) {
        Vec3<int> const cell_ijk = get_faces_adj_cell(face_ijk, axis, dir);
        int const cell = cell_grid.index(cell_ijk);
        U[dir] = eval<NEQ>(soln, block, cell, B, locs[dir], pt);
        real const e = get_eint(U[dir]);
        p[dir] = p_from_rho_e(U[dir][DENS], e, gamma);
        c[dir] = c_from_rho_p(U[dir][DENS], p[dir], gamma);
        F[dir] = get_physical_flux(U[dir], p[dir], axis);
      }
      F_hat = get_hllc_flux(U, F, p, c, axis);
      for (int eq = 0; eq < NEQ; ++eq) {
        flux[block](face, pt, eq) = F_hat[eq];
      }
    }
  };
  std::string const name = fmt::format("fluxes{}", get_axis_name(axis));
  for_each(name, num_blocks, owned_faces, functor);
}

static void compute_fluxes(
    State* state,
    Input const* input,
    int const soln_idx)
{
  for (int axis = 0; axis < state->mesh.dim(); ++axis) {
    compute_fluxes(state, input, soln_idx, axis);
  }
}

static void advance_explicitly(
    State* state,
    int const from_idx,
    int const to_idx,
    real const dt)
{
  Mesh& mesh = state->mesh;
  int const num_blocks = mesh.num_owned_blocks();
  Grid3 const cell_grid = mesh.cell_grid();
  Subgrid3 const owned_cells = get_owned_cells(cell_grid);
  auto const B = mesh.basis();
  auto const block_info = mesh.block_info();
  auto const R = mesh.get_residual("hydro").get();
  auto const from = mesh.get_solution("hydro", from_idx).get();
  auto const to = mesh.get_solution("hydro", to_idx).get();
  auto functor = [=] DGT_DEVICE (
      int const block,
      Vec3<int> const& cell_ijk) DGT_ALWAYS_INLINE
  {
    int const cell = cell_grid.index(cell_ijk);
    real const detJ = block_info.cell_detJs[block];
    for (int mode = 0; mode < B.num_modes; ++mode) {
      real const mass = detJ * B.mass[mode];
      real const fac = dt / mass;
      for (int eq = 0; eq < NEQ; ++eq) {
        real const from_eq = from[block](cell, eq, mode);
        real const R_eq = R[block](cell, eq, mode);
        to[block](cell, eq, mode) = from_eq + fac * R_eq;
      }
    }
  };
  for_each("hydro::advance_explicitly",
      num_blocks, owned_cells, functor);
}

static void compute_volume_integral(
    State* state,
    Input const* input,
    int const soln_idx)
{
  Mesh& mesh = state->mesh;
  int const num_blocks = mesh.num_owned_blocks();
  int const dim = mesh.dim();
  int const loc = basis_locations::CELL;
  Grid3 const cell_grid = mesh.cell_grid();
  Subgrid3 const owned_cells = get_owned_cells(cell_grid);
  real const gamma = input->hydro.gamma;
  auto const B = mesh.basis();
  auto const block_info = mesh.block_info();
  auto const soln = mesh.get_solution("hydro", soln_idx).get();
  auto const R = mesh.get_residual("hydro").get();
  auto functor = [=] DGT_DEVICE (
      int const block,
      Vec3<int> const& cell_ijk) DGT_ALWAYS_INLINE
  {
    Vec<real, NEQ> U, F;
    int const cell = cell_grid.index(cell_ijk);
    real const detJ = block_info.cell_detJs[block];
    Vec3<real> const dx = block_info.cell_dxs[block];
    for (int pt = 0; pt < B.num_cell_pts; ++pt) {
      real const wt = B.cell_weights[pt];
      U = eval<NEQ>(soln, block, cell, B, loc, pt);
      real const e = get_eint(U);
      real const p = p_from_rho_e(U[DENS], e, gamma);
      for (int axis = 0; axis < dim; ++axis) {
        F = get_physical_flux(U, p, axis);
        for (int eq = 0; eq < NEQ; ++eq) {
          for (int mode = 0; mode < B.num_modes; ++mode) {
            real const dphi_dxi = B.modes[loc].grad_phis(axis, pt, mode);
            real const dphi_dx = dphi_dxi * (2./dx[axis]);
            R[block](cell, eq, mode) += F[eq] * dphi_dx * detJ * wt;
          }
        }

      }
    }
  };
  for_each("volume_integral", num_blocks, owned_cells, functor);
}

static void compute_face_integral(State* state)
{
  Mesh& mesh = state->mesh;
  int const num_blocks = mesh.num_owned_blocks();
  int const dim = mesh.dim();
  Grid3 const cell_grid = mesh.cell_grid();
  Subgrid3 const owned_cells = get_owned_cells(cell_grid);
  auto const B = mesh.basis();
  auto const block_info = mesh.block_info();
  auto const R = mesh.get_residual("hydro").get();
  Array<Field<real***>::accessor_t, DIMENSIONS> F;
  Array<Grid3, DIMENSIONS> face_grids;
  for (int axis = 0; axis < dim; ++axis) {
    F[axis] = mesh.get_fluxes("hydro", axis).get();
    face_grids[axis] = get_face_grid(cell_grid, axis);
  }
  auto functor = [=] DGT_DEVICE(
      int const block,
      Vec3<int> const& cell_ijk) DGT_ALWAYS_INLINE
  {
    int const cell = cell_grid.index(cell_ijk);
    for (int axis = 0; axis < dim; ++axis) {
      real const detJ = block_info.face_detJs[axis][block];
      for (int dir = 0; dir < DIRECTIONS; ++dir) {
        int const loc = basis_locations::face(axis, dir);
        Vec3<int> const face_ijk = get_cells_adj_face(cell_ijk, axis, dir);
        int const face = face_grids[axis].index(face_ijk);
        real const sgn = get_dir_sign(dir);
        for (int pt = 0; pt < B.num_face_pts; ++pt) {
          real const wt = B.face_weights[pt];
          for (int mode = 0; mode < B.num_modes; ++mode) {
            real const phi = B.modes[loc].phis(pt, mode);
            for (int eq = 0; eq < NEQ; ++eq) {
              R[block](cell, eq, mode) -=
                sgn * F[axis][block](face, pt, eq) * phi * wt * detJ;
            }
          }
        }
      }
    }
  };
  for_each("face_integral", num_blocks, owned_cells, functor);
}

static void axpby(
    State* state,
    Field<real***>& r,
    real const a,
    Field<real***> const& x,
    real const b,
    Field<real***> const& y)
{
  auto R = r.get();
  auto X = x.get();
  auto Y = y.get();
  Mesh const& mesh = state->mesh;
  Grid3 const cell_grid = mesh.cell_grid();
  Subgrid3 const owned_cells = get_owned_cells(cell_grid);
  int const num_blocks = mesh.num_owned_blocks();
  int const num_modes = mesh.basis().num_modes;
  auto functor = [=] DGT_DEVICE(
      int const block,
      Vec3<int> const& cell_ijk) DGT_ALWAYS_INLINE
  {
    int const cell = cell_grid.index(cell_ijk);
    for (int eq = 0; eq < NEQ; ++eq) {
      for (int mode = 0; mode < num_modes; ++mode) {
        R[block](cell, eq, mode) =
          a * X[block](cell, eq, mode) +
          b * Y[block](cell, eq, mode);
      }
    }
  };
  for_each("axpby", num_blocks, owned_cells, functor);
}

void Hydro::apply_initial_conditions()
{
  app::apply_initial_conditions(m_input, m_state);
}

double Hydro::compute_time_step()
{
  return app::compute_time_step(m_input, m_state);
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
    int const,
    real const,
    real const)
{
  zero_residual(m_state);
  m_state->mesh.ghost("hydro", from);
  compute_fluxes(m_state, m_input, from);
  compute_volume_integral(m_state, m_input, from);
  compute_face_integral(m_state);
}

void Hydro::advance_explicitly(
    int const from,
    int const,
    int const into,
    real const dt,
    real const)
{
  app::advance_explicitly(m_state, from, into, dt);
}

void Hydro::axpby(
    int const r,
    real const a,
    int const x,
    real const b,
    int const y)
{
  Mesh& mesh = m_state->mesh;
  auto& R = mesh.get_solution("hydro", r);
  auto const& X = mesh.get_solution("hydro", x);
  auto const& Y = mesh.get_solution("hydro", y);
  app::axpby(m_state, R, a, X, b, Y);
}

}
