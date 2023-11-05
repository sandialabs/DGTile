#include <dgt_array.hpp>
#include <dgt_cartesian.hpp>
#include <dgt_reduce.hpp>

#include "dgt_print.hpp" // debug

#include "example.hpp"

namespace example {

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

real compute_dt(Input const& in, State const& state)
{
  Mesh const& mesh = state.mesh;
  int const dim = mesh.dim();
  int const nblocks = mesh.num_owned_blocks();
  Grid3 const cell_grid = mesh.cell_grid();
  Subgrid3 const owned_cells = get_owned_cells(cell_grid);
  EoS const eos = state.eos;
  BlockInfo<View> const blocks = mesh.block_info();
  real const factor = 2*mesh.basis().p + 1;
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
    real const e = En/rho - 0.5*dot(v,v);
    real const p = eos.p_from_rho_e(rho, e);
    real const c = eos.c_from_rho_p(rho, p);
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
  real const cfl = in.time.cfl;
  real sim_dt = cfl * dt;
  sim_dt = std::min(sim_dt, in.time.end_time - state.time);
  return sim_dt;
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

void zero_residual(State& state)
{
  Mesh& mesh = state.mesh;
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

static void compute_fluxes(State& state, int const soln_idx, int const axis)
{
  Mesh& mesh = state.mesh;
  int const num_blocks = mesh.num_owned_blocks();
  Grid3 const cell_grid = mesh.cell_grid();
  Grid3 const face_grid = get_face_grid(cell_grid, axis);
  Subgrid3 const owned_faces = get_owned_faces(cell_grid, axis);
  auto const eos = state.eos;
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
        p[dir] = eos.p_from_rho_e(U[dir][DENS], e);
        c[dir] = eos.c_from_rho_p(U[dir][DENS], p[dir]);
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

void compute_fluxes(State& state, int const soln_idx)
{
  for (int axis = 0; axis < state.mesh.dim(); ++axis) {
    compute_fluxes(state, soln_idx, axis);
  }
}

void compute_volume_integral(State& state, int const soln_idx)
{
  Mesh& mesh = state.mesh;
  int const num_blocks = mesh.num_owned_blocks();
  int const dim = mesh.dim();
  int const loc = basis_locations::CELL;
  Grid3 const cell_grid = mesh.cell_grid();
  Subgrid3 const owned_cells = get_owned_cells(cell_grid);
  auto const eos = state.eos;
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
      real const p = eos.p_from_rho_e(U[DENS], e);
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

void compute_face_integral(State& state)
{
  Mesh& mesh = state.mesh;
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

void advance_explicitly(
    State& state,
    int const from_idx,
    int const to_idx,
    real const dt)
{
  Mesh& mesh = state.mesh;
  int const num_blocks = mesh.num_owned_blocks();
  Grid3 const cell_grid = mesh.cell_grid();
  Subgrid3 const owned_cells = get_owned_cells(cell_grid);
  auto const B = mesh.basis();
  auto const block_info = mesh.block_info();
  auto const R = mesh.get_residual("hydro").get();
  auto const from = mesh.get_solution("hydro", from_idx).get();
  auto to = mesh.get_solution("hydro", to_idx).get();
  auto functor = [=] DGT_DEVICE(
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
  for_each("explicit_advance", num_blocks, owned_cells, functor);
}

//TODO: this should probably be added (along with other functionality
// like zero'ing a field) to a file like field_ops.hpp
// that provides some fast mechanism to perform operations on fields
// that is independent of the discretization state
void axpby(
    State& state,
    Field<real***>& r,
    real const a,
    Field<real***> const& x,
    real const b,
    Field<real***> const& y)
{
  auto R = r.get();
  auto X = x.get();
  auto Y = y.get();
  Mesh const& mesh = state.mesh;
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

}
