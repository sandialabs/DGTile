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
  return cfl * dt;
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

}
