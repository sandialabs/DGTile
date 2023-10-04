#include <dgt_cartesian.hpp>
#include <dgt_reduce.hpp>

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
  auto const U = mesh.get_solution("hydro", soln_idx).get();
  auto const F = mesh.get_fluxes("hydro", axis).get();
  auto functor = [=] DGT_DEVICE (
      int const block,
      Vec3<int> const& face_ijk) DGT_ALWAYS_INLINE
  {
    int const face = face_grid.index(face_ijk);
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
