#include <dgt_cartesian.hpp>
#include <dgt_reduce.hpp>

#include "app.hpp"
#include "hydro.hpp"

namespace app {
namespace hydro {

real compute_time_step(Input const* input, State* state)
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

}
}
