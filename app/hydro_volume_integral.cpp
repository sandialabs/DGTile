#include <dgt_cartesian.hpp>
#include <dgt_for_each.hpp>

#include "app.hpp"
#include "hydro.hpp"

namespace app {
namespace hydro {

void compute_volume_integral(
    Input const* input,
    State* state,
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
      real const wt = B.cell_weights(pt);
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

}
}
