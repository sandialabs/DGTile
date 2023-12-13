#include <dgt_cartesian.hpp>
#include <dgt_for_each.hpp>

#include "app.hpp"
#include "hydro.hpp"

namespace app {
namespace hydro {

template <class ModalT>
DGT_METHOD inline real get_min_rho(
    ModalT const U,
    Basis<View> const& B,
    int const block,
    int const cell)
{
  constexpr int EVAL = basis_locations::EVALUATION;
  real rho_min = DBL_MAX;
  for (int pt = 0; pt < B.num_eval_pts; ++pt) {
    real const rho = eval(U, block, cell, DENS, B, EVAL, pt);
    rho_min = std::min(rho_min, rho);
  }
  return rho_min;
}

void preserve_bounds(
    State* state,
    int const soln_idx)
{
  Mesh& mesh = state->mesh;
  if (mesh.basis().p == 0) return;
  constexpr int EVAL = basis_locations::EVALUATION;
  int const num_blocks = mesh.num_owned_blocks();
  Grid3 const cell_grid = mesh.cell_grid();
  Subgrid3 const owned_cells = get_owned_cells(cell_grid);
  auto const B = mesh.basis();
  auto U = mesh.get_solution("hydro", soln_idx).get();
  auto functor = [=] DGT_DEVICE (
      int const block,
      Vec3<int> const& cell_ijk) DGT_ALWAYS_INLINE
  {
    int const cell = cell_grid.index(cell_ijk);
    { // small densities
      if (U[block](cell, DENS, 0) < rho_floor) {
        U[block](cell, DENS, 0) = rho_floor;
        for (int mode = 1; mode < B.num_modes; ++mode) {
          U[block](cell, DENS, mode) = 0.;
        }
      } else {
        real const rho_avg = U[block](cell, DENS, 0);
        real const rho_min = get_min_rho(U, B, block, cell);
        if (rho_min < rho_floor) {
          real const num = std::abs(rho_avg - rho_floor);
          real const den = std::abs(rho_avg - rho_min);
          real theta = 1.;
          if (den != 0.) theta = num / den;
          else theta = 0.;
          theta = std::min(theta, 1.);
          for (int mode = 1; mode < B.num_modes; ++mode) {
            U[block](cell, DENS, mode) *= theta;
          }
        }
      }
    }
    { // small internal energies
      Vec<real, NEQ> U_avg, U_pt;
      U_avg = get_avg(U, block, cell);
      real const eint_avg = get_eint(U_avg);
      real theta_min = 1.;
      for (int pt = 0; pt < B.num_eval_pts; ++pt) {
        U_pt = eval<NEQ>(U, block, cell, B, EVAL, pt);
        real const eint = get_eint(U_pt);
        if (eint < eint_floor) {
          real const num = std::abs(eint_avg - eint_floor);
          real const den = std::abs(eint_avg - eint);
          real theta = 1.;
          if (den != 0.) theta = num / den;
          else theta = 0.;
          theta = std::min(theta, 1.);
          theta_min = std::min(theta_min, theta);
        }
      }
      for (int eq = 0; eq < NEQ; ++eq) {
        for (int mode = 1; mode < B.num_modes; ++mode) {
          U[block](cell, eq, mode) *= theta_min;
        }
      }
    }
  };
  for_each("hydro::preserve_bounds",
      num_blocks, owned_cells, functor);
}

}
}
