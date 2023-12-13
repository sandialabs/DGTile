#include <dgt_cartesian.hpp>
#include <dgt_for_each.hpp>

#include "app.hpp"
#include "hydro.hpp"

namespace app {
namespace hydro {

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

DGT_METHOD inline Vec<real, NEQ> get_llf_flux(
    Vec<real, NEQ> const U[DIRECTIONS],
    Vec<real, NEQ> const F[DIRECTIONS],
    real const p[DIRECTIONS],
    real const c[DIRECTIONS],
    int const axis)
{
  (void)p;
  Vec<real, NEQ> F_hat;
  real c_stab = 0.;
  for (int dir = 0; dir < DIRECTIONS; ++dir) {
    real const rho = U[dir][DENS];
    real const v = U[dir][MMTM + axis] / rho;
    real const cs = std::abs(v) + c[dir];
    c_stab = std::max(c_stab, cs);
  }
  for (int eq = 0; eq < NEQ; ++eq) {
    F_hat[eq] =
      0.5*(F[RIGHT][eq] + F[LEFT][eq]) -
      0.5*c_stab*(U[RIGHT][eq] - U[LEFT][eq]);
  }
  return F_hat;
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
      if ((0)) F_hat = get_llf_flux(U, F, p, c, axis);
      for (int eq = 0; eq < NEQ; ++eq) {
        flux[block](face, pt, eq) = F_hat[eq];
      }
    }
  };
  std::string const name = fmt::format("fluxes{}", get_axis_name(axis));
  for_each(name, num_blocks, owned_faces, functor);
}

void compute_fluxes(Input const* input, State* state, int const from)
{
  for (int axis = 0; axis < state->mesh.dim(); ++axis) {
    compute_fluxes(state, input, from, axis);
  }
}

}
}
