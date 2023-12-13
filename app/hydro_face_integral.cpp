#include <dgt_array.hpp>
#include <dgt_cartesian.hpp>
#include <dgt_for_each.hpp>

#include "app.hpp"
#include "hydro.hpp"

namespace app {
namespace hydro {

void compute_face_integral(State* state)
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
              real const F_eq = F[axis][block](face, pt, eq);
              R[block](cell, eq, mode) -= sgn * F_eq * phi * detJ * wt;
            }
          }
        }
      }
    }
  };
  for_each("face_integral", num_blocks, owned_cells, functor);
}

}
}
