#include <dgt_cartesian.hpp>
#include <dgt_for_each.hpp>

#include "app.hpp"
#include "hydro.hpp"

namespace app {
namespace passive {

void apply_ics(Input const* input, State* state) {
  static constexpr int CELL = basis_locations::CELL;
  Mesh& mesh = state->mesh;
  int const nblocks = mesh.num_owned_blocks();
  Grid3 const cell_grid = mesh.cell_grid();
  Subgrid3 const owned_cells = get_owned_cells(cell_grid);
  int const ncells = generalize(mesh.dim(), cell_grid).size();
  inputs::Passive const& ics = input->passive.value();
  Basis<HostView> const& B = mesh.basis_h();
  int const npassive = ics.values.size();
  BlockInfo<HostView> const& info = mesh.block_info_h();
  HostView<real***> U_host("passive_ic", ncells, npassive, B.num_modes);
  Field<real***> U_field = mesh.get_solution("passive", 0);
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
        for (int eq = 0; eq < npassive; ++eq) {
          real const val = ics.values[eq]->operator()(x);
          for (int mode = 0; mode < B.num_modes; ++mode) {
            real const phi = B.modes[CELL].phis(pt, mode);
            real const M_inv = 1./B.mass(mode);
            U_host(cell, eq, mode) += val * phi * wt * M_inv;
          }
        }
      }
    };
    seq_for_each(owned_cells, functor);
    Kokkos::deep_copy(U_field.get_view(block), U_host);
  }
}

}
}
