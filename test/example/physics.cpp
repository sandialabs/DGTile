#include <dgt_cartesian.hpp>
#include <dgt_reduce.hpp>

#include <dgt_print.hpp> // debug

#include "example.hpp"

namespace example {

real compute_dt(Input const& in, State const& state)
{
  Mesh const& mesh = state.mesh;
  int const num_materials = in.num_materials;
  int const nblocks = mesh.num_owned_blocks();
  Grid3 const cell_grid = mesh.cell_grid();
  Subgrid3 const owned_cells = get_owned_cells(cell_grid);
  auto eos = state.eos;
  auto U = mesh.get_solution("hydro", 0).get();
  real const cfl = in.time.cfl;
  auto functor = [=] DGT_HOST_DEVICE (
      int const block,
      Vec3<int> const& cell_ijk,
      real& dt) DGT_ALWAYS_INLINE
  {
    dt = std::min(dt, block + 13.);
  };
  real dt = DBL_MAX;
  reduce_for_each<real>("dt", nblocks, owned_cells, functor, Kokkos::Min<real>(dt));

  std::cout << dt << "\n";
  return cfl * dt;
}

}
