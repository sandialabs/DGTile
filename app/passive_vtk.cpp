#include <filesystem>

#include <dgt_cartesian.hpp>
#include <dgt_for_each.hpp>
#include <dgt_print.hpp>
#include <dgt_vtk.hpp>

#include "app.hpp"
#include "hydro.hpp"

namespace app {
namespace passive {

struct Data
{
  typename Field<real***>::accessor_t U;
  Basis<View> B;
  real gamma;
  int block;
};

dgt::vtk::VtkView<real> get_variable(
    Input const* input,
    State const* state,
    int const block,
    int const soln_idx,
    int const eq)
{
  static constexpr int CELL_LOC = basis_locations::CELL;
  Mesh const& mesh = state->mesh;
  Basis<View> const& B = mesh.basis();
  Grid3 const cell_grid = mesh.cell_grid();
  Grid3 const inner_grid = tensor_bounds(B.dim, B.q-1);
  Grid3 const ginner_grid = generalize(B.dim, inner_grid);
  Grid3 const viz_cell_grid = vtk::get_viz_cell_grid(cell_grid, B.q);
  Grid3 const gviz_cell_grid = generalize(B.dim, viz_cell_grid);
  Subgrid3 const owned_cells = get_owned_cells(cell_grid);
  Vec3<int> const ghost_offset = owned_cells.lower();
  real const gamma = input->hydro.gamma;
  auto const U = mesh.get_solution("passive", soln_idx).get();
  dgt::vtk::VtkView<real> var;
  Kokkos::resize(var, gviz_cell_grid.size(), 1);
  Data data = {U, B, gamma, block};
  auto functor = [=] DGT_HOST_DEVICE (Vec3<int> const& cell_ijk) {
    int const cell = cell_grid.index(cell_ijk);
    Vec3<int> const owned_ijk = cell_ijk - ghost_offset;
    inner_for_each(ginner_grid,
    [&] (Vec3<int> const& inner_ijk) DGT_ALWAYS_INLINE {
      int const pt = ginner_grid.index(inner_ijk);
      Vec3<int> const viz_cell_ijk = (B.q * owned_ijk) + inner_ijk;
      int const viz_cell = gviz_cell_grid.index(viz_cell_ijk);
      real const val = eval(U, block, cell, eq, B, CELL_LOC, pt);
      var.h_view(viz_cell, 0) = val;
    });
  };
  for_each("vtk::get_variable", owned_cells, functor);
  return var;
}

void write_out(
    std::stringstream& stream,
    Input const* input,
    State const* state,
    int const soln_idx,
    int const block)
{
  int const neqs = input->passive.value().names.size();
  for (int var = 0; var < neqs; ++var) {
    std::string const name = input->passive.value().names[var];
    dgt::vtk::write_vtr_field(stream, name, get_variable(input, state, block, soln_idx, var));
  }
}

}
}
