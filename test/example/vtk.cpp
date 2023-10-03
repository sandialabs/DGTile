#include <filesystem>

#include <dgt_cartesian.hpp>
#include <dgt_for_each.hpp>
#include <dgt_print.hpp>
#include <dgt_vtk.hpp>

#include "example.hpp"

namespace example {

template <class T>
DGT_HOST_DEVICE void assign_variable(
    dgt::vtk::VtkView<real> v, int const idx, T const& val);

template <>
DGT_HOST_DEVICE void assign_variable(
    dgt::vtk::VtkView<real> v, int const idx, real const& val)
{
  v.d_view(idx, 0) = val;
}

template <>
DGT_HOST_DEVICE void assign_variable(
    dgt::vtk::VtkView<real> v, int const idx, Vec3<real> const& val)
{
  v.d_view(idx, X) = val.x();
  v.d_view(idx, Y) = val.y();
  v.d_view(idx, Z) = val.z();
}

struct Data
{
  typename Field<real***>::accessor_t U;
  Basis<View> B;
  Array<EoS, nmax_mat> eos;
  Equations eqs;
  int block;
};

template <class Function>
dgt::vtk::VtkView<real> get_variable(
    Function const& function,
    State const& state,
    int const block,
    int const num_comps,
    int const soln_idx,
    int const mat = 0)
{
  Mesh const& mesh = state.mesh;
  Basis<View> const& B = mesh.basis();
  Grid3 const cell_grid = mesh.cell_grid();
  Grid3 const inner_grid = tensor_bounds(B.dim, B.q-1);
  Grid3 const ginner_grid = generalize(B.dim, inner_grid);
  Grid3 const viz_cell_grid = vtk::get_viz_cell_grid(cell_grid, B.q);
  Grid3 const gviz_cell_grid = generalize(B.dim, viz_cell_grid);
  Subgrid3 const owned_cells = get_owned_cells(cell_grid);
  Vec3<int> const ghost_offset = owned_cells.lower();
  auto const U = mesh.get_solution("hydro", soln_idx).get();
  dgt::vtk::VtkView<real> var;
  Kokkos::resize(var, gviz_cell_grid.size(), num_comps);
  Data data = {U, B, state.eos, state.eqs, block};
  auto functor = [=] DGT_HOST_DEVICE (Vec3<int> const& cell_ijk) {
    int const cell = cell_grid.index(cell_ijk);
    Vec3<int> const owned_ijk = cell_ijk - ghost_offset;
    inner_for_each(ginner_grid,
    [&] (Vec3<int> const& inner_ijk) DGT_ALWAYS_INLINE {
      int const pt = ginner_grid.index(inner_ijk);
      Vec3<int> const viz_cell_ijk = (B.q * owned_ijk) + inner_ijk;
      int const viz_cell = gviz_cell_grid.index(viz_cell_ijk);
      auto const val = function(data, cell, pt, mat);
      assign_variable(var, viz_cell, val);
    });
  };
  for_each("vtk::get_variable", owned_cells, functor);
  return var;
}

static constexpr int CELL = basis_locations::CELL;

struct density
{
  DGT_METHOD inline real operator()(
      Data const& d,
      int const cell,
      int const pt,
      int const mat) const
  {
    return eval(d.U, d.block, cell, d.eqs.rho(mat), d.B, CELL, pt);
  }
};

static void write_mesh(
    std::filesystem::path const& path,
    State const& state,
    int const soln_idx)
{
  Mesh const& mesh = state.mesh;
  int const nblocks = mesh.num_owned_blocks();
  density f_rho;
  for (int block = 0; block < nblocks; ++block) {
    std::stringstream stream;
    std::filesystem::path const block_path = path / fmt::format("{}.vtr", block);
    dgt::vtk::write_vtr_start(stream, block, mesh, state.time, state.step);

    // TODO debug this
    dgt::vtk::write_vtr_field(stream, "rho", get_variable(f_rho, state, block, 1, soln_idx));

    dgt::vtk::write_vtr_end(stream);
    dgt::write_stream(block_path, stream);
  }
  if (state.mesh.comm()->rank() == 0) {
    std::filesystem::path const vtm_path = path / "blocks.vtm";
    std::stringstream stream;
    dgt::vtk::write_vtm(stream, "", mesh.num_total_blocks());
    dgt::write_stream(vtm_path, stream);
  }


  (void)soln_idx;
}

void write_out(Input const& in, State const& state, int soln_idx)
{
  static int ctr = 0;
  std::filesystem::path const out_dir = in.name;
  std::filesystem::path const vtk_dir = out_dir / "vtk";
  std::filesystem::path const path = vtk_dir / std::to_string(ctr);
  std::filesystem::create_directory(vtk_dir);
  std::filesystem::create_directory(path);
  write_mesh(path, state, soln_idx);
  ctr++;
}

}
