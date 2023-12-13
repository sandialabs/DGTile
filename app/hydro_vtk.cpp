#include <filesystem>

#include <dgt_cartesian.hpp>
#include <dgt_for_each.hpp>
#include <dgt_print.hpp>
#include <dgt_vtk.hpp>

#include "app.hpp"
#include "hydro.hpp"

namespace app {
namespace hydro {

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
  real gamma;
  int block;
};

template <class Function>
dgt::vtk::VtkView<real> get_variable(
    Function const& function,
    Input const* input,
    State const* state,
    int const block,
    int const num_comps,
    int const soln_idx)
{
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
  auto const U = mesh.get_solution("hydro", soln_idx).get();
  dgt::vtk::VtkView<real> var;
  Kokkos::resize(var, gviz_cell_grid.size(), num_comps);
  Data data = {U, B, gamma, block};
  auto functor = [=] DGT_HOST_DEVICE (Vec3<int> const& cell_ijk) {
    int const cell = cell_grid.index(cell_ijk);
    Vec3<int> const owned_ijk = cell_ijk - ghost_offset;
    inner_for_each(ginner_grid,
    [&] (Vec3<int> const& inner_ijk) DGT_ALWAYS_INLINE {
      int const pt = ginner_grid.index(inner_ijk);
      Vec3<int> const viz_cell_ijk = (B.q * owned_ijk) + inner_ijk;
      int const viz_cell = gviz_cell_grid.index(viz_cell_ijk);
      auto const val = function(data, cell, pt);
      assign_variable(var, viz_cell, val);
    });
  };
  for_each("vtk::get_variable", owned_cells, functor);
  return var;
}

static constexpr int CELL_LOC = basis_locations::CELL;

struct density
{
  DGT_METHOD inline real operator()(Data const& d, int const cell, int const pt) const
  {
    return eval(d.U, d.block, cell, DENS, d.B, CELL_LOC, pt);
  }
};

struct momentum
{
  DGT_METHOD inline Vec3<real> operator()(Data const& d, int const cell, int const pt) const
  {
    return eval_vec3(d.U, d.block, cell, MMTM, d.B, CELL_LOC, pt);
  }
};

struct total_energy
{
  DGT_METHOD inline real operator()(Data const& d, int const cell, int const pt) const
  {
    return eval(d.U, d.block, cell, ENER, d.B, CELL_LOC, pt);
  }
};

struct velocity
{
  density f_rho;
  momentum f_mmtm;
  DGT_METHOD inline Vec3<real> operator()(Data const& d, int const cell, int const pt) const
  {
    real const rho = f_rho(d, cell, pt);
    Vec3<real> const mmtm = f_mmtm(d, cell, pt);
    return mmtm / rho;
  }
};

struct internal_energy
{
  density f_rho;
  velocity f_v;
  total_energy f_En;
  DGT_METHOD inline real operator()(Data const& d, int const cell, int const pt) const
  {
    real const rho = f_rho(d, cell, pt);
    real const En = f_En(d, cell, pt);
    Vec3<real> const v = f_v(d, cell, pt);
    real const half_v2 = 0.5*dot(v,v);
    return En/rho - half_v2;
  }
};

struct pressure
{
  density f_rho;
  internal_energy f_e;
  DGT_METHOD inline real operator()(Data const& d, int const cell, int const pt) const
  {
    real const rho = f_rho(d, cell, pt);
    real const e = f_e(d, cell, pt);
    return p_from_rho_e(rho, e, d.gamma);
  }
};

void write_out(
    std::stringstream& stream,
    Input const* input,
    State const* state,
    int const soln_idx,
    int const block)
{
  density f_rho;
  momentum f_mmtm;
  total_energy f_En;
  internal_energy f_e;
  pressure f_p;
  velocity f_v;
  dgt::vtk::write_vtr_field(stream, "density", get_variable(f_rho, input, state, block, 1, soln_idx));
  dgt::vtk::write_vtr_field(stream, "momentum", get_variable(f_mmtm, input, state, block, DIMENSIONS, soln_idx));
  dgt::vtk::write_vtr_field(stream, "total_energy", get_variable(f_En, input, state, block, 1, soln_idx));
  dgt::vtk::write_vtr_field(stream, "internal_energy", get_variable(f_e, input, state, block, 1, soln_idx));
  dgt::vtk::write_vtr_field(stream, "pressure", get_variable(f_p, input, state, block, 1, soln_idx));
  dgt::vtk::write_vtr_field(stream, "velocity", get_variable(f_v, input, state, block, DIMENSIONS, soln_idx));
}

}
}
