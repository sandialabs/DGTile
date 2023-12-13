#include <filesystem>

#include <dgt_array.hpp>
#include <dgt_cartesian.hpp>
#include <dgt_for_each.hpp>
#include <dgt_reduce.hpp>
#include <dgt_print.hpp>
#include <dgt_vtk.hpp>

#include "app.hpp"
#include "hydro.hpp"

namespace app {

Hydro::Hydro(Input const* input_in, State* state_in)
{
  m_state = state_in;
  m_input = input_in;
  int const num_stored_solns = m_state->integrator->required_containers();
  m_state->mesh.add_modal({"hydro", num_stored_solns, NEQ, true});
}

static void zero_residual(State* state)
{
  Mesh& mesh = state->mesh;
  int const num_blocks = mesh.num_owned_blocks();
  int const num_modes = mesh.basis().num_modes;
  Grid3 const cell_grid = mesh.cell_grid();
  auto R = mesh.get_residual("hydro").get();
  auto functor = [=] DGT_DEVICE (
      int const block,
      Vec3<int> const& cell_ijk) DGT_ALWAYS_INLINE
  {
    int const cell = cell_grid.index(cell_ijk);
    for (int eq = 0; eq < NEQ; ++eq) {
      for (int mode = 0; mode < num_modes; ++mode) {
        R[block](cell, eq, mode) = 0.;
      }
    }
  };
  for_each("zero_residual", num_blocks, cell_grid, functor);
}

static void advance_explicitly(
    State* state,
    int const from_idx,
    int const to_idx,
    real const dt)
{
  Mesh& mesh = state->mesh;
  int const num_blocks = mesh.num_owned_blocks();
  Grid3 const cell_grid = mesh.cell_grid();
  Subgrid3 const owned_cells = get_owned_cells(cell_grid);
  auto const B = mesh.basis();
  auto const block_info = mesh.block_info();
  auto const R = mesh.get_residual("hydro").get();
  auto const from = mesh.get_solution("hydro", from_idx).get();
  auto to = mesh.get_solution("hydro", to_idx).get();
  auto functor = [=] DGT_DEVICE (
      int const block,
      Vec3<int> const& cell_ijk) DGT_ALWAYS_INLINE
  {
    int const cell = cell_grid.index(cell_ijk);
    real const detJ = block_info.cell_detJs[block];
    for (int mode = 0; mode < B.num_modes; ++mode) {
      real const mass = detJ * B.mass[mode];
      real const fac = dt / mass;
      for (int eq = 0; eq < NEQ; ++eq) {
        real const from_eq = from[block](cell, eq, mode);
        real const R_eq = R[block](cell, eq, mode);
        to[block](cell, eq, mode) = from_eq + fac * R_eq;
      }
    }
  };
  for_each("hydro::advance_explicitly",
      num_blocks, owned_cells, functor);
}

static void axpby(
    State* state,
    Field<real***>& r,
    real const a,
    Field<real***> const& x,
    real const b,
    Field<real***> const& y)
{
  auto R = r.get();
  auto X = x.get();
  auto Y = y.get();
  Mesh const& mesh = state->mesh;
  Grid3 const cell_grid = mesh.cell_grid();
  Subgrid3 const owned_cells = get_owned_cells(cell_grid);
  int const num_blocks = mesh.num_owned_blocks();
  int const num_modes = mesh.basis().num_modes;
  auto functor = [=] DGT_DEVICE(
      int const block,
      Vec3<int> const& cell_ijk) DGT_ALWAYS_INLINE
  {
    int const cell = cell_grid.index(cell_ijk);
    for (int eq = 0; eq < NEQ; ++eq) {
      for (int mode = 0; mode < num_modes; ++mode) {
        R[block](cell, eq, mode) =
          a * X[block](cell, eq, mode) +
          b * Y[block](cell, eq, mode);
      }
    }
  };
  for_each("axpby", num_blocks, owned_cells, functor);
}

void Hydro::apply_initial_conditions()
{
  hydro::apply_ics(m_input, m_state);
}

real Hydro::compute_time_step()
{
  return hydro::compute_time_step(m_input, m_state);
}

void Hydro::handle_visualization()
{
  hydro::write_out(m_input, m_state, 0);
}

void Hydro::handle_history()
{
}

void Hydro::begin_explicit_stage(
    int const from,
    int const,
    int const,
    real const,
    real const)
{
  zero_residual(m_state);
  m_state->mesh.ghost("hydro", from);
}

void Hydro::compute_explicit_residual(
    int const from,
    int const,
    int const,
    real const,
    real const)
{
  hydro::compute_fluxes(m_input, m_state, from);
  hydro::compute_volume_integral(m_input, m_state, from);
  hydro::compute_face_integral(m_state);
}

void Hydro::advance_explicitly(
    int const from,
    int const,
    int const into,
    real const dt,
    real const)
{
  app::advance_explicitly(m_state, from, into, dt);
  hydro::preserve_bounds(m_state, into);
}

void Hydro::axpby(
    int const r,
    real const a,
    int const x,
    real const b,
    int const y)
{
  Mesh& mesh = m_state->mesh;
  auto& R = mesh.get_solution("hydro", r);
  auto const& X = mesh.get_solution("hydro", x);
  auto const& Y = mesh.get_solution("hydro", y);
  app::axpby(m_state, R, a, X, b, Y);
}

}
