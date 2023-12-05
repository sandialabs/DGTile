#include "app.hpp"
#include "hydro.hpp"

namespace app {

Hydro::Hydro(State* state_in, Input const* input_in)
{
  state = state_in;
  input = input_in;
  int const num_stored_solns = state->integrator->required_containers();
  state->mesh.add_modal({"hydro", num_stored_solns, 5, true});
}

void Hydro::apply_initial_conditions()
{
}

double Hydro::compute_time_step()
{
  return 1.;
}

void Hydro::handle_visualization()
{
}

void Hydro::handle_history()
{
}

void Hydro::compute_explicit_residual(
    int const from,
    int const,
    int const into,
    real const dt,
    real const)
{
  (void)from;
  (void)into;
  (void)dt;
}

void Hydro::advance_explicitly(
    int const from,
    int const,
    int const into,
    real const dt,
    real const)
{
  (void)from;
  (void)into;
  (void)dt;
}

void Hydro::axpby(
    int const r,
    real const a,
    int const x,
    real const b,
    int const y)
{
  (void)r;
  (void)a;
  (void)x;
  (void)b;
  (void)y;
}

}
