#include "app.hpp"
#include "passive.hpp"

namespace app {

Passive::Passive(Input const* input_in, State* state_in)
{
  m_input = input_in;
  m_state = state_in;
  m_num_vars = m_input->passive.value().names.size();
  int const num_stored_solns = m_state->integrator->required_containers();
  m_state->mesh.add_modal({"passive", num_stored_solns, m_num_vars, true});
}

void Passive::apply_initial_conditions()
{
  passive::apply_ics(m_input, m_state);
}

void Passive::handle_vtk(std::stringstream& stream, int const block)
{
  passive::write_out(stream, m_input, m_state, 0, block);
}

}
