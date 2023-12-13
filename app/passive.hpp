#pragma once

#include <limits>

#include <dgt_physics.hpp>

namespace app {

class Passive : public PhysicsPackage
{
  private:
    Input const* m_input;
    State* m_state;
    int m_num_vars;
  public:
    Passive(Input const* in, State* state);
    std::string name() const override { return "Passive Scalars"; }
    void apply_initial_conditions() override;
    real compute_time_step() override { return std::numeric_limits<double>::max(); }
    void handle_vtk(std::stringstream&, int) override;
//    void begin_explicit_stage(int, int, int, real, real) override;
};

namespace passive {
void apply_ics(Input const* input, State* state);
void write_out(
    std::stringstream& stream, Input const* in, State const* state,
    int const soln_idx, int const block);
void compute_fluxes(Input const* in, State* state, int const from);
void compute_volume_integral(Input const* in, State* state, int const from);
}

}
