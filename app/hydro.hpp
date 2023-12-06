#pragma once

#include <dgt_physics.hpp>

namespace app {

struct Input;
struct State;

class Hydro : public PhysicsPackage
{
  private:
    State* m_state;
    Input const* m_input;
  public:
    Hydro(State* state, Input const* in);
    std::string name() const override { return "Hydro"; }
    void apply_initial_conditions() override;
    double compute_time_step() override;
    void handle_visualization() override;
    void handle_history() override;
    void compute_explicit_residual(int, int, int, real, real) override;
    void advance_explicitly(int, int, int, real, real) override;
    void axpby(int, real, int, real, int) override;
};

}
