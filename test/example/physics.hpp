#pragma once

namespace example {

using namespace dgt;

class Hydro4EQ : public PhysicsPackage
{
  public:
    void apply_initial_conditions() override;
    double compute_time_step() override;
    void handle_visualization() override;
    void handle_steps() override;
    void handle_history() override;
    void begin_explicit_stage(int, int, int, real) override;
    void compute_explicit_residual(int, int, int, real) override;
    void advance_explicitly(int, int, int, real) override;
    void end_explicit_stage(int, int, int, real) override;
    void axpby(int, real, int, real, int) override;
};

}
