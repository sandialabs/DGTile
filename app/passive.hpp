#pragma once

#include <dgt_physics.hpp>

namespace app {

class Passive : public PhysicsPackage
{
  private:
    Input const* m_input;
    State* m_state;
  public:
    Passive(Input const* in, State* state) { }
    std::string name() const override { return "Passive Scalars"; }
    real compute_time_step() override { return std::numeric_limits<double>::max(); }
};

}
