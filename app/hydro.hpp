#pragma once

#include <dgt_physics.hpp>

#include "app.hpp"

namespace app {

class Hydro : public PhysicsPackage
{
  public:
    Euler();
    std::string name() const { return "Hydro"; }
    void apply_initial_conditions();
};

}
