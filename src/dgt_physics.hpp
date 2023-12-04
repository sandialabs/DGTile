#pragma once

#include <vector>

#include "dgt_defines.hpp"

namespace dgt {

class PhysicsPackage
{

  public:

    PhysicsPackage() = default;
    virtual ~PhysicsPackage() {}
    virtual std::string name() const = 0;

  public:

    virtual void apply_initial_conditions() {}

  public:

    virtual double compute_time_step() = 0;
    virtual void handle_visualization() {}
    virtual void handle_steps() {}
    virtual void handle_history() {}

  public:

    virtual void begin_explicit_stage(
        int const from1,
        int const from2,
        int const into,
        real const dt)
    {
      (void)from1;
      (void)from2;
      (void)into;
      (void)dt;
    }

    virtual void compute_explicit_residual(
        int const from1,
        int const from2,
        int const into,
        real const dt)
    {
      (void)from1;
      (void)from2;
      (void)into;
      (void)dt;
    }

    virtual void advance_explicitly(
        int const from1,
        int const from2,
        int const into,
        real const dt)
    {
      (void)from1;
      (void)from2;
      (void)into;
      (void)dt;
    }

    virtual void end_explicit_stage(
        int const from1,
        int const from2,
        int const into,
        real const dt)
    {
      (void)from1;
      (void)from2;
      (void)into;
      (void)dt;
    }
  
  public:

    virtual void begin_implicit_stage(
        int const from,
        int const into,
        real const dt)
    {
      (void)from;
      (void)into;
      (void)dt;
    }

    virtual void advance_implicitly(
        int const from,
        int const into,
        real const dt)
    {
      (void)from;
      (void)into;
      (void)dt;
    }

    virtual void end_explicit_stage(
        int const from,
        int const into,
        real const dt)
    {
      (void)from;
      (void)into;
      (void)dt;
    }

  public:

    virtual void axpby(
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

};

using Physics = std::vector<PhysicsPackage*>;

}
