#pragma once

#include <sstream>
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

    virtual void handle_vtk(
        std::stringstream& stream,
        int const block)
    {
      (void)stream;
      (void)block;
    }

    virtual void handle_history() {}

  public:

    virtual void begin_explicit_stage(
        int const from1,
        int const from2,
        int const into,
        real const dt,
        real const t)
    {
      (void)from1;
      (void)from2;
      (void)into;
      (void)dt;
      (void)t;
    }

    virtual void compute_explicit_residual(
        int const from1,
        int const from2,
        int const into,
        real const dt,
        real const t)
    {
      (void)from1;
      (void)from2;
      (void)into;
      (void)dt;
      (void)t;
    }

    virtual void advance_explicitly(
        int const from1,
        int const from2,
        int const into,
        real const dt,
        real const t)
    {
      (void)from1;
      (void)from2;
      (void)into;
      (void)dt;
      (void)t;
    }

    virtual void end_explicit_stage(
        int const from1,
        int const from2,
        int const into,
        real const dt,
        real const t)
    {
      (void)from1;
      (void)from2;
      (void)into;
      (void)dt;
      (void)t;
    }
  
  public:

    virtual void begin_implicit_stage(
        int const from,
        int const into,
        real const dt,
        real const t)
    {
      (void)from;
      (void)into;
      (void)dt;
      (void)t;
    }

    virtual void advance_implicitly(
        int const from,
        int const into,
        real const dt,
        real const t)
    {
      (void)from;
      (void)into;
      (void)dt;
      (void)t;
    }

    virtual void end_explicit_stage(
        int const from,
        int const into,
        real const dt,
        real const t)
    {
      (void)from;
      (void)into;
      (void)dt;
      (void)t;
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

using PhysicsPackagePtr = std::shared_ptr<PhysicsPackage>;
using Physics = std::vector<PhysicsPackagePtr>;

}
