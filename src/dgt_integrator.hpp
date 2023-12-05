#pragma once

#include <memory>

#include "dgt_physics.hpp"

namespace dgt {

class Integrator
{
  public:
    Integrator() = default;
    virtual ~Integrator() {}
    virtual std::string name() const = 0;
    virtual int required_containers() const = 0;
    virtual int num_stages() const = 0;
  public:
    virtual void do_stage(
        Physics& p,
        int const stage,
        real const dt,
        real const t)
    {
      (void)p;
      (void)stage;
      (void)dt;
      (void)t;
    }
};

class SSPRK : public Integrator
{
  private:
    int m_order;
  public:
    SSPRK(int const order);
    std::string name() const override;
    int required_containers() const override;
    int num_stages() const override;
    void do_stage(
        Physics& p,
        int const stage,
        real const dt,
        real const t) override;
};

using IntegratorPtr = std::unique_ptr<Integrator>;

IntegratorPtr create_integrator(std::string const& name);

}
