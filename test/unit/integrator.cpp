#include <dgt_integrator.hpp>

#include <gtest/gtest.h>

#include <iomanip>
#include <iostream>

using namespace dgt;

class Spring : public PhysicsPackage
{

  public:

    real k;
    real m;

    real pos[2];
    real vel[2];

  public:

    Spring(
        real const k_in,
        real const m_in,
        real const x0,
        real const v0)
    {
      k = k_in;
      m = m_in;
      pos[0] = x0;
      vel[0] = v0;
    }

    std::string name() const override { return "Spring"; }

    real compute_time_step() override { return 0.; }

    void advance_explicitly(
        int const from,
        int const,
        int const into,
        real const dt) override
    {
      pos[into] = pos[from] + dt * vel[from];
      vel[into] = vel[from] - dt * (k/m) * pos[from];
    }

    void axpby(
        int const r,
        real const a,
        int const x,
        real const b,
        int const y) override
    {
      pos[r] = a*pos[x] + b*pos[y];
      vel[r] = a*vel[x] + b*vel[y];
    }

};

static void test_ssp_spring(int const order, real const dt)
{
  real const k = 10.;
  real const m = 1.;
  real const w = std::sqrt(k/m);
  real const pi = 3.14159265358979323;
  Spring spring(k, m, 1., 0.);
  Physics P = {&spring};
  SSPRK I(order);
  real t = 0.;
  while (t <= 4.*pi/w) {
    for (int stage = 0; stage < I.num_stages(); ++stage) {
      I.do_stage(P, stage, dt);
    }
    t += dt;
  }
  std::cout << std::scientific << std::setprecision(17);
  std::cout << std::abs(1.-spring.pos[0]) << "\n";
}

TEST(integrator, spring_ssprk1)
{
  real start = 0.001;
  test_ssp_spring(1, start);
//  test_ssp_spring(1, start/2.);
//  test_ssp_spring(1, start/4.);
//  test_ssp_spring(1, start/8.);
//  test_ssp_spring(1, start/16.);
//  test_ssp_spring(1, start/32.);
//  test_ssp_spring(1, start/64.);
}
