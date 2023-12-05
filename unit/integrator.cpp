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

    real pos_resid;
    real vel_resid;

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

    void compute_explicit_residual(
        int const from,
        int const,
        int const,
        real const,
        real const) override
    {
      pos_resid = vel[from];
      vel_resid = -(k/m)*pos[from];
    }

    void advance_explicitly(
        int const from,
        int const,
        int const into,
        real const dt,
        real const) override
    {
      pos[into] = pos[from] + dt*pos_resid;
      vel[into] = vel[from] + dt*vel_resid;
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

static real test_ssprk_spring(int const order, real const dt)
{
  real const k = 1.;
  real const m = 1.;
  real const w = std::sqrt(k/m);
  real const pi = 3.14159265358979323;
  auto spring = std::make_shared<Spring>(k, m, 1., 0.);
  Physics P = {spring};
  SSPRK I(order);
  real t = 0.;
  while (t <= 4.*pi/w) {
    for (int stage = 0; stage < I.num_stages(); ++stage) {
      I.do_stage(P, stage, dt, t);
    }
    t += dt;
  }
  real const pos_exact = std::cos(w*t);
  real const vel_exact = -w*std::sin(w*t);
  real const pos_error = pos_exact - spring->pos[0];
  real const vel_error = vel_exact - spring->vel[0];
  real const error = std::sqrt(pos_error*pos_error + vel_error*vel_error);
  return error;
}

static void check_rate(
    std::vector<real> const& errors,
    real const expected,
    real const tolerance)
{
  for (std::size_t i = 0; i < errors.size()-1; ++i) {
    real const reduction = errors[i] / errors[i+1];
    real const rate = std::log(reduction) / std::log(2.);
    EXPECT_NEAR(rate, expected, tolerance);
  }
}

static void check_ssprk_spring(
    int const order,
    real const dt,
    real const tolerance)
{
  int const num_iters = 5;
  std::vector<real> errors(num_iters, 0.);
  for (int i = 1; i <= num_iters; ++i) {
    real const factor = 1 << i;
    errors[i-1] = test_ssprk_spring(order, dt/factor);
  }
  check_rate(errors, order, tolerance);
}

TEST(integrator, ssprk1_spring)
{
  check_ssprk_spring(1, 0.01, 0.02);
}

TEST(integrator, ssprk2_spring)
{
  check_ssprk_spring(2, 0.01, 0.01);
}

TEST(integrator, ssprk3_spring)
{
  check_ssprk_spring(3, 0.1, 0.01);
}
