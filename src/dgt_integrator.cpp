#include <fmt/core.h>

#include "dgt_integrator.hpp"

namespace dgt {

static void verify_stage(
    std::string const& name,
    int const max_stages,
    int const stage)
{
  if ((stage < 0) || (stage >= max_stages)) {
    std::string const msg = fmt::format(
        "dgt::{} -> invalid stage `{}`", name, stage);
    throw std::runtime_error(msg);
  }
}

SSPRK::SSPRK(int const order)
{
  m_order = order;
  if ((order < 1) || (order > 3)) {
    std::string const msg = fmt::format(
        "dgt::SSPRK -> invalid order `{}`", order);
    throw std::runtime_error(msg);
  }
}

std::string SSPRK::name() const
{
  return fmt::format("SSPRK{}", m_order);
}

int SSPRK::required_containers() const
{
  if (m_order == 1) return 1;
  if (m_order == 2) return 2;
  if (m_order == 3) return 2;
  return -1;
}

int SSPRK::num_stages() const
{
  if (m_order == 1) return 1;
  if (m_order == 2) return 2;
  if (m_order == 3) return 3;
  return -1;
}

static constexpr int ssprk_table[3][4] = {
  {0, 0, 0, 0},
  {0, 1, 1, 0},
  {0, 1, 1, 1}
};

struct AXPBY
{
  int r;
  real a;
  int x;
  real b;
  int y;
};

static AXPBY get_ssprk_axpby(int const order, int const stage)
{
  if ((order == 2) && (stage == 1)) return {0, 0.5, 0, 0.5, 1};
  if ((order == 3) && (stage == 1)) return {1, 0.75, 0, 0.25, 1};
  if ((order == 3) && (stage == 2)) return {0, 1./3., 0, 2./3., 1};
  return {-1, 0., -1, 0., -1};
}

static void begin_explicit_stage(
    Physics& physics,
    int const from1,
    int const from2,
    int const into,
    real const dt,
    real const t)
{
  for (auto p : physics) {
    p->begin_explicit_stage(from1, from2, into, t, dt);
  }
}

static void compute_explicit_residual(
    Physics& physics,
    int const from1,
    int const from2,
    int const into,
    real const dt,
    real const t)
{
  for (auto p : physics) {
    p->compute_explicit_residual(from1, from2, into, t, dt);
  }
}

static void advance_explicitly(
    Physics& physics,
    int const from1,
    int const from2,
    int const into,
    real const dt,
    real const t)
{
  for (auto p : physics) {
    p->advance_explicitly(from1, from2, into, dt, t);
  }
}

static void end_explicit_stage(
    Physics& physics,
    int const from1,
    int const from2,
    int const into,
    real const dt,
    real const t)
{
  for (auto p : physics) {
    p->end_explicit_stage(from1, from2, into, dt, t);
  }
}

static void axpby(
    Physics& physics,
    int const r,
    real const a,
    int const x,
    real const b,
    int const y)
{
  for (auto p : physics) {
    p->axpby(r, a, x, b, y);
  }
}

void SSPRK::do_stage(
    Physics& physics,
    int const stage,
    real const t,
    real const dt)
{
  verify_stage(name(), num_stages(), stage);
  int const index = m_order-1;
  int const from = ssprk_table[index][stage];
  int const into = ssprk_table[index][stage+1];
  AXPBY const I = get_ssprk_axpby(m_order, stage);
  dgt::begin_explicit_stage(physics, from, -1, into, t, dt);
  dgt::compute_explicit_residual(physics, from, -1, into, t, dt);
  dgt::advance_explicitly(physics, from, -1, into, t, dt);
  dgt::end_explicit_stage(physics, from, -1, into, t, dt);
  if (stage == 0) return;
  dgt::axpby(physics, I.r, I.a, I.x, I.b, I.y);
}

IntegratorPtr create_integrator(std::string const& name)
{
  if (name == "ssprk1") {
    return std::make_unique<SSPRK>(1);
  } else if (name == "ssprk2") {
    return std::make_unique<SSPRK>(2);
  } else if (name == "ssprk3") {
    return std::make_unique<SSPRK>(3);
  } else {
    throw std::runtime_error("dgt::Integrator - invalid name: " + name);
  }
}

}
