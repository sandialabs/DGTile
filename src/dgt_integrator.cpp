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

struct axpby
{
  int r;
  real a;
  int x;
  real b;
  int y;
};

static axpby get_ssprk_axpby(int const order, int const stage)
{
  if ((order == 2) && (stage == 1)) return {0, 0.5, 0, 0.5, 1};
  if ((order == 3) && (stage == 1)) return {1, 0.75, 0, 0.25, 1};
  if ((order == 3) && (stage == 2)) return {0, 1./3., 0, 2./3., 1};
  return {-1, 0., -1, 0., -1};
}

void SSPRK::do_stage(Physics& physics, int const stage, real const dt)
{
  verify_stage(name(), num_stages(), stage);
  int const index = m_order-1;
  int const from = ssprk_table[index][stage];
  int const into = ssprk_table[index][stage+1];
  axpby const I = get_ssprk_axpby(m_order, stage);
  for (auto p : physics) {
    p->begin_explicit_stage(from, -1, into, dt);
    p->compute_explicit_residual(from, -1, into, dt);
    p->advance_explicitly(from, -1, into, dt);
    p->end_explicit_stage(from, -1, into, dt);
    if (stage == 0) continue;
    else p->axpby(I.r, I.a, I.x, I.b, I.y);
  }
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
