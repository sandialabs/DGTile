#pragma once

#include <dgt_physics.hpp>

namespace app {

struct Input;
struct State;

enum {DENS = 0, MMTM = 1, ENER = 4, NEQ = 5};

static constexpr real rho_floor = 1.e-12;
static constexpr real eint_floor = 1.e-12;

template <class FieldT>
DGT_METHOD inline Vec<real, NEQ> get_avg(
    FieldT const& U,
    int const block,
    int const cell)
{
  Vec<real, NEQ> U_avg;
  for (int eq = 0; eq < NEQ; ++eq) {
    U_avg[eq] = U[block](cell, eq, 0);
  }
  return U_avg;
}

DGT_METHOD inline Vec3<real> get_vec3(
    Vec<real, NEQ> const& U,
    int const eq)
{
  return Vec3<real>(U[eq + X], U[eq + Y], U[eq + Z]);
}

DGT_METHOD inline real get_eint(Vec<real, NEQ> const& U)
{
  real const E = U[ENER];
  real const rho = U[DENS];
  Vec3<real> const v = get_vec3(U, MMTM) / rho;
  real const half_v2 = 0.5 * dot(v,v);
  return E/rho - half_v2;
}

DGT_METHOD inline real e_from_rho_p(
    real const rho,
    real const p,
    real const gamma)
{
  return p/(rho*(gamma-1.));
}

DGT_METHOD inline real p_from_rho_e(
    real const rho,
    real const e,
    real const gamma)
{
  return rho*e*(gamma-1.);
}

DGT_METHOD inline real c_from_rho_p(
    real const rho,
    real const p,
    real const gamma)
{
  return std::sqrt(gamma*p/rho);
}

DGT_METHOD inline Vec<real, NEQ> get_physical_flux(
    Vec<real, NEQ> const& U,
    real const p,
    int const j)
{
  Vec<real, NEQ> F;
  real const rho = U[DENS];
  real const En = U[ENER];
  Vec3<real> const v = get_vec3(U, MMTM) / rho;
  F[DENS]     = v[j] * rho;
  F[MMTM + X] = v[j] * v[X] * rho;
  F[MMTM + Y] = v[j] * v[Y] * rho;
  F[MMTM + Z] = v[j] * v[Z] * rho;
  F[ENER]     = v[j] * (En + p);
  F[MMTM + j] += p;
  return F;
}

class Hydro : public PhysicsPackage
{
  private:
    Input const* m_input;
    State* m_state;
  public:
    Hydro(Input const* in, State* state);
    std::string name() const override { return "Hydro"; }
    void apply_initial_conditions() override;
    double compute_time_step() override;
    void handle_vtk(std::stringstream&, int) override;
    void handle_history() override;
    void begin_explicit_stage(int, int, int, real, real) override;
    void compute_explicit_residual(int, int, int, real, real) override;
    void advance_explicitly(int, int, int, real, real) override;
    void axpby(int, real, int, real, int) override;
};

namespace hydro {
void apply_ics(Input const* in, State* state);
real compute_time_step(Input const* in, State* state);
void write_out(
    std::stringstream& stream, Input const* in, State const* state,
    int const soln_idx, int const block);
void compute_fluxes(Input const* in, State* state, int const from);
void compute_volume_integral(Input const* in, State* state, int const from);
void compute_face_integral(State* state);
void advance_explicitly(State* state, int const from, int const into, real const dt);
void preserve_bounds(State* state, int const into);
}

}
