#pragma once

#include <filesystem>
#include <functional>

#include "p3a_reduce.hpp"
#include "p3a_static_vector.hpp"
#include "p3a_dynamic_array.hpp"

#include "dgt_library.hpp"
#include "dgt_mesh.hpp"
#include "dgt_spatial.hpp"

namespace hydro {

using namespace dgt;

enum {RH=0,MM=1,MX=1,MY=2,MZ=3,EN=4,NEQ=5};
enum {VE=1,VX=1,VY=2,VZ=3,PR=4};

struct State;

using Exact = std::function<void(State&, Block&, View<double***>)>;

struct Input {
  std::string name;
  mpicpp::comm* comm;
  int p = -1;
  bool tensor = true;
  vector3<double> xmin;
  vector3<double> xmax;
  vector3<int> block_grid;
  vector3<int> cell_grid;
  vector3<bool> periodic;
  std::string init_amr = "";
  std::string ics = "";
  std::string amr = "";
  double gamma = -1.;
  double tfinal = -1.;
  double CFL = -1.;
  double M = 1.e+8;
  double beta = 0.5;
  double gravity = 0.;
  int gravity_axis = Y;
  int step_frequency = -1;
  double out_frequency = -1.;
  double amr_frequency = -1.;
  Exact exact_solution = nullptr;
  double error_regression = 0.;
};

struct State {
  Input in;
  Mesh mesh;
  double dt;
  double t;
  int step;
  int ssp_rk_stages;
  device_array<std::int8_t> error_code;
  std::vector<double> out_times;
  View<double***> scratch;
};

template <class T>
[[nodiscard]] P3A_ALWAYS_INLINE P3A_HOST_DEVICE inline
vector3<T> get_vec3(static_vector<T, NEQ> U, int eq) {
  return vector3<T>(U[eq + X], U[eq + Y], U[eq + Z]);
}

template <class T>
[[nodiscard]] P3A_ALWAYS_INLINE P3A_HOST_DEVICE inline
T get_eint(static_vector<T, NEQ> U) {
  T const E = U[EN];
  T const rho = U[RH];
  vector3<T> const v = get_vec3(U, MM)/rho;
  T const half_v2 = 0.5 * dot_product(v, v);
  return E/rho - half_v2;
}

template <class T>
[[nodiscard]] P3A_ALWAYS_INLINE P3A_HOST_DEVICE inline
T get_pressure(static_vector<T, NEQ> U, double gamma) {
  T const rho = U[RH];
  T const eint = get_eint(U);
  return rho * eint * (gamma - 1.);
}

template <class T>
[[nodiscard]] P3A_ALWAYS_INLINE P3A_HOST_DEVICE inline
T get_wave_speed(static_vector<T, NEQ> U, double gamma) {
  T const rho = U[RH];
  T const P = get_pressure(U, gamma);
  return sqrt(gamma * P / rho);
}

template <class T>
[[nodiscard]] P3A_ALWAYS_INLINE P3A_HOST_DEVICE inline
static_vector<T, NEQ> get_flux(
    static_vector<T, NEQ> const& U,
    T const& P,
    int j) {
  static_vector<T, NEQ> Fj;
  T const rho = U[RH];
  T const en = U[EN];
  vector3<T> const v = get_vec3(U, MM) / rho;
  Fj[RH] = rho * v[j];
  Fj[MX] = rho * v[j] * v[X];
  Fj[MY] = rho * v[j] * v[Y];
  Fj[MZ] = rho * v[j] * v[Z];
  Fj[EN] = v[j] * (en + P);
  Fj[MM + j] += P;
  return Fj;
}

template <class T>
[[nodiscard]] P3A_ALWAYS_INLINE P3A_HOST_DEVICE
static_vector<T, NEQ> get_hllc_flux(
    static_vector<T, NEQ> U[ndirs],
    static_vector<T, NEQ> F[ndirs],
    T P[ndirs],
    T c[ndirs],
    int j) {
  T rho[ndirs], v[ndirs], cs[ndirs], qs[ndirs], s[ndirs];
  static_vector<T, NEQ> U_star[ndirs], F_star[ndirs], F_hllc;
  int const mparr = MM + permute(X, j);
  int const mperp1 = MM + permute(Y, j);
  int const mperp2 = MM + permute(Z, j);
  for (int dir = 0; dir < ndirs; ++dir) {
    rho[dir] = U[dir][RH];
    v[dir] = U[dir][mparr] / rho[dir];
  }
  cs[left]  = v[left]  - c[left];
  cs[right] = v[right] + c[right];
  qs[left]  = v[right] - c[right];
  qs[right] = v[left]  + c[left];
  s[left]  = min(cs[left],  qs[left]);
  s[right] = max(cs[right], qs[right]);
  T const s_m_num =
    (rho[right] * v[right] * (s[right] - v[right]) -
     rho[left]  * v[left]  * (s[left]  - v[left])) +
    (P[left] - P[right]);
  T const s_m_den =
    rho[right] * (s[right] - v[right]) -
    rho[left]  * (s[left]  - v[left]);
  T const s_m = condition(
      s_m_den != 0.,
      s_m_num / s_m_den,
      v[right] - v[left]);
  T const P_star =
    rho[left] * (v[left] - s[left]) * (v[left] - s_m) + P[left];
  for (int k = 0; k < ndirs; ++k) {
    U_star[k] = U[k];
    T const sm_k = 1.0/ (s[k] - s_m);
    T const sq_k = s[k] - v[k];
    U_star[k][RH] = sq_k * U[k][RH] * sm_k;
    U_star[k][mparr] = (sq_k * U[k][mparr] + P_star - P[k]) * sm_k;
    U_star[k][mperp1] = sq_k * U[k][mperp1] * sm_k;
    U_star[k][mperp2] = sq_k * U[k][mperp2] * sm_k;
    U_star[k][EN] = (sq_k * U[k][EN] - P[k] * v[k] + P_star * s_m) * sm_k;
    F_star[k] = F[k] + s[k] * (U_star[k] - U[k]);
  }
  for (int eq = 0; eq < NEQ; ++eq) {
    F_hllc[eq] =
      condition(s[left] > 0.,                     F[left][eq],
      condition((s[left] <= 0.) && (s_m > 0.),    F_star[left][eq],
      condition((s_m <= 0.) && (s[right] >= 0.),  F_star[right][eq],
      condition(s[right] < 0.,                    F[right][eq],
      T(0.)))));
  }
  return F_hllc;
}

void parse_input(Input& in, mpicpp::comm* comm, std::string const& name);
void print_input(Input const& in);
void verify_input(Input const& in);

void do_initial_amr(State& state);
void do_amr(State& state);

void set_ics(State& state);
void set_exact(State& state);

double compute_stable_time_step(State& state, Block const& block);
void compute_intr_fluxes(State& state, Block& block, int axis, int soln_idx);
void compute_border_fluxes(State& state, Block& block, int axis, int dir);
void compute_amr_border_fluxes(State& state, Block& block, int axis, int dir);
void compute_vol_integral(State& state, Block& block, int soln_idx);
void compute_side_integral(Block& block, int axis);
void compute_amr_side_integral(Block& block, int axis, int dir);
void compute_gravity_source(Block& block, int soln_idx, double g, int axis);
void advance_explicitly(Block& block, int from_idx, int to_idx, double dt);
void limit(State& state, Block& block, int soln_idx, View<double***> scratch);
void preserve_bounds(State& state, Block& block, int soln_idx);
void preserve_bounds_amr(State& state, Block& block, int axis, int dir, int soln_idx);
void reflect_boundary(Border& border);

double compute_tally(Block& block, int eq);
double compute_L1_error(Block& block, View<double***> U_ex, int eq);
double compute_L2_error(Block& block, View<double***> U_ex, int eq);

void write_mesh(std::filesystem::path const& path, State const& state, int soln_idx);
void write_out(State& state, int soln_idx = 0);
void write_pvd(State& state);
void write_tree_pvd(State& state);

}
