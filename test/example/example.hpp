#include <string>

#include <dgt_array.hpp>
#include <dgt_mesh.hpp>
#include <dgt_when.hpp>

#include <mpicpp.hpp>

namespace example {

using namespace dgt;

namespace inputs {

template <class T>
struct function
{
  virtual T operator()(Vec3<real> const& x) = 0;
  virtual ~function() = default;
};

template <class T>
using function_ptr = std::unique_ptr<function<T>>;

struct Time
{
  real cfl = 0.;
  real end_time = 0.;
  WhenPtr to_terminal;
};

struct Basis
{
  int polynomial_order = -1;
  int quadrature_rule = -1;
  bool tensor_product = true;
};

struct Mesh
{
  Grid3 block_grid = {0,0,0};
  Grid3 cell_grid = {0,0,0};
  Box3<real> domain = {{0,0,0}, {0,0,0}};
  Vec3<bool> periodic = {false, false, false};
};

struct InitialConditions
{
  function_ptr<real> density;
  function_ptr<real> pressure;
  function_ptr<Vec3<real>> velocity;
};

}

enum {DENS=0, MMTM=1, ENER=4, NEQ=5};

struct Input
{
  std::string name = "";
  std::string input_file_name = "";
  int gamma = -1;
  inputs::Basis basis;
  inputs::Time time;
  inputs::Mesh mesh;
  inputs::InitialConditions ics;
};

struct EoS
{
  private:
    real m_gamma = 0.;
  public:
    EoS() = default;
    EoS(real gamma) : m_gamma(gamma) {}
    DGT_METHOD real e_from_rho_p(real const rho, real const p) const
    {
      return p/(rho*(m_gamma-1.));
    }
};

struct State
{
  EoS eos;
  Mesh mesh;
  real time = 0.;
  real dt = 0.;
  int step = 0;
};

void run_lua_file(std::string const& path);
void run(mpicpp::comm* comm, Input const& in);
void setup(State& state, mpicpp::comm* comm, Input const& in);
void write_out(Input const& in, State const& state, int soln_idx);
real compute_dt(Input const& in, State const& state);

}
