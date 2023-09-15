#include <string>

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
  std::vector<function_ptr<real>> densities;
  std::vector<function_ptr<real>> pressures;
  function_ptr<Vec3<real>> velocity;
};

struct Materials
{
  std::vector<real> gammas;
};

}

struct Input
{
  std::string name = "";
  std::string input_file_name = "";
  int num_materials = -1;
  inputs::Basis basis;
  inputs::Time time;
  inputs::Mesh mesh;
  inputs::Materials materials;
  inputs::InitialConditions ics;
};

struct Equations
{
  enum {RHO, MMTM, ENER, NVAR};
  int offsets[NVAR];
  Equations() = default;
  Equations(int const num_mats);
  int rho(int const mat) { return offsets[RHO] + mat; }
  int mmtm() { return offsets[MMTM]; }
  int ener() { return offsets[ENER]; }
  int num_eqs() { return offsets[ENER] + 1; }
};

struct State
{
  Equations eqs;
  Mesh mesh;
};

void run_lua_file(std::string const& path);
void run(mpicpp::comm* comm, Input const& in);
void setup(State& state, mpicpp::comm* comm, Input const& in);

}
