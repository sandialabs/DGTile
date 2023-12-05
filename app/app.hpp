#pragma once

#include <string>

#include <mpicpp.hpp>

#include <dgt_integrator.hpp>
#include <dgt_mesh.hpp>
#include <dgt_when.hpp>

namespace app {

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
  std::string integrator;
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

}

struct Input
{
  std::string name = "";
  std::string input_file_name = "";
  inputs::Time time;
  inputs::Basis basis;
  inputs::Mesh mesh;
};

struct Timer
{
  real previous_time = 0;
  int previous_step = 0;
  int total_cells = 0;
  void update(dgt::Mesh const& mesh);
  real compute(int step);
};

struct State
{
  int step = 0;
  real time = 0.;
  real dt = 0.;
  Mesh mesh;
  Timer timer;
  IntegratorPtr integrator;
};

void run_lua_file(std::string const& path);
void run(mpicpp::comm* comm, Input const& in);
void setup(mpicpp::comm* comm, Input const& in, State& s);

}
