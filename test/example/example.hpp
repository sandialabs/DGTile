#include <string>

#include <dgt_dg.hpp>
#include <dgt_mesh.hpp>
#include <dgt_view.hpp>

#include <mpicpp.hpp>

namespace example {

using namespace dgt;

namespace inputs {

struct Basis
{
  int polynomial_order = -1;
  int quadrature_rule = -1;
  bool tensor_product = true;
};

struct Time
{
  real cfl = 0.;
  real end_time = 0.;
};

struct Mesh
{
  Grid3 cell_grid;
  Box3<real> domain;
  Vec3<bool> periodic;
};

}

struct Input
{
  mpicpp::comm comm;
  std::string name = "";
  std::string input_file_name = "";
  inputs::Basis basis;
  inputs::Time time;
  inputs::Mesh mesh;
};

void run_lua_file(std::string const& path);
void run(Input const& in);

}
