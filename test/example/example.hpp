#include <string>

#include <dgt_box3.hpp>
#include <dgt_grid3.hpp>
#include <dgt_when.hpp>

#include <mpicpp.hpp>

namespace example {

using namespace dgt;

namespace inputs {

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
  Grid3 block_grid;
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
