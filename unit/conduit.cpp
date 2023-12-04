#include <dgt_conduit.hpp>
#include <dgt_mesh.hpp>

#include <gtest/gtest.h>

using namespace dgt;

static Mesh get_example_mesh(int const dim, mpicpp::comm* comm)
{
  Grid3 cell_grid(0,0,0);
  Grid3 block_grid(0,0,0);
  Box3<real> domain({0,0,0},{0,0,0});
  for (int axis = 0; axis < dim; ++axis) {
    domain.upper()[axis] = 1.;
    cell_grid.extents()[axis] = 2;
    block_grid.extents()[axis] = 2;
  }
  Mesh mesh;
  mesh.set_comm(comm);
  mesh.set_domain(domain);
  mesh.set_cell_grid(cell_grid);
  mesh.set_basis(1, 2, true);
  mesh.add_modal({"hydro", 2, 5, true});
  mesh.initialize(block_grid);
  return mesh;
}

TEST(conduit, example)
{
  mpicpp::comm comm = mpicpp::comm::world();
  Mesh mesh = get_example_mesh(2, &comm);
  dgt::silo::dummy(mesh);
}
