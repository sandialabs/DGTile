#include <dgt_mesh.hpp>

#include <gtest/gtest.h>

using namespace dgt;

static Mesh get_example_mesh(int const dim, mpicpp::comm* comm)
{
  Basis<View> basis = build_basis<View>(dim, 1, 2, true);
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
  mesh.set_basis(basis);
  mesh.init(block_grid);
  return mesh;
}

TEST(mesh, init_1D)
{
  mpicpp::comm comm = mpicpp::comm::world();
  Mesh mesh = get_example_mesh(1, &comm);
  EXPECT_EQ(mesh.num_total_blocks(), 2);
  EXPECT_EQ(mesh.num_owned_blocks(), 2);
  EXPECT_EQ(mesh.num_total_cells(), 4);
}

TEST(mesh, init_2D)
{
  mpicpp::comm comm = mpicpp::comm::world();
  Mesh mesh = get_example_mesh(2, &comm);
  EXPECT_EQ(mesh.num_total_blocks(), 4);
  EXPECT_EQ(mesh.num_owned_blocks(), 4);
  EXPECT_EQ(mesh.num_total_cells(), 16);
}

TEST(mesh, init_3D)
{
  mpicpp::comm comm = mpicpp::comm::world();
  Mesh mesh = get_example_mesh(3, &comm);
  EXPECT_EQ(mesh.num_total_blocks(), 8);
  EXPECT_EQ(mesh.num_owned_blocks(), 8);
  EXPECT_EQ(mesh.num_total_cells(), 64);
}
