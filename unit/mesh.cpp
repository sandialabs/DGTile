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

TEST(mesh, init_1D)
{
  mpicpp::comm comm = mpicpp::comm::world();
  Mesh mesh = get_example_mesh(1, &comm);
  EXPECT_EQ(mesh.dim(), 1);
  EXPECT_EQ(mesh.num_total_blocks(), 2);
  EXPECT_EQ(mesh.num_owned_blocks(), 2);
  EXPECT_EQ(mesh.num_total_cells(), 4);
  EXPECT_EQ(mesh.num_owned_cells(), 4);
}

TEST(mesh, init_2D)
{
  mpicpp::comm comm = mpicpp::comm::world();
  Mesh mesh = get_example_mesh(2, &comm);
  EXPECT_EQ(mesh.dim(), 2);
  EXPECT_EQ(mesh.num_total_blocks(), 4);
  EXPECT_EQ(mesh.num_owned_blocks(), 4);
  EXPECT_EQ(mesh.num_total_cells(), 16);
  EXPECT_EQ(mesh.num_owned_cells(), 16);
}

TEST(mesh, init_3D)
{
  mpicpp::comm comm = mpicpp::comm::world();
  Mesh mesh = get_example_mesh(3, &comm);
  EXPECT_EQ(mesh.dim(), 3);
  EXPECT_EQ(mesh.num_total_blocks(), 8);
  EXPECT_EQ(mesh.num_owned_blocks(), 8);
  EXPECT_EQ(mesh.num_total_cells(), 64);
  EXPECT_EQ(mesh.num_owned_cells(), 64);
}

TEST(mesh, add_modal_field)
{
  mpicpp::comm comm = mpicpp::comm::world();
  Mesh mesh = get_example_mesh(2, &comm);
  int const num_blocks = 4;
  int const num_ghost_cells_per_block = 16;
  int const num_ghost_faces_per_block = 20;
  int const num_modes = 4;
  int const num_face_pts = 2;
  int const num_eqs = 5;
  Field<real***> U0 = mesh.get_solution("hydro", 0);
  Field<real***> U1 = mesh.get_solution("hydro", 1);
  Field<real***> R = mesh.get_residual("hydro");
  Field<real***> Fx = mesh.get_fluxes("hydro", X);
  Field<real***> Fy = mesh.get_fluxes("hydro", Y);
  EXPECT_EQ(U0.name(), "hydro_0");
  EXPECT_EQ(U1.name(), "hydro_1");
  EXPECT_EQ(R.name(), "hydro_residual");
  EXPECT_EQ(Fx.name(), "hydro_fluxes_X");
  EXPECT_EQ(Fy.name(), "hydro_fluxes_Y");
  EXPECT_EQ(U0.get().size(), num_blocks);
  EXPECT_EQ(U1.get().size(), num_blocks);
  EXPECT_EQ(R.get().size(), num_blocks);
  EXPECT_EQ(Fx.get().size(), num_blocks);
  EXPECT_EQ(Fy.get().size(), num_blocks);
  for (int block = 0; block < num_blocks; ++block) {
    EXPECT_EQ(U0.get()[block].extent(0), num_ghost_cells_per_block);
    EXPECT_EQ(U1.get()[block].extent(0), num_ghost_cells_per_block);
    EXPECT_EQ(R.get()[block].extent(0), num_ghost_cells_per_block);
    EXPECT_EQ(U0.get()[block].extent(1), num_eqs);
    EXPECT_EQ(U1.get()[block].extent(1), num_eqs);
    EXPECT_EQ(R.get()[block].extent(1), num_eqs);
    EXPECT_EQ(U0.get()[block].extent(2), num_modes);
    EXPECT_EQ(U1.get()[block].extent(2), num_modes);
    EXPECT_EQ(R.get()[block].extent(2), num_modes);
    EXPECT_EQ(Fx.get()[block].extent(0), num_ghost_faces_per_block);
    EXPECT_EQ(Fy.get()[block].extent(0), num_ghost_faces_per_block);
    EXPECT_EQ(Fx.get()[block].extent(1), num_face_pts);
    EXPECT_EQ(Fy.get()[block].extent(1), num_face_pts);
  }
}