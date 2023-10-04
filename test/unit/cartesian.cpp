#include <dgt_cartesian.hpp>

#include <gtest/gtest.h>

using namespace dgt;

TEST(cartesian, grid_definitions)
{
  EXPECT_EQ(child_grid, Grid3(2,2,2));
  EXPECT_EQ(offset_grid, Subgrid3({-1,-1,-1},{2,2,2}));
}

TEST(cartesian, get_dir_sign)
{
  EXPECT_EQ(get_dir_sign(LEFT), -1);
  EXPECT_EQ(get_dir_sign(RIGHT), 1);
}

TEST(cartesian, invert_dir)
{
  EXPECT_EQ(invert_dir(LEFT), RIGHT);
  EXPECT_EQ(invert_dir(RIGHT), LEFT);
}

TEST(cartesian, get_cells_adj_face)
{
  EXPECT_EQ(get_cells_adj_face({2,2,2}, X, LEFT), Vec3<int>(2,2,2));
  EXPECT_EQ(get_cells_adj_face({2,2,2}, X, RIGHT), Vec3<int>(3,2,2));
  EXPECT_EQ(get_cells_adj_face({2,2,2}, Y, LEFT), Vec3<int>(2,2,2));
  EXPECT_EQ(get_cells_adj_face({2,2,2}, Y, RIGHT), Vec3<int>(2,3,2));
  EXPECT_EQ(get_cells_adj_face({2,2,2}, Z, LEFT), Vec3<int>(2,2,2));
  EXPECT_EQ(get_cells_adj_face({2,2,2}, Z, RIGHT), Vec3<int>(2,2,3));
}

TEST(cartesian, permute)
{
  EXPECT_EQ(permute(X, X), X);
  EXPECT_EQ(permute(X, Y), Y);
  EXPECT_EQ(permute(X, Z), Z);
  EXPECT_EQ(permute(Y, X), Y);
  EXPECT_EQ(permute(Y, Y), Z);
  EXPECT_EQ(permute(Y, Z), X);
  EXPECT_EQ(permute(Z, X), Z);
  EXPECT_EQ(permute(Z, Y), X);
  EXPECT_EQ(permute(Z, Z), Y);
}

TEST(cartesian, get_cell_center)
{
  EXPECT_EQ(get_cell_center({0,0,0}, {0.,0.,0.}, {1.,0.,0.}), Vec3<real>(0.5, 0., 0.));
  EXPECT_EQ(get_cell_center({0,0,0}, {0.,0.,0.}, {1.,1.,0.}), Vec3<real>(0.5, 0.5, 0.));
  EXPECT_EQ(get_cell_center({0,0,0}, {0.,0.,0.}, {1.,1.,1.}), Vec3<real>(0.5, 0.5, 0.5));
  EXPECT_EQ(get_cell_center({1,1,1}, {0.,0.,0.}, {1.,1.,1.}), Vec3<real>(1.5, 1.5, 1.5));
}

TEST(cartesian, map_to_physical)
{
  EXPECT_EQ(map_to_physical({0,0,0}, {0.,0.,0.}, {1.,0.,0.}, {0.,0.,0.}), Vec3<real>(0.5,0.,0.));
  EXPECT_EQ(map_to_physical({0,0,0}, {0.,0.,0.}, {1.,1.,0.}, {0.,0.,0.}), Vec3<real>(0.5,0.5,0.));
  EXPECT_EQ(map_to_physical({0,0,0}, {0.,0.,0.}, {1.,1.,1.}, {0.,0.,0.}), Vec3<real>(0.5,0.5,0.5));
}

TEST(cartesian, get_axis_name)
{
  EXPECT_EQ(get_axis_name(X), "X");
  EXPECT_EQ(get_axis_name(Y), "Y");
  EXPECT_EQ(get_axis_name(Z), "Z");
}

TEST(cartesian, get_cells_ghost_1D)
{
  Grid3 const g(8,0,0);
  EXPECT_EQ(get_cells(GHOST, g, {-1,0,0}), Subgrid3({0,0,0}, {1,0,0}));
  EXPECT_EQ(get_cells(GHOST, g, { 1,0,0}), Subgrid3({7,0,0}, {8,0,0}));
}

TEST(cartesian, get_num_cells)
{
  EXPECT_EQ(get_num_cells({5,5,5}), 125);
  EXPECT_EQ(get_num_cells({5,5,0}), 25);
  EXPECT_EQ(get_num_cells({5,0,0}), 5);
  EXPECT_THROW((void)get_num_cells({5,0,5}), std::runtime_error);
  EXPECT_THROW((void)get_num_cells({0,5,5}), std::runtime_error);
}

TEST(cartesian, get_face_grid)
{
  EXPECT_EQ(get_face_grid({5,5,5}, X), Grid3(6,5,5));
  EXPECT_EQ(get_face_grid({5,5,5}, Y), Grid3(5,6,5));
  EXPECT_EQ(get_face_grid({5,5,5}, Z), Grid3(5,5,6));
  EXPECT_EQ(get_face_grid({5,5,0}, X), Grid3(6,5,0));
  EXPECT_EQ(get_face_grid({5,5,0}, Y), Grid3(5,6,0));
  EXPECT_EQ(get_face_grid({5,0,0}, X), Grid3(6,0,0));
}

TEST(cartesian, get_num_faces)
{
  EXPECT_EQ(get_num_faces({5,5,5}, X), 150);
  EXPECT_EQ(get_num_faces({5,5,5}, Y), 150);
  EXPECT_EQ(get_num_faces({5,5,5}, Z), 150);
  EXPECT_EQ(get_num_faces({5,5,0}, X), 30);
  EXPECT_EQ(get_num_faces({5,5,0}, Y), 30);
  EXPECT_EQ(get_num_faces({5,0,0}, X), 6);
  EXPECT_THROW((void)get_num_faces({5,0,5}, X), std::runtime_error);
}

TEST(cartesian, get_owned_cells)
{
  EXPECT_EQ(get_owned_cells({6,0,0}), Subgrid3({1,0,0}, {5,0,0}));
  EXPECT_EQ(get_owned_cells({6,6,0}), Subgrid3({1,1,0}, {5,5,0}));
  EXPECT_EQ(get_owned_cells({6,6,6}), Subgrid3({1,1,1}, {5,5,5}));
}

TEST(cartesian, get_owned_faces)
{
  EXPECT_EQ(get_owned_faces({6,0,0}, X), Subgrid3({1,0,0}, {6,0,0}));
  EXPECT_EQ(get_owned_faces({6,6,0}, X), Subgrid3({1,1,0}, {6,5,0}));
  EXPECT_EQ(get_owned_faces({6,6,0}, Y), Subgrid3({1,1,0}, {5,6,0}));
  EXPECT_EQ(get_owned_faces({6,6,6}, X), Subgrid3({1,1,1}, {6,5,5}));
  EXPECT_EQ(get_owned_faces({6,6,6}, Y), Subgrid3({1,1,1}, {5,6,5}));
  EXPECT_EQ(get_owned_faces({6,6,6}, Z), Subgrid3({1,1,1}, {5,5,6}));
}

TEST(cartesian, get_cells_ghost_2D)
{
  Grid3 const g(8,8,0);
  EXPECT_EQ(get_cells(GHOST, g, {-1,-1,0}), Subgrid3({0,0,0}, {1,1,0}));
  EXPECT_EQ(get_cells(GHOST, g, { 0,-1,0}), Subgrid3({1,0,0}, {7,1,0}));
  EXPECT_EQ(get_cells(GHOST, g, { 1,-1,0}), Subgrid3({7,0,0}, {8,1,0}));
  EXPECT_EQ(get_cells(GHOST, g, {-1, 0,0}), Subgrid3({0,1,0}, {1,7,0}));
  EXPECT_EQ(get_cells(GHOST, g, { 1, 0,0}), Subgrid3({7,1,0}, {8,7,0}));
  EXPECT_EQ(get_cells(GHOST, g, {-1, 1,0}), Subgrid3({0,7,0}, {1,8,0}));
  EXPECT_EQ(get_cells(GHOST, g, { 0, 1,0}), Subgrid3({1,7,0}, {7,8,0}));
  EXPECT_EQ(get_cells(GHOST, g, { 1, 1,0}), Subgrid3({7,7,0}, {8,8,0}));
}

TEST(cartesian, get_cells_ghost_3D)
{
  Grid3 const g(8,8,8);
  EXPECT_EQ(get_cells(GHOST, g, {-1,-1,-1}), Subgrid3({0,0,0}, {1,1,1}));
  EXPECT_EQ(get_cells(GHOST, g, { 0,-1,-1}), Subgrid3({1,0,0}, {7,1,1}));
  EXPECT_EQ(get_cells(GHOST, g, { 1,-1,-1}), Subgrid3({7,0,0}, {8,1,1}));
  EXPECT_EQ(get_cells(GHOST, g, {-1, 0,-1}), Subgrid3({0,1,0}, {1,7,1}));
  EXPECT_EQ(get_cells(GHOST, g, { 0, 0,-1}), Subgrid3({1,1,0}, {7,7,1}));
  EXPECT_EQ(get_cells(GHOST, g, { 1, 0,-1}), Subgrid3({7,1,0}, {8,7,1}));
  EXPECT_EQ(get_cells(GHOST, g, {-1, 1,-1}), Subgrid3({0,7,0}, {1,8,1}));
  EXPECT_EQ(get_cells(GHOST, g, { 0, 1,-1}), Subgrid3({1,7,0}, {7,8,1}));
  EXPECT_EQ(get_cells(GHOST, g, { 1, 1,-1}), Subgrid3({7,7,0}, {8,8,1}));
  EXPECT_EQ(get_cells(GHOST, g, {-1,-1, 0}), Subgrid3({0,0,1}, {1,1,7}));
  EXPECT_EQ(get_cells(GHOST, g, { 0,-1, 0}), Subgrid3({1,0,1}, {7,1,7}));
  EXPECT_EQ(get_cells(GHOST, g, { 1,-1, 0}), Subgrid3({7,0,1}, {8,1,7}));
  EXPECT_EQ(get_cells(GHOST, g, {-1, 0, 0}), Subgrid3({0,1,1}, {1,7,7}));
  EXPECT_EQ(get_cells(GHOST, g, { 1, 0, 0}), Subgrid3({7,1,1}, {8,7,7}));
  EXPECT_EQ(get_cells(GHOST, g, {-1, 1, 0}), Subgrid3({0,7,1}, {1,8,7}));
  EXPECT_EQ(get_cells(GHOST, g, { 0, 1, 0}), Subgrid3({1,7,1}, {7,8,7}));
  EXPECT_EQ(get_cells(GHOST, g, { 1, 1, 0}), Subgrid3({7,7,1}, {8,8,7}));
  EXPECT_EQ(get_cells(GHOST, g, {-1,-1, 1}), Subgrid3({0,0,7}, {1,1,8}));
  EXPECT_EQ(get_cells(GHOST, g, { 0,-1, 1}), Subgrid3({1,0,7}, {7,1,8}));
  EXPECT_EQ(get_cells(GHOST, g, { 1,-1, 1}), Subgrid3({7,0,7}, {8,1,8}));
  EXPECT_EQ(get_cells(GHOST, g, {-1, 0, 1}), Subgrid3({0,1,7}, {1,7,8}));
  EXPECT_EQ(get_cells(GHOST, g, { 0, 0, 1}), Subgrid3({1,1,7}, {7,7,8}));
  EXPECT_EQ(get_cells(GHOST, g, { 1, 0, 1}), Subgrid3({7,1,7}, {8,7,8}));
  EXPECT_EQ(get_cells(GHOST, g, {-1, 1, 1}), Subgrid3({0,7,7}, {1,8,8}));
  EXPECT_EQ(get_cells(GHOST, g, { 0, 1, 1}), Subgrid3({1,7,7}, {7,8,8}));
  EXPECT_EQ(get_cells(GHOST, g, { 1, 1, 1}), Subgrid3({7,7,7}, {8,8,8}));
}

TEST(cartesian, get_cells_owned_1D)
{
  Grid3 const g(8,0,0);
  EXPECT_EQ(get_cells(OWNED, g, {-1,0,0}), Subgrid3({1,0,0}, {2,0,0}));
  EXPECT_EQ(get_cells(OWNED, g, { 1,0,0}), Subgrid3({6,0,0}, {7,0,0}));
}

TEST(cartesian, get_cells_owned_2D)
{
  Grid3 const g(8,8,0);
  EXPECT_EQ(get_cells(OWNED, g, {-1,-1,0}), Subgrid3({1,1,0},{2,2,0}));
  EXPECT_EQ(get_cells(OWNED, g, { 0,-1,0}), Subgrid3({1,1,0},{7,2,0}));
  EXPECT_EQ(get_cells(OWNED, g, { 1,-1,0}), Subgrid3({6,1,0},{7,2,0}));
  EXPECT_EQ(get_cells(OWNED, g, {-1, 0,0}), Subgrid3({1,1,0},{2,7,0}));
  EXPECT_EQ(get_cells(OWNED, g, { 1, 0,0}), Subgrid3({6,1,0},{7,7,0}));
  EXPECT_EQ(get_cells(OWNED, g, {-1, 1,0}), Subgrid3({1,6,0},{2,7,0}));
  EXPECT_EQ(get_cells(OWNED, g, { 0, 1,0}), Subgrid3({1,6,0},{7,7,0}));
  EXPECT_EQ(get_cells(OWNED, g, { 1, 1,0}), Subgrid3({6,6,0},{7,7,0}));
}

TEST(cartesian, get_cells_owned_3D)
{
  Grid3 const g(8,8,8);
  EXPECT_EQ(get_cells(OWNED, g, {-1,-1,-1}), Subgrid3({1,1,1},{2,2,2}));
  EXPECT_EQ(get_cells(OWNED, g, { 0,-1,-1}), Subgrid3({1,1,1},{7,2,2}));
  EXPECT_EQ(get_cells(OWNED, g, { 1,-1,-1}), Subgrid3({6,1,1},{7,2,2}));
  EXPECT_EQ(get_cells(OWNED, g, {-1, 0,-1}), Subgrid3({1,1,1},{2,7,2}));
  EXPECT_EQ(get_cells(OWNED, g, { 1, 0,-1}), Subgrid3({6,1,1},{7,7,2}));
  EXPECT_EQ(get_cells(OWNED, g, {-1, 1,-1}), Subgrid3({1,6,1},{2,7,2}));
  EXPECT_EQ(get_cells(OWNED, g, { 0, 1,-1}), Subgrid3({1,6,1},{7,7,2}));
  EXPECT_EQ(get_cells(OWNED, g, { 1, 1,-1}), Subgrid3({6,6,1},{7,7,2}));
  EXPECT_EQ(get_cells(OWNED, g, {-1,-1, 0}), Subgrid3({1,1,1},{2,2,7}));
  EXPECT_EQ(get_cells(OWNED, g, { 0,-1, 0}), Subgrid3({1,1,1},{7,2,7}));
  EXPECT_EQ(get_cells(OWNED, g, { 1,-1, 0}), Subgrid3({6,1,1},{7,2,7}));
  EXPECT_EQ(get_cells(OWNED, g, {-1, 0, 0}), Subgrid3({1,1,1},{2,7,7}));
  EXPECT_EQ(get_cells(OWNED, g, { 1, 0, 0}), Subgrid3({6,1,1},{7,7,7}));
  EXPECT_EQ(get_cells(OWNED, g, {-1, 1, 0}), Subgrid3({1,6,1},{2,7,7}));
  EXPECT_EQ(get_cells(OWNED, g, { 0, 1, 0}), Subgrid3({1,6,1},{7,7,7}));
  EXPECT_EQ(get_cells(OWNED, g, { 1, 1, 0}), Subgrid3({6,6,1},{7,7,7}));
  EXPECT_EQ(get_cells(OWNED, g, {-1,-1, 1}), Subgrid3({1,1,6},{2,2,7}));
  EXPECT_EQ(get_cells(OWNED, g, { 0,-1, 1}), Subgrid3({1,1,6},{7,2,7}));
  EXPECT_EQ(get_cells(OWNED, g, { 1,-1, 1}), Subgrid3({6,1,6},{7,2,7}));
  EXPECT_EQ(get_cells(OWNED, g, {-1, 0, 1}), Subgrid3({1,1,6},{2,7,7}));
  EXPECT_EQ(get_cells(OWNED, g, { 1, 0, 1}), Subgrid3({6,1,6},{7,7,7}));
  EXPECT_EQ(get_cells(OWNED, g, {-1, 1, 1}), Subgrid3({1,6,6},{2,7,7}));
  EXPECT_EQ(get_cells(OWNED, g, { 0, 1, 1}), Subgrid3({1,6,6},{7,7,7}));
  EXPECT_EQ(get_cells(OWNED, g, { 1, 1, 1}), Subgrid3({6,6,6},{7,7,7}));
}

TEST(cartesian, get_fine_to_coarse_cells_owned_1D)
{
  Grid3 const g(8,0,0);
  EXPECT_EQ(get_fine_to_coarse_cells(OWNED, g, {-1,0,0}), Subgrid3({1,0,0}, {3,0,0}));
  EXPECT_EQ(get_fine_to_coarse_cells(OWNED, g, { 1,0,0}), Subgrid3({5,0,0}, {7,0,0}));
}

TEST(cartesian, get_fine_to_coarse_cells_owned_2D)
{
  Grid3 const g(8,8,0);
  EXPECT_EQ(get_fine_to_coarse_cells(OWNED, g, {-1,-1,0}), Subgrid3({1,1,0}, {3,3,0}));
  EXPECT_EQ(get_fine_to_coarse_cells(OWNED, g, { 0,-1,0}), Subgrid3({1,1,0}, {7,3,0}));
  EXPECT_EQ(get_fine_to_coarse_cells(OWNED, g, { 1,-1,0}), Subgrid3({5,1,0}, {7,3,0}));
  EXPECT_EQ(get_fine_to_coarse_cells(OWNED, g, {-1, 0,0}), Subgrid3({1,1,0}, {3,7,0}));
  EXPECT_EQ(get_fine_to_coarse_cells(OWNED, g, { 1, 0,0}), Subgrid3({5,1,0}, {7,7,0}));
  EXPECT_EQ(get_fine_to_coarse_cells(OWNED, g, {-1, 1,0}), Subgrid3({1,5,0}, {3,7,0}));
  EXPECT_EQ(get_fine_to_coarse_cells(OWNED, g, { 0, 1,0}), Subgrid3({1,5,0}, {7,7,0}));
  EXPECT_EQ(get_fine_to_coarse_cells(OWNED, g, { 1, 1,0}), Subgrid3({5,5,0}, {7,7,0}));
}

//TODO: get_fine_to_coarse_cells_owned_3D

TEST(cartesian, get_coarse_to_fine_cells_ghost_1D)
{
  Grid3 const g(8,0,0);
  EXPECT_EQ(get_coarse_to_fine_cells(GHOST, g, {1,0,0}, {-1,0,0}), Subgrid3({0,0,0}, {1,0,0}));
  EXPECT_EQ(get_coarse_to_fine_cells(GHOST, g, {0,0,0}, { 1,0,0}), Subgrid3({7,0,0}, {8,0,0}));
}

TEST(cartesian, get_coarse_to_fine_cells_ghost_2D)
{
  Grid3 const g(8,8,0);
  // XMIN face for the two children near it
  EXPECT_EQ(get_coarse_to_fine_cells(GHOST, g, {1,0,0}, {-1,0,0}), Subgrid3({0,1,0}, {1,4,0}));
  EXPECT_EQ(get_coarse_to_fine_cells(GHOST, g, {1,1,0}, {-1,0,0}), Subgrid3({0,4,0}, {1,7,0}));
  //TODO: fill this in
}

TEST(cartesian, get_coarse_to_fine_cells_ghost_3D)
{
  Grid3 const g(8,8,8);
  // XMIN face for the four children near it
  EXPECT_EQ(get_coarse_to_fine_cells(GHOST, g, {1,0,0}, {-1,0,0}), Subgrid3({0,1,1}, {1,4,4}));
  EXPECT_EQ(get_coarse_to_fine_cells(GHOST, g, {1,1,0}, {-1,0,0}), Subgrid3({0,4,1}, {1,7,4}));
  EXPECT_EQ(get_coarse_to_fine_cells(GHOST, g, {1,0,1}, {-1,0,0}), Subgrid3({0,1,4}, {1,4,7}));
  EXPECT_EQ(get_coarse_to_fine_cells(GHOST, g, {1,1,1}, {-1,0,0}), Subgrid3({0,4,4}, {1,7,7}));
  //TODO: fill this in
}

//TODO: get_coarse_to_fine_cells_owned_1D
//TODO: get_coarse_to_fine_cells_owned_2D
//TODO: get_coarse_to_fine_cells_owned_3D
