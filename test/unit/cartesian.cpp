#include <dgt_cartesian.hpp>
#include <dgt_tree.hpp>

#include <gtest/gtest.h>

#include <dgt_print.hpp> // debug

using namespace dgt;

TEST(cartesian, grid_definitions)
{
  EXPECT_EQ(child_grid, Grid3(2,2,2));
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

TEST(cartesian, is_cell_ghost)
{
  EXPECT_EQ(is_cell_ghost(1, {10,0,0}, {0,0,0}), true);
  EXPECT_EQ(is_cell_ghost(1, {10,0,0}, {3,0,0}), false);
  EXPECT_EQ(is_cell_ghost(1, {10,0,0}, {9,0,0}), true);
  EXPECT_EQ(is_cell_ghost(2, {10,10,0}, {0,0,0}), true);
  EXPECT_EQ(is_cell_ghost(2, {10,10,0}, {0,1,0}), true);
  EXPECT_EQ(is_cell_ghost(2, {10,10,0}, {9,9,0}), true);
  EXPECT_EQ(is_cell_ghost(2, {10,10,0}, {9,6,0}), true);
  EXPECT_EQ(is_cell_ghost(2, {10,10,0}, {1,1,0}), false);
  EXPECT_EQ(is_cell_ghost(2, {10,10,0}, {5,5,0}), false);
  EXPECT_EQ(is_cell_ghost(3, {10,10,10}, {0,0,0}), true);
  EXPECT_EQ(is_cell_ghost(3, {10,10,10}, {0,1,8}), true);
  EXPECT_EQ(is_cell_ghost(3, {10,10,10}, {9,9,9}), true);
  EXPECT_EQ(is_cell_ghost(3, {10,10,10}, {1,1,1}), false);
  EXPECT_EQ(is_cell_ghost(3, {10,10,10}, {3,4,5}), false);
}

TEST(cartesian, get_child_ijk)
{
  EXPECT_EQ(get_child_ijk(0), Vec3<int>(0,0,0));
  EXPECT_EQ(get_child_ijk(1), Vec3<int>(1,0,0));
  EXPECT_EQ(get_child_ijk(2), Vec3<int>(0,1,0));
  EXPECT_EQ(get_child_ijk(3), Vec3<int>(1,1,0));
  EXPECT_EQ(get_child_ijk(4), Vec3<int>(0,0,1));
  EXPECT_EQ(get_child_ijk(5), Vec3<int>(1,0,1));
  EXPECT_EQ(get_child_ijk(6), Vec3<int>(0,1,1));
  EXPECT_EQ(get_child_ijk(7), Vec3<int>(1,1,1));
}

TEST(cartesian, get_which_child)
{
  EXPECT_EQ(get_which_child({0,0,0}), 0);
  EXPECT_EQ(get_which_child({1,0,0}), 1);
  EXPECT_EQ(get_which_child({0,1,0}), 2);
  EXPECT_EQ(get_which_child({1,1,0}), 3);
  EXPECT_EQ(get_which_child({0,0,1}), 4);
  EXPECT_EQ(get_which_child({1,0,1}), 5);
  EXPECT_EQ(get_which_child({0,1,1}), 6);
  EXPECT_EQ(get_which_child({1,1,1}), 7);
}

TEST(cartesian, get_axis_name)
{
  EXPECT_EQ(get_axis_name(X), "X");
  EXPECT_EQ(get_axis_name(Y), "Y");
  EXPECT_EQ(get_axis_name(Z), "Z");
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

TEST(cartesian, get_cells_equal_ghost_1D)
{
  Grid3 const cell_grid(6,0,0);
  EXPECT_EQ(get_cells(cell_grid, GHOST, tree::EQUAL, X, LEFT), Subgrid3({0,0,0}, {1,0,0}));
  EXPECT_EQ(get_cells(cell_grid, GHOST, tree::EQUAL, X, RIGHT), Subgrid3({5,0,0}, {6,0,0}));
}

TEST(cartesian, get_cells_equal_ghost_2D)
{
  Grid3 const cell_grid(6,6,0);
  EXPECT_EQ(get_cells(cell_grid, GHOST, tree::EQUAL, X, LEFT), Subgrid3({0,1,0}, {1,5,0}));
  EXPECT_EQ(get_cells(cell_grid, GHOST, tree::EQUAL, X, RIGHT), Subgrid3({5,1,0}, {6,5,0}));
  EXPECT_EQ(get_cells(cell_grid, GHOST, tree::EQUAL, Y, LEFT), Subgrid3({1,0,0}, {5,1,0}));
  EXPECT_EQ(get_cells(cell_grid, GHOST, tree::EQUAL, Y, RIGHT), Subgrid3({1,5,0}, {5,6,0}));
}

TEST(cartesian, get_cells_equal_ghost_3D)
{
  Grid3 const cell_grid(6,6,6);
  EXPECT_EQ(get_cells(cell_grid, GHOST, tree::EQUAL, X, LEFT), Subgrid3({0,1,1}, {1,5,5}));
  EXPECT_EQ(get_cells(cell_grid, GHOST, tree::EQUAL, X, RIGHT), Subgrid3({5,1,1}, {6,5,5}));
  EXPECT_EQ(get_cells(cell_grid, GHOST, tree::EQUAL, Y, LEFT), Subgrid3({1,0,1}, {5,1,5}));
  EXPECT_EQ(get_cells(cell_grid, GHOST, tree::EQUAL, Y, RIGHT), Subgrid3({1,5,1}, {5,6,5}));
  EXPECT_EQ(get_cells(cell_grid, GHOST, tree::EQUAL, Z, LEFT), Subgrid3({1,1,0}, {5,5,1}));
  EXPECT_EQ(get_cells(cell_grid, GHOST, tree::EQUAL, Z, RIGHT), Subgrid3({1,1,5}, {5,5,6}));
}

TEST(cartesian, get_cells_equal_owned_1D)
{
  Grid3 const cell_grid(6,0,0);
  EXPECT_EQ(get_cells(cell_grid, OWNED, tree::EQUAL, X, LEFT), Subgrid3({1,0,0}, {2,0,0}));
  EXPECT_EQ(get_cells(cell_grid, OWNED, tree::EQUAL, X, RIGHT), Subgrid3({4,0,0}, {5,0,0}));
}

TEST(cartesian, get_cells_equal_owned_2D)
{
  Grid3 const cell_grid(6,6,0);
  EXPECT_EQ(get_cells(cell_grid, OWNED, tree::EQUAL, X, LEFT), Subgrid3({1,1,0}, {2,5,0}));
  EXPECT_EQ(get_cells(cell_grid, OWNED, tree::EQUAL, X, RIGHT), Subgrid3({4,1,0}, {5,5,0}));
  EXPECT_EQ(get_cells(cell_grid, OWNED, tree::EQUAL, Y, LEFT), Subgrid3({1,1,0}, {5,2,0}));
  EXPECT_EQ(get_cells(cell_grid, OWNED, tree::EQUAL, Y, RIGHT), Subgrid3({1,4,0}, {5,5,0}));
}

TEST(cartesian, get_cells_equal_owned_3D)
{
  Grid3 const cell_grid(6,6,6);
  EXPECT_EQ(get_cells(cell_grid, OWNED, tree::EQUAL, X, LEFT), Subgrid3({1,1,1}, {2,5,5}));
  EXPECT_EQ(get_cells(cell_grid, OWNED, tree::EQUAL, X, RIGHT), Subgrid3({4,1,1}, {5,5,5}));
  EXPECT_EQ(get_cells(cell_grid, OWNED, tree::EQUAL, Y, LEFT), Subgrid3({1,1,1}, {5,2,5}));
  EXPECT_EQ(get_cells(cell_grid, OWNED, tree::EQUAL, Y, RIGHT), Subgrid3({1,4,1}, {5,5,5}));
  EXPECT_EQ(get_cells(cell_grid, OWNED, tree::EQUAL, Z, LEFT), Subgrid3({1,1,1}, {5,5,2}));
  EXPECT_EQ(get_cells(cell_grid, OWNED, tree::EQUAL, Z, RIGHT), Subgrid3({1,1,4}, {5,5,5}));
}
