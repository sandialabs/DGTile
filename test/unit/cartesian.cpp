#include <dgt_cartesian.hpp>

#include <gtest/gtest.h>

using namespace dgt;

TEST(cartesian, grid_definitions)
{
  EXPECT_EQ(child_grid, Grid3(2,2,2));
  EXPECT_EQ(offset_grid, Subgrid3({-1,-1,-1},{2,2,2}));
  EXPECT_EQ(fine_offset_grid, Subgrid3({-1,-1,-1},{3,3,3}));
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

TEST(cartesian, get_adj_meta_ijk_equal_1D)
{
  int const dim = 1;
  EXPECT_EQ(get_adj_meta_ijk(dim, {0,0,0}, 0), Vec3<std::int8_t>(2,0,0));
  EXPECT_EQ(get_adj_meta_ijk(dim, {2,0,0}, 0), Vec3<std::int8_t>(0,0,0));
}

TEST(cartesian, get_adj_meta_ijk_equal_2D)
{
  int const dim = 2;
  EXPECT_EQ(get_adj_meta_ijk(dim, {0,0,0}, 0), Vec3<std::int8_t>(2,2,0));
  EXPECT_EQ(get_adj_meta_ijk(dim, {1,0,0}, 0), Vec3<std::int8_t>(1,2,0));
  EXPECT_EQ(get_adj_meta_ijk(dim, {2,0,0}, 0), Vec3<std::int8_t>(0,2,0));
  EXPECT_EQ(get_adj_meta_ijk(dim, {0,1,0}, 0), Vec3<std::int8_t>(2,1,0));
  EXPECT_EQ(get_adj_meta_ijk(dim, {2,1,0}, 0), Vec3<std::int8_t>(0,1,0));
  EXPECT_EQ(get_adj_meta_ijk(dim, {0,2,0}, 0), Vec3<std::int8_t>(2,0,0));
  EXPECT_EQ(get_adj_meta_ijk(dim, {1,2,0}, 0), Vec3<std::int8_t>(1,0,0));
  EXPECT_EQ(get_adj_meta_ijk(dim, {2,2,0}, 0), Vec3<std::int8_t>(0,0,0));
}

TEST(cartesian, get_adj_meta_ijk_equal_3D)
{
  int const dim = 3;
  EXPECT_EQ(get_adj_meta_ijk(dim, {0,0,0}, 0), Vec3<std::int8_t>(2,2,2));
  EXPECT_EQ(get_adj_meta_ijk(dim, {1,0,0}, 0), Vec3<std::int8_t>(1,2,2));
  EXPECT_EQ(get_adj_meta_ijk(dim, {2,0,0}, 0), Vec3<std::int8_t>(0,2,2));
  EXPECT_EQ(get_adj_meta_ijk(dim, {0,1,0}, 0), Vec3<std::int8_t>(2,1,2));
  EXPECT_EQ(get_adj_meta_ijk(dim, {1,1,0}, 0), Vec3<std::int8_t>(1,1,2));
  EXPECT_EQ(get_adj_meta_ijk(dim, {2,1,0}, 0), Vec3<std::int8_t>(0,1,2));
  EXPECT_EQ(get_adj_meta_ijk(dim, {0,2,0}, 0), Vec3<std::int8_t>(2,0,2));
  EXPECT_EQ(get_adj_meta_ijk(dim, {1,2,0}, 0), Vec3<std::int8_t>(1,0,2));
  EXPECT_EQ(get_adj_meta_ijk(dim, {2,2,0}, 0), Vec3<std::int8_t>(0,0,2));
  EXPECT_EQ(get_adj_meta_ijk(dim, {0,0,1}, 0), Vec3<std::int8_t>(2,2,1));
  EXPECT_EQ(get_adj_meta_ijk(dim, {1,0,1}, 0), Vec3<std::int8_t>(1,2,1));
  EXPECT_EQ(get_adj_meta_ijk(dim, {2,0,1}, 0), Vec3<std::int8_t>(0,2,1));
  EXPECT_EQ(get_adj_meta_ijk(dim, {0,1,1}, 0), Vec3<std::int8_t>(2,1,1));
  EXPECT_EQ(get_adj_meta_ijk(dim, {2,1,1}, 0), Vec3<std::int8_t>(0,1,1));
  EXPECT_EQ(get_adj_meta_ijk(dim, {0,2,1}, 0), Vec3<std::int8_t>(2,0,1));
  EXPECT_EQ(get_adj_meta_ijk(dim, {1,2,1}, 0), Vec3<std::int8_t>(1,0,1));
  EXPECT_EQ(get_adj_meta_ijk(dim, {2,2,1}, 0), Vec3<std::int8_t>(0,0,1));
  EXPECT_EQ(get_adj_meta_ijk(dim, {0,0,2}, 0), Vec3<std::int8_t>(2,2,0));
  EXPECT_EQ(get_adj_meta_ijk(dim, {1,0,2}, 0), Vec3<std::int8_t>(1,2,0));
  EXPECT_EQ(get_adj_meta_ijk(dim, {2,0,2}, 0), Vec3<std::int8_t>(0,2,0));
  EXPECT_EQ(get_adj_meta_ijk(dim, {0,1,2}, 0), Vec3<std::int8_t>(2,1,0));
  EXPECT_EQ(get_adj_meta_ijk(dim, {1,1,2}, 0), Vec3<std::int8_t>(1,1,0));
  EXPECT_EQ(get_adj_meta_ijk(dim, {2,1,2}, 0), Vec3<std::int8_t>(0,1,0));
  EXPECT_EQ(get_adj_meta_ijk(dim, {0,2,2}, 0), Vec3<std::int8_t>(2,0,0));
  EXPECT_EQ(get_adj_meta_ijk(dim, {1,2,2}, 0), Vec3<std::int8_t>(1,0,0));
  EXPECT_EQ(get_adj_meta_ijk(dim, {2,2,2}, 0), Vec3<std::int8_t>(0,0,0));
}

TEST(cartesian, get_cells_ghost_equal_1D)
{
  Grid3 const g(8,0,0);
  EXPECT_EQ(get_cells(GHOST, g, {0,0,0}, 0), Subgrid3({0,0,0}, {1,0,0}));
  EXPECT_EQ(get_cells(GHOST, g, {2,0,0}, 0), Subgrid3({7,0,0}, {8,0,0}));
}

TEST(cartesian, get_cells_owned_equal_1D)
{
  Grid3 const g(8,0,0);
  EXPECT_EQ(get_cells(OWNED, g, {0,0,0}, 0), Subgrid3({1,0,0}, {2,0,0}));
  EXPECT_EQ(get_cells(OWNED, g, {2,0,0}, 0), Subgrid3({6,0,0}, {7,0,0}));
}

TEST(cartesian, get_cells_ghost_equal_2D)
{
  Grid3 const g(8,8,0);
  EXPECT_EQ(get_cells(GHOST, g, {0,0,0}, 0), Subgrid3({0,0,0}, {1,1,0}));
  EXPECT_EQ(get_cells(GHOST, g, {1,0,0}, 0), Subgrid3({1,0,0}, {7,1,0}));
  EXPECT_EQ(get_cells(GHOST, g, {2,0,0}, 0), Subgrid3({7,0,0}, {8,1,0}));
  EXPECT_EQ(get_cells(GHOST, g, {0,1,0}, 0), Subgrid3({0,1,0}, {1,7,0}));
  EXPECT_EQ(get_cells(GHOST, g, {2,1,0}, 0), Subgrid3({7,1,0}, {8,7,0}));
  EXPECT_EQ(get_cells(GHOST, g, {0,2,0}, 0), Subgrid3({0,7,0}, {1,8,0}));
  EXPECT_EQ(get_cells(GHOST, g, {1,2,0}, 0), Subgrid3({1,7,0}, {7,8,0}));
  EXPECT_EQ(get_cells(GHOST, g, {2,2,0}, 0), Subgrid3({7,7,0}, {8,8,0}));
}

TEST(cartesian, get_cells_owned_equal_2D)
{
  Grid3 const g(8,8,0);
  EXPECT_EQ(get_cells(OWNED, g, {0,0,0}, 0), Subgrid3({1,1,0}, {2,2,0}));
  EXPECT_EQ(get_cells(OWNED, g, {1,0,0}, 0), Subgrid3({2,1,0}, {6,2,0}));
  EXPECT_EQ(get_cells(OWNED, g, {2,0,0}, 0), Subgrid3({6,1,0}, {7,2,0}));
  EXPECT_EQ(get_cells(OWNED, g, {0,1,0}, 0), Subgrid3({1,2,0}, {2,6,0}));
  EXPECT_EQ(get_cells(OWNED, g, {2,1,0}, 0), Subgrid3({6,2,0}, {7,6,0}));
  EXPECT_EQ(get_cells(OWNED, g, {0,2,0}, 0), Subgrid3({1,6,0}, {2,7,0}));
  EXPECT_EQ(get_cells(OWNED, g, {1,2,0}, 0), Subgrid3({2,6,0}, {6,7,0}));
  EXPECT_EQ(get_cells(OWNED, g, {2,2,0}, 0), Subgrid3({6,6,0}, {7,7,0}));
}

TEST(cartesian, get_cells_ghost_equal_3D)
{
  Grid3 const g(8,8,8);
  EXPECT_EQ(get_cells(GHOST, g, {0,0,0}, 0), Subgrid3({0,0,0}, {1,1,1}));
  EXPECT_EQ(get_cells(GHOST, g, {1,0,0}, 0), Subgrid3({1,0,0}, {7,1,1}));
  EXPECT_EQ(get_cells(GHOST, g, {2,0,0}, 0), Subgrid3({7,0,0}, {8,1,1}));
  EXPECT_EQ(get_cells(GHOST, g, {0,1,0}, 0), Subgrid3({0,1,0}, {1,7,1}));
  EXPECT_EQ(get_cells(GHOST, g, {1,1,0}, 0), Subgrid3({1,1,0}, {7,7,1}));
  EXPECT_EQ(get_cells(GHOST, g, {2,1,0}, 0), Subgrid3({7,1,0}, {8,7,1}));
  EXPECT_EQ(get_cells(GHOST, g, {0,2,0}, 0), Subgrid3({0,7,0}, {1,8,1}));
  EXPECT_EQ(get_cells(GHOST, g, {1,2,0}, 0), Subgrid3({1,7,0}, {7,8,1}));
  EXPECT_EQ(get_cells(GHOST, g, {2,2,0}, 0), Subgrid3({7,7,0}, {8,8,1}));
  EXPECT_EQ(get_cells(GHOST, g, {0,0,1}, 0), Subgrid3({0,0,1}, {1,1,7}));
  EXPECT_EQ(get_cells(GHOST, g, {1,0,1}, 0), Subgrid3({1,0,1}, {7,1,7}));
  EXPECT_EQ(get_cells(GHOST, g, {2,0,1}, 0), Subgrid3({7,0,1}, {8,1,7}));
  EXPECT_EQ(get_cells(GHOST, g, {0,1,1}, 0), Subgrid3({0,1,1}, {1,7,7}));
  EXPECT_EQ(get_cells(GHOST, g, {2,1,1}, 0), Subgrid3({7,1,1}, {8,7,7}));
  EXPECT_EQ(get_cells(GHOST, g, {0,2,1}, 0), Subgrid3({0,7,1}, {1,8,7}));
  EXPECT_EQ(get_cells(GHOST, g, {1,2,1}, 0), Subgrid3({1,7,1}, {7,8,7}));
  EXPECT_EQ(get_cells(GHOST, g, {2,2,1}, 0), Subgrid3({7,7,1}, {8,8,7}));
  EXPECT_EQ(get_cells(GHOST, g, {0,0,2}, 0), Subgrid3({0,0,7}, {1,1,8}));
  EXPECT_EQ(get_cells(GHOST, g, {1,0,2}, 0), Subgrid3({1,0,7}, {7,1,8}));
  EXPECT_EQ(get_cells(GHOST, g, {2,0,2}, 0), Subgrid3({7,0,7}, {8,1,8}));
  EXPECT_EQ(get_cells(GHOST, g, {0,1,2}, 0), Subgrid3({0,1,7}, {1,7,8}));
  EXPECT_EQ(get_cells(GHOST, g, {1,1,2}, 0), Subgrid3({1,1,7}, {7,7,8}));
  EXPECT_EQ(get_cells(GHOST, g, {2,1,2}, 0), Subgrid3({7,1,7}, {8,7,8}));
  EXPECT_EQ(get_cells(GHOST, g, {0,2,2}, 0), Subgrid3({0,7,7}, {1,8,8}));
  EXPECT_EQ(get_cells(GHOST, g, {1,2,2}, 0), Subgrid3({1,7,7}, {7,8,8}));
  EXPECT_EQ(get_cells(GHOST, g, {2,2,2}, 0), Subgrid3({7,7,7}, {8,8,8}));
}
