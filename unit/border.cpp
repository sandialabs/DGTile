#include "gtest/gtest.h"

#include "dgt_border.hpp"

using namespace dgt;

TEST(grid, owned_adj_cells) {
  ASSERT_EQ(get_adj_cells(owned, grid3(20,20,20), X, left).lower(), vector3<int>(1,1,1));
  ASSERT_EQ(get_adj_cells(owned, grid3(20,20,20), X, left).upper(), vector3<int>(2,21,21));
  ASSERT_EQ(get_adj_cells(owned, grid3(20,20,20), X, right).lower(), vector3<int>(20,1,1));
  ASSERT_EQ(get_adj_cells(owned, grid3(20,20,20), X, right).upper(), vector3<int>(21,21,21));
  ASSERT_EQ(get_adj_cells(owned, grid3(20,20,20), Y, left).lower(), vector3<int>(1,1,1));
  ASSERT_EQ(get_adj_cells(owned, grid3(20,20,20), Y, left).upper(), vector3<int>(21,2,21));
  ASSERT_EQ(get_adj_cells(owned, grid3(20,20,20), Y, right).lower(), vector3<int>(1,20,1));
  ASSERT_EQ(get_adj_cells(owned, grid3(20,20,20), Y, right).upper(), vector3<int>(21,21,21));
  ASSERT_EQ(get_adj_cells(owned, grid3(20,20,20), Z, left).lower(), vector3<int>(1,1,1));
  ASSERT_EQ(get_adj_cells(owned, grid3(20,20,20), Z, left).upper(), vector3<int>(21,21,2));
  ASSERT_EQ(get_adj_cells(owned, grid3(20,20,20), Z, right).lower(), vector3<int>(1,1,20));
  ASSERT_EQ(get_adj_cells(owned, grid3(20,20,20), Z, right).upper(), vector3<int>(21,21,21));
}

TEST(grid, ghost_adj_cells) {
  ASSERT_EQ(get_adj_cells(ghost, grid3(20,20,20), X, left).lower(), vector3<int>(0,1,1));
  ASSERT_EQ(get_adj_cells(ghost, grid3(20,20,20), X, left).upper(), vector3<int>(1,21,21));
  ASSERT_EQ(get_adj_cells(ghost, grid3(20,20,20), X, right).lower(), vector3<int>(21,1,1));
  ASSERT_EQ(get_adj_cells(ghost, grid3(20,20,20), X, right).upper(), vector3<int>(22,21,21));
  ASSERT_EQ(get_adj_cells(ghost, grid3(20,20,20), Y, left).lower(), vector3<int>(1,0,1));
  ASSERT_EQ(get_adj_cells(ghost, grid3(20,20,20), Y, left).upper(), vector3<int>(21,1,21));
  ASSERT_EQ(get_adj_cells(ghost, grid3(20,20,20), Y, right).lower(), vector3<int>(1,21,1));
  ASSERT_EQ(get_adj_cells(ghost, grid3(20,20,20), Y, right).upper(), vector3<int>(21,22,21));
  ASSERT_EQ(get_adj_cells(ghost, grid3(20,20,20), Z, left).lower(), vector3<int>(1,1,0));
  ASSERT_EQ(get_adj_cells(ghost, grid3(20,20,20), Z, left).upper(), vector3<int>(21,21,1));
  ASSERT_EQ(get_adj_cells(ghost, grid3(20,20,20), Z, right).lower(), vector3<int>(1,1,21));
  ASSERT_EQ(get_adj_cells(ghost, grid3(20,20,20), Z, right).upper(), vector3<int>(21,21,22));
}
