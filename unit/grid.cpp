#include "gtest/gtest.h"

#include "dgt_grid.hpp"
#include "dgt_tree.hpp"

TEST(grid, get_dim) {
  ASSERT_EQ(dgt::get_dim(p3a::grid3(2,0,0)), 1);
  ASSERT_EQ(dgt::get_dim(p3a::grid3(2,2,0)), 2);
  ASSERT_EQ(dgt::get_dim(p3a::grid3(2,2,2)), 3);
}

TEST(grid, grid_generalize) {
  ASSERT_EQ(dgt::generalize(p3a::grid3(4,0,0)).extents(), p3a::vector3<int>(4,1,1));
  ASSERT_EQ(dgt::generalize(p3a::grid3(4,4,0)).extents(), p3a::vector3<int>(4,4,1));
  ASSERT_EQ(dgt::generalize(p3a::grid3(4,4,4)).extents(), p3a::vector3<int>(4,4,4));
}

TEST(grid, subgrid_generalize) {
  ASSERT_EQ(dgt::generalize(p3a::subgrid3({2,0,0}, {10,0,0})).lower(), p3a::vector3<int>(2,0,0));
  ASSERT_EQ(dgt::generalize(p3a::subgrid3({2,0,0}, {10,0,0})).upper(), p3a::vector3<int>(10,1,1));
  ASSERT_EQ(dgt::generalize(p3a::subgrid3({2,2,0}, {10,10,0})).lower(), p3a::vector3<int>(2,2,0));
  ASSERT_EQ(dgt::generalize(p3a::subgrid3({2,2,0}, {10,10,0})).upper(), p3a::vector3<int>(10,10,1));
  ASSERT_EQ(dgt::generalize(p3a::subgrid3({2,2,2}, {10,10,10})).lower(), p3a::vector3<int>(2,2,2));
  ASSERT_EQ(dgt::generalize(p3a::subgrid3({2,2,2}, {10,10,10})).upper(), p3a::vector3<int>(10,10,10));
}

TEST(grid, child_grid) {
  ASSERT_EQ(dgt::get_child_grid(1).extents(), p3a::vector3<int>(2,0,0));
  ASSERT_EQ(dgt::get_child_grid(2).extents(), p3a::vector3<int>(2,2,0));
  ASSERT_EQ(dgt::get_child_grid(3).extents(), p3a::vector3<int>(2,2,2));
}

TEST(grid, side_grid) {
  ASSERT_EQ(dgt::get_side_grid(p3a::grid3(10,0,0), dgt::X).extents(), p3a::vector3<int>(11,0,0));
  ASSERT_EQ(dgt::get_side_grid(p3a::grid3(10,10,0), dgt::X).extents(), p3a::vector3<int>(11,10,0));
  ASSERT_EQ(dgt::get_side_grid(p3a::grid3(10,10,0), dgt::Y).extents(), p3a::vector3<int>(10,11,0));
  ASSERT_EQ(dgt::get_side_grid(p3a::grid3(10,10,10), dgt::X).extents(), p3a::vector3<int>(11,10,10));
  ASSERT_EQ(dgt::get_side_grid(p3a::grid3(10,10,10), dgt::Y).extents(), p3a::vector3<int>(10,11,10));
  ASSERT_EQ(dgt::get_side_grid(p3a::grid3(10,10,10), dgt::Z).extents(), p3a::vector3<int>(10,10,11));
}

TEST(grid, intr_sides) {
  ASSERT_EQ(dgt::get_intr_sides(p3a::grid3(10,0,0), dgt::X).lower(), p3a::vector3<int>(1,0,0));
  ASSERT_EQ(dgt::get_intr_sides(p3a::grid3(10,0,0), dgt::X).upper(), p3a::vector3<int>(10,0,0));
  ASSERT_EQ(dgt::get_intr_sides(p3a::grid3(10,10,0), dgt::X).lower(), p3a::vector3<int>(1,0,0));
  ASSERT_EQ(dgt::get_intr_sides(p3a::grid3(10,10,0), dgt::X).upper(), p3a::vector3<int>(10,10,0));
  ASSERT_EQ(dgt::get_intr_sides(p3a::grid3(10,10,0), dgt::Y).lower(), p3a::vector3<int>(0,1,0));
  ASSERT_EQ(dgt::get_intr_sides(p3a::grid3(10,10,0), dgt::Y).upper(), p3a::vector3<int>(10,10,0));
  ASSERT_EQ(dgt::get_intr_sides(p3a::grid3(10,10,10), dgt::X).lower(), p3a::vector3<int>(1,0,0));
  ASSERT_EQ(dgt::get_intr_sides(p3a::grid3(10,10,10), dgt::X).upper(), p3a::vector3<int>(10,10,10));
  ASSERT_EQ(dgt::get_intr_sides(p3a::grid3(10,10,10), dgt::Y).lower(), p3a::vector3<int>(0,1,0));
  ASSERT_EQ(dgt::get_intr_sides(p3a::grid3(10,10,10), dgt::Y).upper(), p3a::vector3<int>(10,10,10));
  ASSERT_EQ(dgt::get_intr_sides(p3a::grid3(10,10,10), dgt::Z).lower(), p3a::vector3<int>(0,0,1));
  ASSERT_EQ(dgt::get_intr_sides(p3a::grid3(10,10,10), dgt::Z).upper(), p3a::vector3<int>(10,10,10));
}

TEST(grid, viz_cell_grid) {
  ASSERT_EQ(dgt::get_viz_cell_grid({2,0,0}, 0).extents(), p3a::vector3<int>(2,0,0));
  ASSERT_EQ(dgt::get_viz_cell_grid({2,0,0}, 1).extents(), p3a::vector3<int>(4,0,0));
  ASSERT_EQ(dgt::get_viz_cell_grid({2,0,0}, 2).extents(), p3a::vector3<int>(6,0,0));
  ASSERT_EQ(dgt::get_viz_cell_grid({2,2,0}, 0).extents(), p3a::vector3<int>(2,2,0));
  ASSERT_EQ(dgt::get_viz_cell_grid({2,2,0}, 1).extents(), p3a::vector3<int>(4,4,0));
  ASSERT_EQ(dgt::get_viz_cell_grid({2,2,0}, 2).extents(), p3a::vector3<int>(6,6,0));
  ASSERT_EQ(dgt::get_viz_cell_grid({2,2,2}, 0).extents(), p3a::vector3<int>(2,2,2));
  ASSERT_EQ(dgt::get_viz_cell_grid({2,2,2}, 1).extents(), p3a::vector3<int>(4,4,4));
  ASSERT_EQ(dgt::get_viz_cell_grid({2,2,2}, 2).extents(), p3a::vector3<int>(6,6,6));
}

TEST(grid, viz_point_grid) {
  ASSERT_EQ(dgt::get_viz_point_grid({2,0,0}, 0).extents(), p3a::vector3<int>(3,0,0));
  ASSERT_EQ(dgt::get_viz_point_grid({2,0,0}, 1).extents(), p3a::vector3<int>(5,0,0));
  ASSERT_EQ(dgt::get_viz_point_grid({2,0,0}, 2).extents(), p3a::vector3<int>(7,0,0));
  ASSERT_EQ(dgt::get_viz_point_grid({2,2,0}, 0).extents(), p3a::vector3<int>(3,3,0));
  ASSERT_EQ(dgt::get_viz_point_grid({2,2,0}, 1).extents(), p3a::vector3<int>(5,5,0));
  ASSERT_EQ(dgt::get_viz_point_grid({2,2,0}, 2).extents(), p3a::vector3<int>(7,7,0));
  ASSERT_EQ(dgt::get_viz_point_grid({2,2,2}, 0).extents(), p3a::vector3<int>(3,3,3));
  ASSERT_EQ(dgt::get_viz_point_grid({2,2,2}, 1).extents(), p3a::vector3<int>(5,5,5));
  ASSERT_EQ(dgt::get_viz_point_grid({2,2,2}, 2).extents(), p3a::vector3<int>(7,7,7));
}

TEST(grid, contains_subgrid) {
  ASSERT_EQ(dgt::contains(p3a::grid3(4,4,4), p3a::subgrid3({0,0,0}, {2,2,2})), true);
}

TEST(grid, contains) {
  ASSERT_EQ(dgt::contains(2, p3a::subgrid3({0,0,0},{2,0,0}), {2, {1,0,0}}), true);
  ASSERT_EQ(dgt::contains(2, p3a::subgrid3({0,0,0},{2,0,0}), {2, {2,0,0}}), false);
  ASSERT_EQ(dgt::contains(2, p3a::subgrid3({0,0,0},{2,0,0}), {2, {1,0,1}}), false);
  ASSERT_EQ(dgt::contains(2, p3a::subgrid3({0,0,0},{2,2,0}), {2, {1,0,0}}), true);
  ASSERT_EQ(dgt::contains(2, p3a::subgrid3({0,0,0},{2,2,0}), {2, {2,0,0}}), false);
  ASSERT_EQ(dgt::contains(2, p3a::subgrid3({0,0,0},{2,2,0}), {2, {1,0,1}}), false);
  ASSERT_EQ(dgt::contains(2, p3a::subgrid3({0,0,0},{2,2,2}), {2, {1,0,0}}), true);
  ASSERT_EQ(dgt::contains(2, p3a::subgrid3({0,0,0},{2,2,2}), {2, {0,1,0}}), true);
  ASSERT_EQ(dgt::contains(2, p3a::subgrid3({0,0,0},{2,2,2}), {2, {2,0,1}}), false);
}

TEST(grid, num_local) {
  ASSERT_EQ(dgt::get_num_local(3, 2, 0), 2);
  ASSERT_EQ(dgt::get_num_local(3, 2, 1), 1);
}

TEST(grid, local_offset) {
  ASSERT_EQ(dgt::get_local_offset(3, 2, 0), 0);
  ASSERT_EQ(dgt::get_local_offset(3, 2, 1), 2);
}

TEST(grid, pt_equality) {
  dgt::Point a = {1, {1,2,3}};
  dgt::Point b = {1, {1,2,3}};
  dgt::Point c = {2, {1,2,3}};
  ASSERT_TRUE(a == b);
  ASSERT_TRUE(a != c);
}
