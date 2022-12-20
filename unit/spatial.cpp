#include "gtest/gtest.h"

#include "dgt_spatial.hpp"

using namespace dgt;

TEST(spatial, dir_sign) {
  ASSERT_EQ(get_dir_sign(left), -1);
  ASSERT_EQ(get_dir_sign(right), 1);
}

TEST(spatial, invert_dir) {
  ASSERT_EQ(invert_dir(left), right);
  ASSERT_EQ(invert_dir(right), left);
}

TEST(spatial, sides_adj_cell) {
  ASSERT_EQ(get_sides_adj_cell({4,4,4}, X, left), vector3<int>(3,4,4));
  ASSERT_EQ(get_sides_adj_cell({4,4,4}, X, right), vector3<int>(4,4,4));
  ASSERT_EQ(get_sides_adj_cell({4,4,4}, Y, left), vector3<int>(4,3,4));
  ASSERT_EQ(get_sides_adj_cell({4,4,4}, Y, right), vector3<int>(4,4,4));
  ASSERT_EQ(get_sides_adj_cell({4,4,4}, Z, left), vector3<int>(4,4,3));
  ASSERT_EQ(get_sides_adj_cell({4,4,4}, Z, right), vector3<int>(4,4,4));
}

TEST(spatial, cells_adj_side) {
  ASSERT_EQ(get_cells_adj_side({4,4,4}, X, left), vector3<int>(4,4,4));
  ASSERT_EQ(get_cells_adj_side({4,4,4}, X, right), vector3<int>(5,4,4));
  ASSERT_EQ(get_cells_adj_side({4,4,4}, Y, left), vector3<int>(4,4,4));
  ASSERT_EQ(get_cells_adj_side({4,4,4}, Y, right), vector3<int>(4,5,4));
  ASSERT_EQ(get_cells_adj_side({4,4,4}, Z, left), vector3<int>(4,4,4));
  ASSERT_EQ(get_cells_adj_side({4,4,4}, Z, right), vector3<int>(4,4,5));
}

TEST(spatial, cells_adj_cell) {
  ASSERT_EQ(get_cells_adj_cell({4,4,4}, X, left), vector3<int>(3,4,4));
  ASSERT_EQ(get_cells_adj_cell({4,4,4}, X, right), vector3<int>(5,4,4));
  ASSERT_EQ(get_cells_adj_cell({4,4,4}, Y, left), vector3<int>(4,3,4));
  ASSERT_EQ(get_cells_adj_cell({4,4,4}, Y, right), vector3<int>(4,5,4));
  ASSERT_EQ(get_cells_adj_cell({4,4,4}, Z, left), vector3<int>(4,4,3));
  ASSERT_EQ(get_cells_adj_cell({4,4,4}, Z, right), vector3<int>(4,4,5));
}

TEST(spatial, border_ijk) {
  ASSERT_EQ(get_border_ijk({4,4,4}, X), vector3<int>(0,4,4));
  ASSERT_EQ(get_border_ijk({4,4,4}, Y), vector3<int>(4,0,4));
  ASSERT_EQ(get_border_ijk({4,4,4}, Z), vector3<int>(4,4,0));
}

TEST(spatial, dx) {
  ASSERT_EQ(get_dx({{0.,0.,0.},{1.,0.,0.}}, {10,0,0}), vector3<double>(.1,0.,0.));
  ASSERT_EQ(get_dx({{0.,0.,0.},{1.,1.,0.}}, {10,10,0}), vector3<double>(.1,.1,0.));
  ASSERT_EQ(get_dx({{0.,0.,0.},{1.,1.,1.}}, {10,10,10}), vector3<double>(.1,.1,.1));
}

TEST(spatial, block_domain_1D) {
  Point const base{1, {2,0,0}};
  Point const node{2, {2,0,0}};
  box3<double> const domain(vector3<double>(0,0,0), vector3<double>(1,0,0));
  box3<double> const exact(vector3<double>(.5,0,0), vector3<double>(.75,0,0));
  ASSERT_EQ(get_block_domain(base, node, domain), exact);
}

TEST(spatial, block_domain_2D) {
  Point const base{1, {2,2,0}};
  Point const node{2, {2,2,0}};
  box3<double> const domain(vector3<double>(0,0,0), vector3<double>(1,1,0));
  box3<double> const exact(vector3<double>(.5,.5,0), vector3<double>(.75,.75,0));
  ASSERT_EQ(get_block_domain(base, node, domain), exact);
}

TEST(spatial, block_domain_3D) {
  Point const base{1, {2,2,2}};
  Point const node{2, {2,2,2}};
  box3<double> const domain(vector3<double>(0,0,0), vector3<double>(1,1,1));
  box3<double> const exact(vector3<double>(.5,.5,.5), vector3<double>(.75,.75,.75));
  ASSERT_EQ(get_block_domain(base, node, domain), exact);
}

TEST(spatial, cell_center) {
  ASSERT_EQ(get_cell_center({2,0,0}, vector3<double>(0,0,0), vector3<double>(1,0,0)), vector3<double>(2.5,0,0));
  ASSERT_EQ(get_cell_center({2,2,0}, vector3<double>(0,0,0), vector3<double>(1,1,0)), vector3<double>(2.5,2.5,0));
  ASSERT_EQ(get_cell_center({2,2,2}, vector3<double>(0,0,0), vector3<double>(1,1,1)), vector3<double>(2.5,2.5,2.5));
}

TEST(spatial, x) {
  ASSERT_EQ(get_x({2,0,0}, vector3<double>(0,0,0), vector3<double>(1,0,0), vector3<double>(-1,0,0)), vector3<double>(2,0,0));
  ASSERT_EQ(get_x({2,0,0}, vector3<double>(0,0,0), vector3<double>(1,0,0), vector3<double>(1,0,0)), vector3<double>(3,0,0));
  ASSERT_EQ(get_x({2,2,0}, vector3<double>(0,0,0), vector3<double>(1,1,0), vector3<double>(-1,-1,0)), vector3<double>(2,2,0));
  ASSERT_EQ(get_x({2,2,0}, vector3<double>(0,0,0), vector3<double>(1,1,0), vector3<double>(1,1,0)), vector3<double>(3,3,0));
  ASSERT_EQ(get_x({2,2,2}, vector3<double>(0,0,0), vector3<double>(1,1,1), vector3<double>(-1,-1,-1)), vector3<double>(2,2,2));
  ASSERT_EQ(get_x({2,2,2}, vector3<double>(0,0,0), vector3<double>(1,1,1), vector3<double>(1,1,1)), vector3<double>(3,3,3));
}

TEST(spatial, volume) {
  ASSERT_EQ(get_volume(1, vector3<double>(10,0,0)), 10.);
  ASSERT_EQ(get_volume(2, vector3<double>(10,10,0)), 100.);
  ASSERT_EQ(get_volume(3, vector3<double>(10,10,10)), 1000.);
}

TEST(spatial, cell_detJ) {
  ASSERT_EQ(get_cell_detJ(1, vector3<double>(1,0,0)), 0.5);
  ASSERT_EQ(get_cell_detJ(2, vector3<double>(1,1,0)), 0.25);
  ASSERT_EQ(get_cell_detJ(3, vector3<double>(1,1,1)), 0.125);
}

TEST(spatial, side_detJ) {
  ASSERT_EQ(get_side_detJ(1, X, vector3<double>(1,0,0)), 1.);
  ASSERT_EQ(get_side_detJ(2, X, vector3<double>(1,1,0)), 0.5);
  ASSERT_EQ(get_side_detJ(2, Y, vector3<double>(1,1,0)), 0.5);
  ASSERT_EQ(get_side_detJ(3, X, vector3<double>(1,1,1)), 0.25);
  ASSERT_EQ(get_side_detJ(3, Y, vector3<double>(1,1,1)), 0.25);
  ASSERT_EQ(get_side_detJ(3, Z, vector3<double>(1,1,1)), 0.25);
}

TEST(spatial, amr_side_detJ) {
  ASSERT_EQ(get_amr_side_detJ(1, X, vector3<double>(1,0,0)), 1.);
  ASSERT_EQ(get_amr_side_detJ(2, X, vector3<double>(1,1,0)), 0.25);
  ASSERT_EQ(get_amr_side_detJ(2, Y, vector3<double>(1,1,0)), 0.25);
  ASSERT_EQ(get_amr_side_detJ(3, X, vector3<double>(1,1,1)), 0.0625);
  ASSERT_EQ(get_amr_side_detJ(3, Y, vector3<double>(1,1,1)), 0.0625);
  ASSERT_EQ(get_amr_side_detJ(3, Z, vector3<double>(1,1,1)), 0.0625);
}

TEST(spatial, permute) {
  ASSERT_EQ(permute(X, X), X);
  ASSERT_EQ(permute(X, Y), Y);
  ASSERT_EQ(permute(X, Z), Z);
  ASSERT_EQ(permute(Y, X), Y);
  ASSERT_EQ(permute(Y, Y), Z);
  ASSERT_EQ(permute(Y, Z), X);
  ASSERT_EQ(permute(Z, X), Z);
  ASSERT_EQ(permute(Z, Y), X);
  ASSERT_EQ(permute(Z, Z), Y);
}
