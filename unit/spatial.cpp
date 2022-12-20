#include "gtest/gtest.h"

#include "dgt_spatial.hpp"

TEST(spatial, dir_sign) {
  ASSERT_EQ(dgt::get_dir_sign(dgt::left), -1);
  ASSERT_EQ(dgt::get_dir_sign(dgt::right), 1);
}

TEST(spatial, invert_dir) {
  ASSERT_EQ(dgt::invert_dir(dgt::left), dgt::right);
  ASSERT_EQ(dgt::invert_dir(dgt::right), dgt::left);
}

TEST(spatial, sides_adj_cell) {
  ASSERT_EQ(dgt::get_sides_adj_cell({4,4,4}, dgt::X, dgt::left), p3a::vector3<int>(3,4,4));
  ASSERT_EQ(dgt::get_sides_adj_cell({4,4,4}, dgt::X, dgt::right), p3a::vector3<int>(4,4,4));
  ASSERT_EQ(dgt::get_sides_adj_cell({4,4,4}, dgt::Y, dgt::left), p3a::vector3<int>(4,3,4));
  ASSERT_EQ(dgt::get_sides_adj_cell({4,4,4}, dgt::Y, dgt::right), p3a::vector3<int>(4,4,4));
  ASSERT_EQ(dgt::get_sides_adj_cell({4,4,4}, dgt::Z, dgt::left), p3a::vector3<int>(4,4,3));
  ASSERT_EQ(dgt::get_sides_adj_cell({4,4,4}, dgt::Z, dgt::right), p3a::vector3<int>(4,4,4));
}

TEST(spatial, cells_adj_side) {
  ASSERT_EQ(dgt::get_cells_adj_side({4,4,4}, dgt::X, dgt::left), p3a::vector3<int>(4,4,4));
  ASSERT_EQ(dgt::get_cells_adj_side({4,4,4}, dgt::X, dgt::right), p3a::vector3<int>(5,4,4));
  ASSERT_EQ(dgt::get_cells_adj_side({4,4,4}, dgt::Y, dgt::left), p3a::vector3<int>(4,4,4));
  ASSERT_EQ(dgt::get_cells_adj_side({4,4,4}, dgt::Y, dgt::right), p3a::vector3<int>(4,5,4));
  ASSERT_EQ(dgt::get_cells_adj_side({4,4,4}, dgt::Z, dgt::left), p3a::vector3<int>(4,4,4));
  ASSERT_EQ(dgt::get_cells_adj_side({4,4,4}, dgt::Z, dgt::right), p3a::vector3<int>(4,4,5));
}

TEST(spatial, cells_adj_cell) {
  ASSERT_EQ(dgt::get_cells_adj_cell({4,4,4}, dgt::X, dgt::left), p3a::vector3<int>(3,4,4));
  ASSERT_EQ(dgt::get_cells_adj_cell({4,4,4}, dgt::X, dgt::right), p3a::vector3<int>(5,4,4));
  ASSERT_EQ(dgt::get_cells_adj_cell({4,4,4}, dgt::Y, dgt::left), p3a::vector3<int>(4,3,4));
  ASSERT_EQ(dgt::get_cells_adj_cell({4,4,4}, dgt::Y, dgt::right), p3a::vector3<int>(4,5,4));
  ASSERT_EQ(dgt::get_cells_adj_cell({4,4,4}, dgt::Z, dgt::left), p3a::vector3<int>(4,4,3));
  ASSERT_EQ(dgt::get_cells_adj_cell({4,4,4}, dgt::Z, dgt::right), p3a::vector3<int>(4,4,5));
}

TEST(spatial, border_ijk) {
  ASSERT_EQ(dgt::get_border_ijk({4,4,4}, dgt::X), p3a::vector3<int>(0,4,4));
  ASSERT_EQ(dgt::get_border_ijk({4,4,4}, dgt::Y), p3a::vector3<int>(4,0,4));
  ASSERT_EQ(dgt::get_border_ijk({4,4,4}, dgt::Z), p3a::vector3<int>(4,4,0));
}

TEST(spatial, dx) {
  ASSERT_EQ(dgt::get_dx({{0.,0.,0.},{1.,0.,0.}}, {10,0,0}), p3a::vector3<double>(.1,0.,0.));
  ASSERT_EQ(dgt::get_dx({{0.,0.,0.},{1.,1.,0.}}, {10,10,0}), p3a::vector3<double>(.1,.1,0.));
  ASSERT_EQ(dgt::get_dx({{0.,0.,0.},{1.,1.,1.}}, {10,10,10}), p3a::vector3<double>(.1,.1,.1));
}

TEST(spatial, block_domain_1D) {
  dgt::Point const base{1, {2,0,0}};
  dgt::Point const node{2, {2,0,0}};
  p3a::box3<double> const domain(p3a::vector3<double>(0,0,0), p3a::vector3<double>(1,0,0));
  p3a::box3<double> const exact(p3a::vector3<double>(.5,0,0), p3a::vector3<double>(.75,0,0));
  ASSERT_EQ(dgt::get_block_domain(base, node, domain), exact);
}

TEST(spatial, block_domain_2D) {
  dgt::Point const base{1, {2,2,0}};
  dgt::Point const node{2, {2,2,0}};
  p3a::box3<double> const domain(p3a::vector3<double>(0,0,0), p3a::vector3<double>(1,1,0));
  p3a::box3<double> const exact(p3a::vector3<double>(.5,.5,0), p3a::vector3<double>(.75,.75,0));
  ASSERT_EQ(dgt::get_block_domain(base, node, domain), exact);
}

TEST(spatial, block_domain_3D) {
  dgt::Point const base{1, {2,2,2}};
  dgt::Point const node{2, {2,2,2}};
  p3a::box3<double> const domain(p3a::vector3<double>(0,0,0), p3a::vector3<double>(1,1,1));
  p3a::box3<double> const exact(p3a::vector3<double>(.5,.5,.5), p3a::vector3<double>(.75,.75,.75));
  ASSERT_EQ(dgt::get_block_domain(base, node, domain), exact);
}

TEST(spatial, cell_center) {
  ASSERT_EQ(dgt::get_cell_center({2,0,0}, p3a::vector3<double>(0,0,0), p3a::vector3<double>(1,0,0)), p3a::vector3<double>(2.5,0,0));
  ASSERT_EQ(dgt::get_cell_center({2,2,0}, p3a::vector3<double>(0,0,0), p3a::vector3<double>(1,1,0)), p3a::vector3<double>(2.5,2.5,0));
  ASSERT_EQ(dgt::get_cell_center({2,2,2}, p3a::vector3<double>(0,0,0), p3a::vector3<double>(1,1,1)), p3a::vector3<double>(2.5,2.5,2.5));
}

TEST(spatial, x) {
  ASSERT_EQ(dgt::get_x({2,0,0}, p3a::vector3<double>(0,0,0), p3a::vector3<double>(1,0,0), p3a::vector3<double>(-1,0,0)), p3a::vector3<double>(2,0,0));
  ASSERT_EQ(dgt::get_x({2,0,0}, p3a::vector3<double>(0,0,0), p3a::vector3<double>(1,0,0), p3a::vector3<double>(1,0,0)), p3a::vector3<double>(3,0,0));
  ASSERT_EQ(dgt::get_x({2,2,0}, p3a::vector3<double>(0,0,0), p3a::vector3<double>(1,1,0), p3a::vector3<double>(-1,-1,0)), p3a::vector3<double>(2,2,0));
  ASSERT_EQ(dgt::get_x({2,2,0}, p3a::vector3<double>(0,0,0), p3a::vector3<double>(1,1,0), p3a::vector3<double>(1,1,0)), p3a::vector3<double>(3,3,0));
  ASSERT_EQ(dgt::get_x({2,2,2}, p3a::vector3<double>(0,0,0), p3a::vector3<double>(1,1,1), p3a::vector3<double>(-1,-1,-1)), p3a::vector3<double>(2,2,2));
  ASSERT_EQ(dgt::get_x({2,2,2}, p3a::vector3<double>(0,0,0), p3a::vector3<double>(1,1,1), p3a::vector3<double>(1,1,1)), p3a::vector3<double>(3,3,3));
}

TEST(spatial, volume) {
  ASSERT_EQ(dgt::get_volume(1, p3a::vector3<double>(10,0,0)), 10.);
  ASSERT_EQ(dgt::get_volume(2, p3a::vector3<double>(10,10,0)), 100.);
  ASSERT_EQ(dgt::get_volume(3, p3a::vector3<double>(10,10,10)), 1000.);
}

TEST(spatial, cell_detJ) {
  ASSERT_EQ(dgt::get_cell_detJ(1, p3a::vector3<double>(1,0,0)), 0.5);
  ASSERT_EQ(dgt::get_cell_detJ(2, p3a::vector3<double>(1,1,0)), 0.25);
  ASSERT_EQ(dgt::get_cell_detJ(3, p3a::vector3<double>(1,1,1)), 0.125);
}

TEST(spatial, side_detJ) {
  ASSERT_EQ(dgt::get_side_detJ(1, dgt::X, p3a::vector3<double>(1,0,0)), 1.);
  ASSERT_EQ(dgt::get_side_detJ(2, dgt::X, p3a::vector3<double>(1,1,0)), 0.5);
  ASSERT_EQ(dgt::get_side_detJ(2, dgt::Y, p3a::vector3<double>(1,1,0)), 0.5);
  ASSERT_EQ(dgt::get_side_detJ(3, dgt::X, p3a::vector3<double>(1,1,1)), 0.25);
  ASSERT_EQ(dgt::get_side_detJ(3, dgt::Y, p3a::vector3<double>(1,1,1)), 0.25);
  ASSERT_EQ(dgt::get_side_detJ(3, dgt::Z, p3a::vector3<double>(1,1,1)), 0.25);
}

TEST(spatial, amr_side_detJ) {
  ASSERT_EQ(dgt::get_amr_side_detJ(1, dgt::X, p3a::vector3<double>(1,0,0)), 1.);
  ASSERT_EQ(dgt::get_amr_side_detJ(2, dgt::X, p3a::vector3<double>(1,1,0)), 0.25);
  ASSERT_EQ(dgt::get_amr_side_detJ(2, dgt::Y, p3a::vector3<double>(1,1,0)), 0.25);
  ASSERT_EQ(dgt::get_amr_side_detJ(3, dgt::X, p3a::vector3<double>(1,1,1)), 0.0625);
  ASSERT_EQ(dgt::get_amr_side_detJ(3, dgt::Y, p3a::vector3<double>(1,1,1)), 0.0625);
  ASSERT_EQ(dgt::get_amr_side_detJ(3, dgt::Z, p3a::vector3<double>(1,1,1)), 0.0625);
}

TEST(spatial, permute) {
  ASSERT_EQ(dgt::permute(dgt::X, dgt::X), dgt::X);
  ASSERT_EQ(dgt::permute(dgt::X, dgt::Y), dgt::Y);
  ASSERT_EQ(dgt::permute(dgt::X, dgt::Z), dgt::Z);
  ASSERT_EQ(dgt::permute(dgt::Y, dgt::X), dgt::Y);
  ASSERT_EQ(dgt::permute(dgt::Y, dgt::Y), dgt::Z);
  ASSERT_EQ(dgt::permute(dgt::Y, dgt::Z), dgt::X);
  ASSERT_EQ(dgt::permute(dgt::Z, dgt::X), dgt::Z);
  ASSERT_EQ(dgt::permute(dgt::Z, dgt::Y), dgt::X);
  ASSERT_EQ(dgt::permute(dgt::Z, dgt::Z), dgt::Y);
}
