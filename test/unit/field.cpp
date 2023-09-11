#include <dgt_cartesian.hpp>
#include <dgt_field.hpp>

#include <fmt/core.h>

#include <gtest/gtest.h>

using namespace dgt;

TEST(field, construct_modal)
{
  Grid3 cell_grid(4,4,4);
  int const num_blocks = 4;
  int const num_eqs = 5;
  int const num_modes = 8;
  field::Modal f("modal", cell_grid, num_blocks, num_eqs, num_modes);
  EXPECT_EQ(f.get().size(), num_blocks);
  for (int block = 0; block < num_blocks; ++block) {
    EXPECT_EQ(f.get()[block].extent(0), 64);
    EXPECT_EQ(f.get()[block].extent(1), num_eqs);
    EXPECT_EQ(f.get()[block].extent(2), num_modes);
    EXPECT_EQ(f.get_view(block).extent(0), 64);
    EXPECT_EQ(f.get_view(block).extent(1), num_eqs);
    EXPECT_EQ(f.get_view(block).extent(2), num_modes);
    EXPECT_EQ(f.get_view(block).label(), fmt::format("modal[{}]", block));
  }
}

TEST(field, construct_cell)
{
  Grid3 cell_grid(4,4,4);
  int const num_blocks = 4;
  int const num_eqs = 5;
  field::Cell f("cell", cell_grid, num_blocks, num_eqs);
  EXPECT_EQ(f.get().size(), num_blocks);
  for (int block = 0; block < num_blocks; ++block) {
    EXPECT_EQ(f.get()[block].extent(0), 64);
    EXPECT_EQ(f.get()[block].extent(1), num_eqs);
    EXPECT_EQ(f.get_view(block).extent(0), 64);
    EXPECT_EQ(f.get_view(block).extent(1), num_eqs);
    EXPECT_EQ(f.get_view(block).label(), fmt::format("cell[{}]", block));
  }
}

TEST(field, construct_cell_points)
{
  Grid3 const cell_grid(4,4,4);
  int const num_blocks = 4;
  int const num_points = 8;
  int const num_eqs = 5;
  field::CellPoints f("cell_points", cell_grid, num_blocks, num_points, num_eqs);
  EXPECT_EQ(f.get().size(), num_blocks);
  for (int block = 0; block < num_blocks; ++block) {
    EXPECT_EQ(f.get()[block].extent(0), 64);
    EXPECT_EQ(f.get()[block].extent(1), num_points);
    EXPECT_EQ(f.get()[block].extent(2), num_eqs);
    EXPECT_EQ(f.get_view(block).extent(0), 64);
    EXPECT_EQ(f.get_view(block).extent(1), num_points);
    EXPECT_EQ(f.get_view(block).extent(2), num_eqs);
    EXPECT_EQ(f.get_view(block).label(), fmt::format("cell_points[{}]", block));
  }
}

TEST(field, construct_face)
{
  Grid3 cell_grid(4,4,4);
  int const num_blocks = 4;
  int const num_eqs = 5;
  field::Face f("face", cell_grid, num_blocks, num_eqs);
  for (int axis = 0; axis < DIMENSIONS; ++axis) {
    auto const axis_name = get_axis_name(axis);
    for (int block = 0; block < num_blocks; ++block) {
      EXPECT_EQ(f.get()[axis][block].extent(0), 80);
      EXPECT_EQ(f.get()[axis][block].extent(1), num_eqs);
      EXPECT_EQ(f.get_view(axis, block).extent(0), 80);
      EXPECT_EQ(f.get_view(axis, block).extent(1), num_eqs);
      EXPECT_EQ(f.get_view(axis, block).extent(1), num_eqs);
      EXPECT_EQ(f.get_view(axis, block).label(), fmt::format("face{}[{}]", axis_name, block));
      EXPECT_EQ(f.get_view(axis, block).label(), fmt::format("face{}[{}]", axis_name, block));
      EXPECT_EQ(f.get_view(axis, block).label(), fmt::format("face{}[{}]", axis_name, block));
    }
  }
}

TEST(field, construct_face_points)
{
  Grid3 cell_grid(4,4,4);
  int const num_blocks = 4;
  int const num_points = 4;
  int const num_eqs = 5;
  field::FacePoints f("face_points", cell_grid, num_blocks, num_points, num_eqs);
  for (int axis = 0; axis < DIMENSIONS; ++axis) {
    auto const axis_name = get_axis_name(axis);
    for (int block = 0; block < num_blocks; ++block) {
      EXPECT_EQ(f.get()[axis][block].extent(0), 80);
      EXPECT_EQ(f.get()[axis][block].extent(1), num_points);
      EXPECT_EQ(f.get()[axis][block].extent(2), num_eqs);
      EXPECT_EQ(f.get_view(axis, block).extent(0), 80);
      EXPECT_EQ(f.get_view(axis, block).extent(1), num_points);
      EXPECT_EQ(f.get_view(axis, block).extent(2), num_eqs);
      EXPECT_EQ(f.get_view(axis, block).label(), fmt::format("face_points{}[{}]", axis_name, block));
      EXPECT_EQ(f.get_view(axis, block).label(), fmt::format("face_points{}[{}]", axis_name, block));
      EXPECT_EQ(f.get_view(axis, block).label(), fmt::format("face_points{}[{}]", axis_name, block));
    }
  }
}
