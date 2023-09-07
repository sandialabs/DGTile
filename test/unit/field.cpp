#include <dgt_field.hpp>

#include <gtest/gtest.h>

using namespace dgt;

TEST(field, construct)
{
  int const num_blocks = 4;
  int const num_cells = 20;
  int const num_eqs = 5;
  BasisInfo basis = {3, 1, 2, true};
  DGField f("hydro_0", num_blocks, num_cells, num_eqs, basis);
  EXPECT_EQ(f.name(), "hydro_0");
}

TEST(field, get)
{
  int const num_blocks = 4;
  int const num_cells = 20;
  int const num_eqs = 5;
  BasisInfo basis = {3, 1, 2, true};
  DGField f("hydro_0", num_blocks, num_cells, num_eqs, basis);
  DGField::accessor_t hydro = f.get();
  EXPECT_EQ(hydro.size(), num_blocks);
}

TEST(field, get_by_index)
{
  int const num_blocks = 4;
  int const num_cells = 20;
  int const num_eqs = 5;
  int const num_modes = 8;
  BasisInfo basis = {3, 1, 2, true};
  DGField f("hydro_0", num_blocks, num_cells, num_eqs, basis);
  for (int block = 0; block < num_blocks; ++block) {
    DGField::view_t hydro_block = f.get(block);
    EXPECT_EQ(hydro_block.extent(0), num_cells);
    EXPECT_EQ(hydro_block.extent(1), num_eqs);
    EXPECT_EQ(hydro_block.extent(2), num_modes);
  }
}
