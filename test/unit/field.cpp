#include <dgt_field.hpp>

#include <gtest/gtest.h>

using namespace dgt;

TEST(field, construct)
{
  Field f("hydro_0", 4, 20, 5, 8);
  EXPECT_EQ(f.name(), "hydro_0");
}

TEST(field, get)
{
  Field f("hydro_0", 4, 20, 5, 8);
  Field::accessor_t hydro = f.get();
  EXPECT_EQ(hydro.size(), 4);
}

TEST(field, get_by_index)
{
  Field f("hydro_0", 4, 20, 5, 8);
  for (int block = 0; block < 4; ++block) {
    Field::view_t hydro = f.get(block);
    EXPECT_EQ(hydro.extent(0), 20);
    EXPECT_EQ(hydro.extent(1), 5);
    EXPECT_EQ(hydro.extent(2), 8);
  }
}

TEST(field, kinds)
{
  EXPECT_EQ(Field::MODAL, 0);
  EXPECT_EQ(Field::FLUX, 1);
  EXPECT_EQ(Field::RESIDUAL, 2);
  EXPECT_EQ(Field::CELLS, 3);
  EXPECT_EQ(Field::FACES, 4);
  EXPECT_EQ(Field::CELL_POINTS, 5);
  EXPECT_EQ(Field::FACE_POINTS, 6);
  EXPECT_EQ(Field::KINDS, 7);
}
