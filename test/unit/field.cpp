#include <dgt_field.hpp>

#include <gtest/gtest.h>

using namespace dgt;

TEST(field, create_real1)
{
  Field<real*> f;
  f.create("f", 2, 10);
  EXPECT_EQ(f.name(), "f");
  EXPECT_EQ(f.get().size(), 2);
  for (int i = 0; i < 2; ++i) {
    EXPECT_EQ(f.get_view(i).label(), fmt::format("f[{}]", i));
    EXPECT_EQ(f.get_view(i).extent(0), 10);
    EXPECT_EQ(f.get()[i].extent(0), 10);
  }
}

TEST(field, create_real2)
{
  Field<real**> f;
  f.create("f", 2, 10, 5);
  EXPECT_EQ(f.name(), "f");
  EXPECT_EQ(f.get().size(), 2);
  for (int i = 0; i < 2; ++i) {
    EXPECT_EQ(f.get_view(i).label(), fmt::format("f[{}]", i));
    EXPECT_EQ(f.get_view(i).extent(0), 10);
    EXPECT_EQ(f.get_view(i).extent(1), 5);
    EXPECT_EQ(f.get()[i].extent(0), 10);
    EXPECT_EQ(f.get()[i].extent(1), 5);
  }
}

TEST(field, create_real3)
{
  Field<real***> f;
  f.create("f", 2, 10, 5, 8);
  EXPECT_EQ(f.name(), "f");
  EXPECT_EQ(f.get().size(), 2);
  for (int i = 0; i < 2; ++i) {
    EXPECT_EQ(f.get_view(i).label(), fmt::format("f[{}]", i));
    EXPECT_EQ(f.get_view(i).extent(0), 10);
    EXPECT_EQ(f.get_view(i).extent(1), 5);
    EXPECT_EQ(f.get_view(i).extent(2), 8);
    EXPECT_EQ(f.get()[i].extent(0), 10);
    EXPECT_EQ(f.get()[i].extent(1), 5);
    EXPECT_EQ(f.get()[i].extent(2), 8);
  }
}

TEST(field, create_no_blocks)
{
  Field<real***> f;
  f.create("f", 0, 10, 5, 8);
  EXPECT_EQ(f.name(), "f");
  EXPECT_EQ(f.get().size(), 0);
}
