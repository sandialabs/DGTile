#include <dgt_partitioning.hpp>

#include <gtest/gtest.h>

using namespace dgt;

TEST(partitioning, linear_get_num_local)
{
  using namespace dgt::linear_partitioning;
  EXPECT_EQ(get_num_local(3,2,0), 2);
  EXPECT_EQ(get_num_local(3,2,1), 1);
  for (int i = 0; i < 4; ++i) {
    EXPECT_EQ(get_num_local(16,4,i), 4);
  }
  EXPECT_EQ(get_num_local(9,4,0), 3);
  EXPECT_EQ(get_num_local(9,4,1), 2);
  EXPECT_EQ(get_num_local(9,4,2), 2);
  EXPECT_EQ(get_num_local(9,4,3), 2);
}

TEST(partitioning, linear_get_local_offset)
{
  using namespace dgt::linear_partitioning;
  EXPECT_EQ(get_local_offset(3,2,0), 0);
  EXPECT_EQ(get_local_offset(3,2,1), 2);
  EXPECT_EQ(get_local_offset(9,4,0), 0);
  EXPECT_EQ(get_local_offset(9,4,1), 3);
  EXPECT_EQ(get_local_offset(9,4,2), 5);
  EXPECT_EQ(get_local_offset(9,4,3), 7);
}
