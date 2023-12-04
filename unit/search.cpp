#include <dgt_search.hpp>

#include <gtest/gtest.h>

using namespace dgt;

TEST(search, binary_search)
{
  std::vector<double> a = {1., 2.12, 4.87, 7.5, 7.6, 278.1, 401.2};
  EXPECT_EQ(binary_search(0., a), 0);
  EXPECT_EQ(binary_search(1.5, a), 0);
  EXPECT_EQ(binary_search(3.2111, a), 1);
  EXPECT_EQ(binary_search(7.1, a), 2);
  EXPECT_EQ(binary_search(7.5, a), 3);
  EXPECT_EQ(binary_search(7.555, a), 3);
  EXPECT_EQ(binary_search(101.555, a), 4);
  EXPECT_EQ(binary_search(281.1, a), 5);
  EXPECT_EQ(binary_search(450., a), 5);
}
