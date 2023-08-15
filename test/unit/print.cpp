#include <dgt_print.hpp>

#include <gtest/gtest.h>

using namespace dgt;

TEST(print, print_vec3_int)
{
  Vec3<int> v(1,2,3);
  std::cout << v << "\n";
}

TEST(print, print_vec3_real)
{
  Vec3<int> v(1.,2.,3.);
  std::cout << v << "\n";
}
