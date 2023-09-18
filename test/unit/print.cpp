#include <dgt_print.hpp>

#include <gtest/gtest.h>

using namespace dgt;

template <class T>
void test_print_vec3(Vec3<T> const& v)
{
  std::cout << v << "\n";
}

template <class T>
void test_print_box3(Box3<T> const& b)
{
  std::cout << b << "\n";
}

template <class T, int N>
void test_print_vec(Vec<T, N> const& v)
{
  std::cout << v << "\n";
}

TEST(print, vec3)
{
  test_print_vec3(Vec3<int>(1,1,1));
  test_print_vec3(Vec3<real>(1.2,2.2,1.2));
}

TEST(print, box3)
{
  test_print_box3(Box3<int>({0,0,0},{2,2,2}));
  test_print_box3(Box3<real>({1.2,2.2,1.2},{2.,2.,2.}));
}

TEST(print, grid3)
{
  Grid3 const g(2,2,2);
  std::cout << g << "\n";
}

TEST(print, subgrid3)
{
  Subgrid3 const s({1,1,1}, {5,5,5});
  std::cout << s << "\n";
}

TEST(print, vec)
{
  Vec<int, 2> v1;
  Vec<real, 2> v2;
  v1[0] = 1;
  v1[1] = 2;
  v2[0] = 1.1;
  v2[1] = 2.1;
  test_print_vec(v1);
  test_print_vec(v2);
}

TEST(print, point)
{
  tree::Point const pt(4, {4,5,6});
  std::cout << pt << "\n";
}

TEST(print, write_stream)
{
  std::stringstream a;
  a << "hi\n";
  a << "bye\n";
  write_stream("out_stream", a);
}
