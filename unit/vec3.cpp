#include <dgt_vec3.hpp>

#include <gtest/gtest.h>

using namespace dgt;

template <class T>
void test_construct()
{
  Vec3<T> a(T(1), T(1), T(1));
  EXPECT_EQ(a, Vec3<T>::ones());
}

template <class T>
void test_construct_list()
{
  Vec3<T> a = {T(1), T(1), T(1)};
  EXPECT_EQ(a, Vec3<T>::ones());
}

template <class T>
void test_construct_copy()
{
  Vec3<T> const a = {T(1), T(1), T(1)};
  Vec3<T> const b(a);
  EXPECT_EQ(a, b);
}

template <class T>
void test_zero()
{
  Vec3<T> const a = Vec3<T>::zero();
  EXPECT_EQ(a.x(), T(0));
  EXPECT_EQ(a.y(), T(0));
  EXPECT_EQ(a.z(), T(0));
}

template <class T>
void test_ones()
{
  Vec3<T> const a = Vec3<T>::ones();
  EXPECT_EQ(a.x(), T(1));
  EXPECT_EQ(a.y(), T(1));
  EXPECT_EQ(a.z(), T(1));
}

template <class T>
void test_axis()
{
  Vec3<T> const x = Vec3<T>::axis(X);
  Vec3<T> const y = Vec3<T>::axis(Y);
  Vec3<T> const z = Vec3<T>::axis(Z);
  EXPECT_EQ(x, Vec3<T>(T(1), T(0), T(0)));
  EXPECT_EQ(y, Vec3<T>(T(0), T(1), T(0)));
  EXPECT_EQ(z, Vec3<T>(T(0), T(0), T(1)));
}

template <class T>
void test_assignment()
{
  Vec3<T> const a = Vec3<T>::ones();
  Vec3<T> const b = a;
  EXPECT_EQ(a, b);
}

template <class T>
void test_unary_addition()
{
  Vec3<T> a = Vec3<T>::ones();
  a += Vec3<T>::zero();
  EXPECT_EQ(a, Vec3<T>::ones());
}

template <class T>
void test_unary_subtraction()
{
  Vec3<T> a = Vec3<T>::ones();
  a -= Vec3<T>::zero();
  EXPECT_EQ(a, Vec3<T>::ones());
}

template <class T>
void test_unary_multiplication()
{
  Vec3<T> a = Vec3<T>::ones();
  a *= T(1);
  EXPECT_EQ(a, Vec3<T>::ones());
}

template <class T>
void test_unary_division()
{
  Vec3<T> a = Vec3<T>::ones();
  a /= T(1);
  EXPECT_EQ(a, Vec3<T>::ones());
}

template <class T>
void test_binary_addition()
{
  Vec3<T> const a = Vec3<T>::ones();
  Vec3<T> const b = Vec3<T>::ones();
  Vec3<T> const c = a + b;
  EXPECT_EQ(c.x(), T(2));
  EXPECT_EQ(c.y(), T(2));
  EXPECT_EQ(c.z(), T(2));
}

template <class T>
void test_binary_subtraction()
{
  Vec3<T> const a = Vec3<T>::ones();
  Vec3<T> const b = Vec3<T>::ones();
  Vec3<T> const c = a - b;
  EXPECT_EQ(c, Vec3<T>::zero());
}

template <class T>
void test_binary_multiplication()
{
  Vec3<T> const a = Vec3<T>::ones();
  T const b = T(1);
  Vec3<T> const c = a * b;
  EXPECT_EQ(c, Vec3<T>::ones());
}

template <class T>
void test_binary_division()
{
  Vec3<T> const a = Vec3<T>::ones();
  T const b = T(1);
  Vec3<T> const c = a / b;
  EXPECT_EQ(c, Vec3<T>::ones());
}

template <class T>
void test_size()
{
  Vec3<T> const a(T(1), T(2), T(3));
  EXPECT_EQ(a.size(), T(6));
}

template <class T>
void test_abs()
{
  Vec3<T> const a = -Vec3<T>::ones();
  Vec3<T> const b = abs(a);
  EXPECT_EQ(b, Vec3<T>::ones());
}

template <class T>
void test_min()
{
  Vec3<T> const a(T(1), T(2), T(3));
  EXPECT_EQ(min(a), T(1));
}

template <class T>
void test_max()
{
  Vec3<T> const a(T(1), T(2), T(3));
  EXPECT_EQ(max(a), T(3));
}

template <class T>
void test_negation()
{
  Vec3<T> const a = Vec3<T>::ones();
  Vec3<T> const b = -a;
  EXPECT_EQ(b.x(), T(-1));
  EXPECT_EQ(b.y(), T(-1));
  EXPECT_EQ(b.z(), T(-1));
}

template <class T>
void test_left_multiply()
{
  T const a = T(1);
  Vec3<T> const b = Vec3<T>::ones();
  Vec3<T> const c = a * b;
  EXPECT_EQ(c, Vec3<T>::ones());
}

template <class T>
void test_comp_product()
{
  Vec3<T> const a = Vec3<T>(T(4), T(4), T(4));
  Vec3<T> const b = Vec3<T>(T(2), T(2), T(2));
  EXPECT_EQ(comp_product(a,b), Vec3<T>(T(8), T(8), T(8)));
}

template <class T>
void test_comp_division()
{
  Vec3<T> const a = Vec3<T>(T(4), T(4), T(4));
  Vec3<T> const b = Vec3<T>(T(2), T(2), T(2));
  EXPECT_EQ(comp_division(a,b), Vec3<T>(T(2), T(2), T(2)));
}

TEST(vec3, construct)
{
  test_construct<int>();
  test_construct<real>();
}

TEST(vec3, construct_list)
{
  test_construct_list<int>();
  test_construct_list<real>();
}

TEST(vec3, construct_copy)
{
  test_construct_copy<int>();
  test_construct_copy<real>();
}

TEST(vec3, zero)
{
  test_zero<int>();
  test_zero<real>();
}

TEST(vec3, ones)
{
  test_ones<int>();
  test_ones<real>();
}

TEST(vec3, axis)
{
  test_axis<int>();
  test_axis<real>();
}

TEST(vec3, assignment)
{
  test_assignment<int>();
  test_assignment<real>();
}

TEST(vec3, unary_addition)
{
  test_unary_addition<int>();
  test_unary_addition<real>();
}

TEST(vec3, unary_subtraction)
{
  test_unary_subtraction<int>();
  test_unary_subtraction<real>();
}

TEST(vec3, unary_multiplication)
{
  test_unary_multiplication<int>();
  test_unary_multiplication<real>();
}

TEST(vec3, unary_division)
{
  test_unary_division<int>();
  test_unary_division<real>();
}

TEST(vec3, binary_addition)
{
  test_binary_addition<int>();
  test_binary_addition<real>();
}

TEST(vec3, binary_subtraction)
{
  test_binary_subtraction<int>();
  test_binary_subtraction<real>();
}

TEST(vec3, binary_multiplication)
{
  test_binary_multiplication<int>();
  test_binary_multiplication<real>();
}

TEST(vec3, binary_division)
{
  test_binary_division<int>();
  test_binary_division<real>();
}

TEST(vec3, size)
{
  test_size<int>();
  test_size<real>();
}

TEST(vec3, negation)
{
  test_negation<int>();
  test_negation<real>();
}

TEST(vec3, left_multiply)
{
  test_left_multiply<int>();
  test_left_multiply<real>();
}

TEST(vec3, comp_product)
{
  test_comp_product<int>();
  test_comp_product<real>();
}

TEST(vec3, comp_division)
{
  test_comp_division<int>();
  test_comp_division<real>();
}

TEST(vec3, abs)
{
  test_abs<int>();
  test_abs<real>();
}

TEST(vec3, min)
{
  test_min<int>();
  test_min<real>();
}

TEST(vec3, max)
{
  test_max<int>();
  test_max<real>();
}
