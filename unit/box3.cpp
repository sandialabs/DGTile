#include <dgt_box3.hpp>

#include <gtest/gtest.h>

using namespace dgt;

template <class T>
void test_construct_upper()
{
  Box3<T> const a({T(1),T(1),T(1)});
  EXPECT_EQ(a.lower(), Vec3<T>(T(0),T(0),T(0)));
  EXPECT_EQ(a.upper(), Vec3<T>(T(1),T(1),T(1)));
}

template <class T>
void test_construct_upper_lower()
{
  Vec3<T> const a(T(-1), T(-1), T(-1));
  Vec3<T> const b(T(1), T(1), T(1));
  Box3<T> const c(a, b);
  EXPECT_EQ(c.lower(), a);
  EXPECT_EQ(c.upper(), b);
}

template <class T>
void test_lower_assignment()
{
  Box3<T> a({T(2), T(2), T(2)});
  a.lower().x() = T(1);
  a.lower().y() = T(1);
  a.lower().z() = T(1);
  EXPECT_EQ(a.lower(), Vec3<T>(T(1), T(1), T(1)));
}

template <class T>
void test_upper_assignment()
{
  Box3<T> a({T(2), T(2), T(2)});
  a.upper().x() = T(3);
  a.upper().y() = T(3);
  a.upper().z() = T(3);
  EXPECT_EQ(a.upper(), Vec3<T>(T(3), T(3), T(3)));
}

template <class T>
void test_extents()
{
  Box3<T> const a({T(-1), T(-1), T(-1)}, {T(1), T(1), T(1)});
  EXPECT_EQ(a.extents(), Vec3<T>(T(2), T(2), T(2)));
}

template <class T>
void test_volume()
{
  Box3<T> const a({T(-1), T(-1), T(-1)}, {T(1), T(1), T(1)});
  EXPECT_EQ(a.volume(), T(8));
}

template <class T>
void test_equality()
{
  Box3<T> const a({T(1), T(1), T(1)});
  Box3<T> const b({T(1), T(1), T(1)});
  EXPECT_EQ(a, b);
}

TEST(box3, construct_with_upper)
{
  test_construct_upper<int>();
  test_construct_upper<real>();
}

TEST(box3, construct_with_upper_lower)
{
  test_construct_upper_lower<int>();
  test_construct_upper_lower<real>();
}

TEST(box3, lower_assignment)
{
  test_lower_assignment<int>();
  test_lower_assignment<real>();
}

TEST(box3, upper_assignment)
{
  test_upper_assignment<int>();
  test_upper_assignment<real>();
}

TEST(box3, extents)
{
  test_extents<int>();
  test_extents<real>();
}

TEST(box3, volume)
{
  test_volume<int>();
  test_volume<real>();
}

TEST(box3, equality)
{
  test_equality<int>();
  test_equality<real>();
}
