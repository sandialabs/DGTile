#pragma once

#include <algorithm>

#include "dgt_defines.hpp"
#include "dgt_macros.hpp"

namespace dgt {

template <class T>
class Vec3
{

  private:

    T m_x;
    T m_y;
    T m_z;

  public:

    using reference = T&;
    using const_reference = T const&;

    DGT_METHOD Vec3() = default;

    DGT_VOID_METHOD constexpr Vec3(T const& a, T const& b, T const& c)
      :m_x(a)
      ,m_y(b)
      ,m_z(c)
    {
    }

    DGT_VOID_METHOD constexpr Vec3(Vec3<T> const& other)
      :m_x(other.x())
      ,m_y(other.y())
      ,m_z(other.z())
    {
    }

    DGT_METHOD constexpr reference x() { return m_x; }
    DGT_METHOD constexpr reference y() { return m_y; }
    DGT_METHOD constexpr reference z() { return m_z; }

    DGT_METHOD constexpr const_reference x() const { return m_x; }
    DGT_METHOD constexpr const_reference y() const { return m_y; }
    DGT_METHOD constexpr const_reference z() const { return m_z; }

    DGT_METHOD constexpr reference operator[](int const axis)
    {
      if (axis == X) return m_x;
      if (axis == Y) return m_y;
      return m_z;
    }

    DGT_METHOD constexpr const_reference operator[](int const axis) const
    {
      if (axis == X) return m_x;
      if (axis == Y) return m_y;
      return m_z;
    }


    DGT_METHOD static constexpr Vec3 zero()
    { 
      return Vec3<T>(T(0), T(0), T(0));
    }

    DGT_METHOD static constexpr Vec3 ones()
    {
      Vec3<T> result(T(1), T(1), T(1));
      return result;
    }

    DGT_METHOD static constexpr Vec3 axis(int const axis)
    {
      if (axis == X) return Vec3(T(1), T(0), T(0));
      if (axis == Y) return Vec3(T(0), T(1), T(0));
      return Vec3(T(0), T(0), T(1));
    }

    DGT_VOID_METHOD constexpr Vec3& operator=(Vec3 const& other)
    {
      m_x = other.x();
      m_y = other.y();
      m_z = other.z();
      return *this;
    }

    DGT_METHOD constexpr bool operator==(Vec3 const& other) const
    {
      return 
        (m_x == other.m_x) &&
        (m_y == other.m_y) &&
        (m_z == other.m_z);
    }

    DGT_METHOD constexpr bool operator!=(Vec3 const& other) const
    {
      return
        (m_x != other.m_x) ||
        (m_y != other.m_y) ||
        (m_z != other.m_z);
    }

    DGT_VOID_METHOD constexpr void operator+=(Vec3 const& other)
    {
      m_x += other.m_x;
      m_y += other.m_y;
      m_z += other.m_z;
    }

    DGT_VOID_METHOD constexpr void operator-=(Vec3 const& other)
    {
      m_x -= other.m_x;
      m_y -= other.m_y;
      m_z -= other.m_z;
    }

    DGT_VOID_METHOD constexpr void operator*=(T const& scalar)
    {
      m_x *= scalar;
      m_y *= scalar;
      m_z *= scalar;
    }

    DGT_VOID_METHOD constexpr void operator/=(T const& scalar)
    {
      m_x /= scalar;
      m_y /= scalar;
      m_z /= scalar;
    }

    DGT_METHOD constexpr Vec3 operator+(Vec3 const& other) const
    {
      return Vec3(m_x + other.m_x, m_y + other.m_y, m_z + other.m_z);
    }

    DGT_METHOD constexpr Vec3 operator-(Vec3 const& other) const
    {
      return Vec3(m_x - other.m_x, m_y - other.m_y, m_z - other.m_z);
    }

    DGT_METHOD constexpr Vec3 operator/(T const& scalar) const
    {
      return Vec3(m_x / scalar, m_y / scalar, m_z / scalar);
    }

    DGT_METHOD constexpr Vec3 operator*(T const& scalar) const
    {
      return Vec3(m_x * scalar, m_y * scalar, m_z * scalar);
    }

    DGT_METHOD constexpr T volume() const
    {
      return m_x * m_y * m_z;
    }

};

template <class T>
DGT_METHOD constexpr Vec3<T> operator-(Vec3<T> const& v)
{
  return Vec3<T>(-v.x(), -v.y(), -v.z());
}

template <class T>
DGT_METHOD constexpr Vec3<T> operator*(T const& scalar, Vec3<T> const& v)
{
  return v * scalar;
}

template <class T>
DGT_METHOD constexpr Vec3<T> comp_product(Vec3<T> const& a, Vec3<T> const& b)
{
  return Vec3<T>(a.x() * b.x(), a.y() * b.y(), a.z() * b.z());
}

template <class T>
DGT_METHOD constexpr Vec3<T> comp_division(Vec3<T> const& a, Vec3<T> const& b)
{
  return Vec3<T>(a.x() / b.x(), a.y() / b.y(), a.z() / b.z());
}

template <class T>
DGT_METHOD constexpr T dot(Vec3<T> const& a, Vec3<T> const& b)
{
  return a.x() * b.x() + a.y() * b.y() + a.z() * b.z();
}

template <class T>
DGT_METHOD constexpr Vec3<T> abs(Vec3<T> const& v)
{
  return Vec3<T>(std::abs(v.x()), std::abs(v.y()), std::abs(v.z()));
}

template <class T>
DGT_METHOD constexpr T min(Vec3<T> const& v)
{
  return std::min(v.x(), std::min(v.y(), v.z()));
}

template <class T>
DGT_METHOD constexpr T max(Vec3<T> const& v)
{
  return std::max(v.x(), std::max(v.y(), v.z()));
}

DGT_METHOD constexpr int infer_dimension(Vec3<int> const& v)
{
  bool const has_x = (v.x() > 0);
  bool const has_y = (v.y() > 0);
  bool const has_z = (v.z() > 0);
  if (has_x && has_y && has_z) return 3;
  if (has_x && has_y && !has_z) return 2;
  if (has_x && !has_y && !has_z) return 1;
  return -1;
}

template <class T>
DGT_METHOD constexpr Vec3<T> dimensionalize(int const dim, Vec3<T> const& v)
{
  Vec3<T> d_v = v;
  for (int axis = dim; axis < DIMENSIONS; ++axis) {
    d_v[axis] = T(0);
  }
  return d_v;
}

template <class T>
DGT_METHOD constexpr Vec3<T> generalize(int const dim, Vec3<T> const v)
{
  Vec3<T> g_v = v;
  for (int axis = dim; axis < DIMENSIONS; ++axis) {
    g_v[axis] = T(1);
  }
  return g_v;
}

}
