#pragma once

#include <algorithm>

#include "dgt_defines.hpp"
#include "dgt_macros.hpp"

namespace dgt {

template <class T>
class vec3
{

  private:

    T m_x;
    T m_y;
    T m_z;

  public:

    using reference = T&;
    using const_reference = T const&;

    DGT_METHOD vec3() = default;

    DGT_VOID_METHOD constexpr vec3(T const& a, T const& b, T const& c)
      :m_x(a)
      ,m_y(b)
      ,m_z(c)
    {
    }

    DGT_VOID_METHOD constexpr vec3(vec3<T> const& other)
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

    DGT_METHOD static constexpr vec3 zero()
    { 
      return vec3<T>(T(0), T(0), T(0));
    }

    DGT_METHOD static constexpr vec3 ones()
    {
      vec3<T> result(T(1), T(1), T(1));
      return result;
    }

    DGT_METHOD static constexpr vec3 axis(int const axis)
    {
      if (axis == X) return vec3(T(1), T(0), T(0));
      if (axis == Y) return vec3(T(0), T(1), T(0));
      return vec3(T(0), T(0), T(1));
    }

    DGT_VOID_METHOD constexpr vec3& operator=(vec3 const& other)
    {
      m_x = other.x();
      m_y = other.y();
      m_z = other.z();
      return *this;
    }

    DGT_METHOD constexpr bool operator==(vec3 const& other) const
    {
      return 
        (m_x == other.m_x) &&
        (m_y == other.m_y) &&
        (m_z == other.m_z);
    }

    DGT_METHOD constexpr bool operator!=(vec3 const& other) const
    {
      return
        (m_x != other.m_x) ||
        (m_y != other.m_y) ||
        (m_z != other.m_z);
    }

    DGT_VOID_METHOD constexpr void operator+=(vec3 const& other)
    {
      m_x += other.m_x;
      m_y += other.m_y;
      m_z += other.m_z;
    }

    DGT_VOID_METHOD constexpr void operator-=(vec3 const& other)
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

    DGT_METHOD constexpr vec3 operator+(vec3 const& other) const
    {
      return vec3(m_x + other.m_x, m_y + other.m_y, m_z + other.m_z);
    }

    DGT_METHOD constexpr vec3 operator-(vec3 const& other) const
    {
      return vec3(m_x - other.m_x, m_y - other.m_y, m_z - other.m_z);
    }

    DGT_METHOD constexpr vec3 operator/(T const& scalar) const
    {
      return vec3(m_x / scalar, m_y / scalar, m_z / scalar);
    }

    DGT_METHOD constexpr vec3 operator*(T const& scalar) const
    {
      return vec3(m_x * scalar, m_y * scalar, m_z * scalar);
    }

    DGT_METHOD constexpr T size() const
    {
      return m_x * m_y * m_z;
    }

};

template <class T>
DGT_METHOD constexpr vec3<T> operator-(vec3<T> const& v)
{
  return vec3<T>(-v.x(), -v.y(), -v.z());
}

template <class T>
DGT_METHOD constexpr vec3<T> operator*(T const& scalar, vec3<T> const& v)
{
  return v * scalar;
}

template <class T>
DGT_METHOD constexpr vec3<T> comp_product(vec3<T> const& a, vec3<T> const& b)
{
  return vec3<T>(a.x() * b.x(), a.y() * b.y(), a.z() * b.z());
}

template <class T>
DGT_METHOD constexpr vec3<T> comp_division(vec3<T> const& a, vec3<T> const& b)
{
  return vec3<T>(a.x() / b.x(), a.y() / b.y(), a.z() / b.z());
}

template <class T>
DGT_METHOD constexpr vec3<T> abs(vec3<T> const& v)
{
  return vec3<T>(std::abs(v.x()), std::abs(v.y()), std::abs(v.z()));
}

template <class T>
DGT_METHOD constexpr T min(vec3<T> const& v)
{
  return std::min(v.x(), std::min(v.y(), v.z()));
}

template <class T>
DGT_METHOD constexpr T max(vec3<T> const& v)
{
  return std::max(v.x(), std::max(v.y(), v.z()));
}

}
