#pragma once

#include "dgt_defines.hpp"
#include "dgt_macros.hpp"

namespace dgt {

template <class T, int N>
class Vec
{

  private:

    T m_data[N] = {};

  public:

    using reference = T&;
    using const_reference = T const&;

    DGT_ALWAYS_INLINE explicit constexpr Vec() = default;

    DGT_METHOD constexpr reference operator[](int const i)
    {
      return m_data[i];
    }

    DGT_METHOD constexpr const_reference operator[](int const i) const
    {
      return m_data[i];
    }

    DGT_METHOD static Vec zero()
    {
      Vec result;
      for (int i = 0; i < N; ++i) {
        result.m_data[i] = T(0);
      }
      return result;
    }

    DGT_METHOD static Vec ones()
    {
      Vec result;
      for (int i = 0; i < N; ++i) {
        result.m_data[i] = T(1);
      }
      return result;
    }

    DGT_VOID_METHOD constexpr void operator+=(Vec const& other)
    {
      for (int i = 0; i < N; ++i) {
        m_data[i] += other.m_data[i];
      }
    }

    DGT_VOID_METHOD constexpr void operator-=(Vec const& other)
    {
      for (int i = 0; i < N; ++i) {
        m_data[i] -= other.m_data[i];
      }
    }

    DGT_VOID_METHOD constexpr void operator*=(T const& scalar)
    {
      for (int i = 0; i < N; ++i) {
        m_data[i] *= scalar;
      }
    }

    DGT_VOID_METHOD constexpr void operator/=(T const& scalar)
    {
      for (int i = 0; i < N; ++i) {
        m_data[i] /= scalar;
      }
    }

    DGT_METHOD constexpr Vec operator+(Vec const& other)
    {
      Vec result;
      for (int i = 0; i < N; ++i) {
        result.m_data[i] = m_data[i] + other.m_data[i];
      }
      return result;
    }

    DGT_METHOD constexpr Vec operator-(Vec const& other)
    {
      Vec result;
      for (int i = 0; i < N; ++i) {
        result.m_data[i] = m_data[i] - other.m_data[i];
      }
      return result;
    }

    DGT_METHOD constexpr Vec operator*(T const& scalar)
    {
      Vec result;
      for (int i = 0; i < N; ++i) {
        result.m_data[i] = m_data[i] * scalar;
      }
      return result;
    }

    DGT_METHOD constexpr Vec operator/(T const& scalar)
    {
      Vec result;
      for (int i = 0; i < N; ++i) {
        result.m_data[i] = m_data[i] / scalar;
      }
      return result;
    }

};

template <class T, int N>
DGT_METHOD constexpr Vec<T, N> operator-(Vec<T, N> const& v)
{
  Vec<T, N> result;
  for (int i = 0; i < N; ++i) {
    result[i] = -v[i];
  }
  return result;
}

template <class T, int N>
DGT_METHOD constexpr Vec<T, N> operator*(T const& scalar, Vec<T, N> const& v)
{
  return v * scalar;
}

}
