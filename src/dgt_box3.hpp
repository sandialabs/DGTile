#pragma once

#include "dgt_vec3.hpp"

namespace dgt {

template <class T>
class Box3
{

  private:

    Vec3<T> m_lower;
    Vec3<T> m_upper;

  public:

    DGT_ALWAYS_INLINE Box3() = default;

    DGT_METHOD constexpr Box3(Vec3<T> const& upper)
      :m_lower(Vec3<T>::zero())
      ,m_upper(upper)
    {
    }

    DGT_METHOD constexpr Box3(Vec3<T> const& lower, Vec3<T> const& upper)
      :m_lower(lower)
      ,m_upper(upper)
    {
    }

    DGT_METHOD constexpr Vec3<T>& lower()
    {
      return m_lower;
    }

    DGT_METHOD constexpr Vec3<T>& upper()
    {
      return m_upper;
    }

    DGT_METHOD constexpr Vec3<T> const& lower() const
    {
      return m_lower;
    }

    DGT_METHOD constexpr Vec3<T> const& upper() const
    {
      return m_upper;
    }

    DGT_METHOD constexpr Vec3<T> extents() const
    {
      return m_upper - m_lower;
    }

    DGT_METHOD constexpr T volume() const
    {
      return extents().volume();
    }

    DGT_METHOD constexpr bool operator==(Box3 const& other) const
    {
      return
        (m_lower == other.m_lower) &&
        (m_upper == other.m_upper);
    }

    DGT_METHOD constexpr bool operator!=(Box3 const& other) const
    {
      return !operator==(other);
    }

};

}
