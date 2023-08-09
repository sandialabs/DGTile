#pragma once

#include "dgt_vec3.hpp"

namespace dgt {

class Grid3
{

  private:

    Vec3<int> m_extents;

  public:

  DGT_ALWAYS_INLINE Grid3() = default;

  DGT_VOID_METHOD constexpr Grid3(int const a, int const b, int const c)
    :m_extents(a, b, c)
  {
  }

  DGT_VOID_METHOD constexpr Grid3(Vec3<int> const& extents)
    :m_extents(extents)
  {
  }

  DGT_METHOD constexpr Vec3<int> const& extents() const
  {
    return m_extents;
  }

  DGT_METHOD constexpr int dimension() const
  {
    return m_extents.dimension();
  }

  DGT_METHOD constexpr int size() const
  {
    return m_extents.volume();
  }

  DGT_METHOD constexpr int index(Vec3<int> const& ijk) const
  {
    return ijk.x() + m_extents.x()*(ijk.y() + m_extents.y() * ijk.z());
  }

  DGT_METHOD constexpr Vec3<int> ijk(int const index) const
  {
    int i(0), j(0), k(0);
    int idx = index;
    i = idx % m_extents.x();
    if (dimension() > 1) {
      idx /= m_extents.x();
      j = idx % m_extents.y();
    }
    if (dimension() > 2) {
      idx /= m_extents.y();
      k = idx;
    }
    return Vec3<int>(i, j ,k);
  }

  DGT_METHOD constexpr bool operator==(Grid3 const& other) const
  {
    return m_extents == other.m_extents;
  }

  DGT_METHOD constexpr bool operator!=(Grid3 const& other) const
  {
    return !operator==(other);
  }

};

}
