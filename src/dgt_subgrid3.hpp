#pragma once

#include "dgt_box3.hpp"
#include "dgt_grid3.hpp"

namespace dgt {

class Subgrid3
{

  private:

    Box3<int> m_box;

  public:

    DGT_ALWAYS_INLINE Subgrid3() = default;

    DGT_METHOD constexpr Subgrid3(Vec3<int> const& lower, Vec3<int> const& upper)
      :m_box(lower, upper)
    {
    }

    DGT_METHOD constexpr Subgrid3(Grid3 const& grid)
      :Subgrid3(Vec3<int>::zero(), grid.extents())
    {
    }

    DGT_METHOD constexpr Vec3<int>& lower()
    {
      return m_box.lower();
    }

    DGT_METHOD constexpr Vec3<int>& upper()
    {
      return m_box.upper();
    }

    DGT_METHOD constexpr Vec3<int> const& lower() const
    {
      return m_box.lower();
    }

    DGT_METHOD constexpr Vec3<int> const& upper() const
    {
      return m_box.upper();
    }

    DGT_METHOD constexpr Vec3<int> extents() const
    {
      return m_box.upper() - m_box.lower();
    }

    DGT_METHOD constexpr int size() const
    {
      if (m_box.upper().volume() == -1) return -1;
      return extents().volume();
    }

    DGT_METHOD constexpr int index(Vec3<int> const& ijk) const
    {
      return Grid3(extents()).index(ijk - lower());
    }

    DGT_METHOD constexpr bool operator==(Subgrid3 const& other) const
    {
      return m_box == other.m_box;
    }

    DGT_METHOD constexpr bool operator!=(Subgrid3 const& other) const
    {
      return !operator==(other);
    }

};


}
