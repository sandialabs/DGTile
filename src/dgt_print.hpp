#pragma once

#include <iostream>

#include "dgt_grid3.hpp"
#include "dgt_subgrid3.hpp"
#include "dgt_tree.hpp"
#include "dgt_vec.hpp"

namespace dgt {

template <class T>
inline std::ostream& operator<<(std::ostream& os, Vec3<T> const& v)
{
  os << "[" << v.x() << ", " << v.y() << ", " << v.z() << "]";
  return os;
}

template <class T>
inline std::ostream& operator<<(std::ostream& os, Box3<T> const& b)
{
  os << b.lower() << "<->" << b.upper();
  return os;
}

inline std::ostream& operator<<(std::ostream& os, Grid3 const& g)
{
  os << g.extents();
  return os;
}

inline std::ostream& operator<<(std::ostream& os, Subgrid3 const& s)
{
  os << s.lower() << "<->" << s.upper();
  return os;
}

template <class T, int N>
inline std::ostream& operator<<(std::ostream& os, Vec<T, N> const& v)
{
  os << "[";
  for (int i = 0; i < N-1; ++i) {
    os << v[i] << ", ";
  }
  os << v[N-1] << "]";
  return os;
}

inline std::ostream& operator<<(std::ostream& os, tree::Point const& pt) {
  os << "[" << pt.level << "], " << pt.ijk;
  return os;
}

}
