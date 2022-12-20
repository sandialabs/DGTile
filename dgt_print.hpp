#pragma once

#include "p3a_grid3.hpp"
#include "p3a_static_vector.hpp"

#include "dgt_point.hpp"

namespace p3a {

template <class T>
inline std::ostream& operator<<(std::ostream& os, vector3<T> const& v) {
  os << "[" << v.x() << ", " << v.y() << ", " << v.z() << "]";
  return os;
}

template <class T>
inline std::ostream& operator<<(std::ostream& os, box3<T> const& b) {
  os << b.lower() << "<->" << b.upper();
  return os;
}

inline std::ostream& operator<<(std::ostream& os, grid3 const& g) {
  os << g.extents();
  return os;
}

inline std::ostream& operator<<(std::ostream& os, subgrid3 const& s) {
  os << s.lower() << "<->" << s.upper();
  return os;
}

template <class T, int N>
inline std::ostream& operator<<(std::ostream& os, static_vector<T, N> const& v) {
  os << "[";
  for (int i = 0; i < N-1; ++i) {
    os << v[i] << ", ";
  }
  os << v[N-1] << "]";
  return os;
}

}

namespace dgt {

inline std::ostream& operator<<(std::ostream& os, Point const& pt) {
  os << "["<< pt.depth << "], " << pt.ijk;
  return os;
}

}
