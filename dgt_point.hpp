#pragma once

#include "p3a_vector3.hpp"

namespace dgt {

struct Point {
  int depth;
  p3a::vector3<int> ijk;
  [[nodiscard]] P3A_ALWAYS_INLINE P3A_HOST_DEVICE constexpr
  bool operator==(Point const& other) const {
    return (depth == other.depth) && (ijk == other.ijk);
  }
  [[nodiscard]] P3A_ALWAYS_INLINE P3A_HOST_DEVICE constexpr
  bool operator!=(Point const& other) const {
    return !operator==(other);
  }
};

}
