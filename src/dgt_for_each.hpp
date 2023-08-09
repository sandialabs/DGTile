#pragma once

#include <Kokkos_Core.hpp>

#include "dgt_grid3.hpp"
#include "dgt_subgrid3.hpp"
#include "dgt_vec3.hpp"

namespace dgt {

template <class Functor>
DGT_ALWAYS_INLINE inline constexpr void seq_for_each(
    Subgrid3 const& subgrid,
    Functor const& functor)
{
  if (subgrid.size() < 0) return;
  Subgrid3 s = subgrid;
  if (s.upper().dimension() < 3) s.upper().z() = 1;
  if (s.upper().dimension() < 2) s.upper().y() = 1;
  for (int k = s.lower().z(); k < s.upper().z(); ++k) {
    for (int j = s.lower().y(); j < s.upper().y(); ++j) {
      for (int i = s.lower().x(); i < s.upper().x(); ++i) {
        functor(Vec3<int>(i, j, k));
      }
    }
  }
}

}
