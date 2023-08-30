#pragma once

#include <Kokkos_Core.hpp>

#include "dgt_grid3.hpp"
#include "dgt_subgrid3.hpp"

namespace dgt {

template <class Functor>
DGT_ALWAYS_INLINE inline constexpr void seq_for_each(
    Subgrid3 const& subgrid,
    Functor const& functor)
{
  int const dim = infer_dimension(subgrid.upper());
  Subgrid3 const s = generalize(dim, subgrid);
  if (s.size() == 0) return;
  for (int k = s.lower().z(); k < s.upper().z(); ++k) {
    for (int j = s.lower().y(); j < s.upper().y(); ++j) {
      for (int i = s.lower().x(); i < s.upper().x(); ++i) {
        functor(Vec3<int>(i, j, k));
      }
    }
  }
}

}
