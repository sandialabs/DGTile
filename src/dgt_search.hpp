#pragma once

#include <algorithm>

#include "dgt_defines.hpp"
#include "dgt_macros.hpp"

namespace dgt {

template <class ArrayT>
DGT_METHOD inline int binary_search(real const x, ArrayT const& array)
{
  int const n = array.size();
  if (n == 1) return 0;
  int begin_idx = 0;
  int end_idx = n-1;
  while (std::abs(end_idx - begin_idx) > 1) {
    int const mid_idx = (end_idx + begin_idx) / 2;
    if (x < array[mid_idx]) end_idx = mid_idx;
    else begin_idx = mid_idx;
  }
  return begin_idx;
}

}
