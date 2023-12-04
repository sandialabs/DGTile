#pragma once

#include <algorithm>

#include "dgt_defines.hpp"
#include "dgt_macros.hpp"

namespace dgt {

template <class ArrayT>
DGT_METHOD inline int binary_search(ArrayT const& array, real const x)
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

namespace bilinear {

template <class ArrayT>
DGT_METHOD inline real interpolate(
    ArrayT const& f_values,
    ArrayT const& x_axis,
    ArrayT const& y_axis,
    real const x,
    real const y,
    int const i,
    int const j)
{
  int const nx = x_axis.size();
  real const dx = x_axis[i+1] - x_axis[i];
  real const dy = y_axis[j+1] - y_axis[j];
  real const qx = (x - x_axis[i])/dx;
  real const qy = (y - y_axis[j])/dy;
  real const rx = 1. - qx;
  real const ry = 1. - qy;
  real const f_i_j = f_values[j*nx + i];
  real const f_i_jp1 = f_values[(j+1)*nx + i];
  real const f_ip1_j = f_values[j*nx + (i+1)];
  real const f_ip1_jp1 = f_values[(j+1)*nx + (i+1)];
  real const f =
    f_i_j * rx * ry +
    f_i_jp1 * rx * qy +
    f_ip1_j * qx * ry +
    f_ip1_jp1 * qx * qy;
  return f;
}

}

}
