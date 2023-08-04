#pragma once

#include "dgt_vec3.hpp"

namespace dgt {

DGT_METHOD inline int get_dir_sign(int const dir)
{
  static constexpr int sign[DIRECTIONS] = {-1, 1};
  return sign[dir];
}

DGT_METHOD inline int invert_dir(int const dir)
{
  static constexpr int inverse[DIRECTIONS] = {RIGHT, LEFT};
  return inverse[dir];
}

DGT_METHOD inline vec3<int> get_cells_adj_face(
    vec3<int> const& cell_ijk,
    int const axis,
    int const dir)
{
  static constexpr int offset[DIRECTIONS] = {0, 1};
  return cell_ijk + offset[dir] * vec3<int>::axis(axis);
}

DGT_METHOD inline vec3<int> get_faces_adj_cell(
    vec3<int> const& face_ijk,
    int const axis,
    int const dir)
{
  static constexpr int offset[DIRECTIONS] = {-1, 0};
  return face_ijk + offset[dir] * vec3<int>::axis(axis);
}

}
