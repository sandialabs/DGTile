#pragma once

#include <cstdint>
#include <string>

#include "dgt_subgrid3.hpp"

namespace dgt {

static constexpr Grid3 child_grid = {2,2,2};

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

DGT_METHOD inline Vec3<int> get_cells_adj_face(
    Vec3<int> const& cell_ijk,
    int const axis,
    int const dir)
{
  static constexpr int offset[DIRECTIONS] = {0, 1};
  return cell_ijk + offset[dir] * Vec3<int>::axis(axis);
}

DGT_METHOD inline Vec3<int> get_faces_adj_cell(
    Vec3<int> const& face_ijk,
    int const axis,
    int const dir)
{
  static constexpr int offset[DIRECTIONS] = {-1, 0};
  return face_ijk + offset[dir] * Vec3<int>::axis(axis);
}

DGT_METHOD inline int permute(int const ijk, int const axis)
{
  return
    (ijk == X) ? axis :
    (ijk == Y) ? (axis+1) % DIMENSIONS :
    (ijk == Z) ? (axis+2) % DIMENSIONS :
    -1;
}

DGT_METHOD inline Vec3<real> get_cell_center(
    Vec3<int> const& cell_ijk,
    Vec3<real> const& origin,
    Vec3<real> const& dx)
{
  return Vec3<real>(
      origin.x() + (cell_ijk.x() + 0.5) * dx.x(),
      origin.y() + (cell_ijk.y() + 0.5) * dx.y(),
      origin.z() + (cell_ijk.z() + 0.5) * dx.z());
}

DGT_METHOD inline Vec3<real> map_to_physical(
    Vec3<int> const& cell_ijk,
    Vec3<real> const& origin,
    Vec3<real> const& dx,
    Vec3<real> const& xi)
{
  Vec3<real> const center = get_cell_center(cell_ijk, origin, dx);
  return Vec3<real>(
      center.x() + 0.5 * xi.x() * dx.x(),
      center.y() + 0.5 * xi.y() * dx.y(),
      center.z() + 0.5 * xi.z() * dx.z());
}

DGT_METHOD inline bool is_cell_ghost(
    int const dim,
    Grid3 const& cell_grid,
    Vec3<int> const& cell_ijk)
{
  for (int axis = 0; axis < dim; ++axis) {
    int const min = 0;
    int const max = cell_grid.extents()[axis] - 1;
    if (cell_ijk[axis] == min) return true;
    if (cell_ijk[axis] == max) return true;
  }
  return false;
}

DGT_METHOD inline Vec3<int> get_child_ijk(int const which_child)
{
  Vec3<int> child_ijk;
  int const a = which_child % 2;
  int const b = (which_child - a) / 2;
  child_ijk.x() = a;
  child_ijk.y() = b % 2;
  child_ijk.z() = b / 2;
  return child_ijk;
}

DGT_METHOD inline int get_which_child(Vec3<int> const& child_ijk)
{
  return child_ijk.x() + 2 * (child_ijk.y() + 2 * child_ijk.z());
}

std::string inline get_axis_name(int const axis)
{
  return
    (axis == X) ? "X" :
    (axis == Y) ? "Y" :
    (axis == Z) ? "Z" :
    "";
}

[[nodiscard]] int get_num_cells(Grid3 const& cell_grid);

[[nodiscard]] Grid3 get_face_grid(Grid3 const& cell_grid, int const axis);
[[nodiscard]] int get_num_faces(Grid3 const& cell_grid, int const axis);

[[nodiscard]] Subgrid3 get_owned_cells(Grid3 const& cell_grid);
[[nodiscard]] Subgrid3 get_owned_faces(Grid3 const& cell_grid, int const axis);

[[nodiscard]] Subgrid3 get_cells(
    Grid3 const& cell_grid,
    int const ownership,
    int const adjacency_kind,
    int const axis,
    int const dir,
    int const which_child = 0);

}
