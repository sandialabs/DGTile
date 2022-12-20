#pragma once

#include <cmath>

#include "p3a_grid3.hpp"

#include "dgt_defines.hpp"
#include "dgt_point.hpp"

namespace dgt {

using namespace p3a;

[[nodiscard]] P3A_ALWAYS_INLINE P3A_HOST_DEVICE inline
int get_dir_sign(int dir) {
  return (dir == left) ? -1 : 1;
}

[[nodiscard]] P3A_ALWAYS_INLINE P3A_HOST_DEVICE inline
int invert_dir(int dir) {
  return (dir == left) ? right : left;
}

[[nodiscard]] P3A_ALWAYS_INLINE P3A_HOST_DEVICE inline
vector3<int> get_sides_adj_cell(vector3<int> const& side, int axis, int dir) {
  static constexpr int offset[ndirs] = {-1,0};
  return side + offset[dir]*vector3<int>::axis(axis);
}

[[nodiscard]] P3A_ALWAYS_INLINE P3A_HOST_DEVICE inline
vector3<int> get_cells_adj_side(vector3<int> const& cell, int axis, int dir) {
  static constexpr int offset[ndirs] = {0,1};
  return cell + offset[dir]*vector3<int>::axis(axis);
}

[[nodiscard]] P3A_ALWAYS_INLINE P3A_HOST_DEVICE inline
vector3<int> get_cells_adj_cell(vector3<int> const& cell, int axis, int dir) {
  static constexpr int offset[ndirs] = {-1,1};
  return cell + offset[dir]*vector3<int>::axis(axis);
}

[[nodiscard]] P3A_ALWAYS_INLINE P3A_HOST_DEVICE inline
vector3<int> get_border_ijk(vector3<int> const& ijk, int axis) {
  vector3<int> border_ijk = ijk;
  border_ijk[axis] = 0;
  return border_ijk;
}

[[nodiscard]] P3A_ALWAYS_INLINE P3A_HOST_DEVICE inline
vector3<double> get_dx(box3<double> const& domain, grid3 const& nparts) {
  vector3<double> dx(0,0,0);
  vector3<double> const l = domain.extents();
  vector3<int> const n = nparts.extents();
  for (int axis = 0; axis < DIMS; ++axis) {
    if (n[axis] == 0) continue;
    dx[axis] = l[axis] / n[axis];
  }
  return dx;
}

[[nodiscard]] P3A_ALWAYS_INLINE P3A_HOST_DEVICE inline
box3<double> get_block_domain(
    Point const& base,
    Point const& node,
    box3<double> const& domain) {
  int const diff = node.depth - base.depth;
  vector3<int> const nbase = base.ijk;
  vector3<double> const base_dx = get_dx(domain, nbase);
  vector3<double> const dx = base_dx / std::pow(2, diff);
  vector3<double> const xmin = domain.lower();
  vector3<double> const start = xmin + hadamard_product(node.ijk, dx);
  vector3<double> const end = start + dx;
  return box3<double>(start, end);
}

[[nodiscard]] P3A_ALWAYS_INLINE P3A_HOST_DEVICE inline
vector3<double> get_cell_center(
    vector3<int> const& cell_ijk,
    vector3<double> const& origin,
    vector3<double> const& dx) {
  return vector3<double>(
      origin.x() + (cell_ijk.x() + 0.5) * dx.x(),
      origin.y() + (cell_ijk.y() + 0.5) * dx.y(),
      origin.z() + (cell_ijk.z() + 0.5) * dx.z());
}

[[nodiscard]] P3A_ALWAYS_INLINE P3A_HOST_DEVICE inline
vector3<double> get_x(
    vector3<int> const& cell_ijk,
    vector3<double> const& origin,
    vector3<double> const& dx,
    vector3<double> const& xi) {
  vector3<double> const center = get_cell_center(cell_ijk, origin, dx);
  return vector3<double>(
      center.x() + 0.5 * xi.x() * dx.x(),
      center.y() + 0.5 * xi.y() * dx.y(),
      center.z() + 0.5 * xi.z() * dx.z());
}

[[nodiscard]] P3A_ALWAYS_INLINE P3A_HOST_DEVICE inline
double get_volume(int dim, vector3<double> const& v) {
  double volume = v.x();
  if (dim > 1) volume *= v.y();
  if (dim > 2) volume *= v.z();
  return volume;
}

[[nodiscard]] P3A_ALWAYS_INLINE P3A_HOST_DEVICE inline
double get_cell_detJ(int dim, vector3<double> const& dx) {
  return std::pow(0.5, dim) * get_volume(dim, dx);
}

[[nodiscard]] P3A_ALWAYS_INLINE P3A_HOST_DEVICE inline
double get_side_detJ(int dim, int axis, vector3<double> const& dx) {
  return std::pow(0.5, dim-1) * get_volume(dim, dx) / dx[axis];
}

[[nodiscard]] P3A_ALWAYS_INLINE P3A_HOST_DEVICE inline
double get_amr_side_detJ(int dim, int axis, vector3<double> const& dx) {
  return std::pow(0.5, dim-1) * get_side_detJ(dim, axis, dx);
}

[[nodiscard]] P3A_ALWAYS_INLINE P3A_HOST_DEVICE inline
int permute(int ijk, int axis) {
  return
    (ijk == X) ? axis :
    (ijk == Y) ? (axis + 1) % DIMS :
    (ijk == Z) ? (axis + 2) % DIMS :
    -1;
}

}
