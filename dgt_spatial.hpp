#pragma once

#include <cmath>

#include "p3a_grid3.hpp"

#include "dgt_defines.hpp"
#include "dgt_point.hpp"
#include "dgt_basis.hpp"

namespace dgt {

[[nodiscard]] P3A_ALWAYS_INLINE P3A_HOST_DEVICE inline
int get_dir_sign(int dir) {
  return (dir == left) ? -1 : 1;
}

[[nodiscard]] P3A_ALWAYS_INLINE P3A_HOST_DEVICE inline
int invert_dir(int dir) {
  return (dir == left) ? right : left;
}

[[nodiscard]] P3A_ALWAYS_INLINE P3A_HOST_DEVICE inline
p3a::vector3<int> get_sides_adj_cell(
    p3a::vector3<int> const& side, int axis, int dir) {
  static constexpr int offset[ndirs] = {-1,0};
  return side + offset[dir]*p3a::vector3<int>::axis(axis);
}

[[nodiscard]] P3A_ALWAYS_INLINE P3A_HOST_DEVICE inline
p3a::vector3<int> get_cells_adj_side(
    p3a::vector3<int> const& cell, int axis, int dir) {
  static constexpr int offset[ndirs] = {0,1};
  return cell + offset[dir]*p3a::vector3<int>::axis(axis);
}

[[nodiscard]] P3A_ALWAYS_INLINE P3A_HOST_DEVICE inline
p3a::vector3<int> get_cells_adj_cell(
    p3a::vector3<int> const& cell, int axis, int dir) {
  static constexpr int offset[ndirs] = {-1,1};
  return cell + offset[dir]*p3a::vector3<int>::axis(axis);
}

[[nodiscard]] P3A_ALWAYS_INLINE P3A_HOST_DEVICE inline
p3a::vector3<int> get_border_ijk(p3a::vector3<int> const& ijk, int axis) {
  p3a::vector3<int> border_ijk = ijk;
  border_ijk[axis] = 0;
  return border_ijk;
}

[[nodiscard]] P3A_ALWAYS_INLINE P3A_HOST_DEVICE inline
p3a::vector3<double> get_dx(
    p3a::box3<double> const& domain, p3a::grid3 const& nparts) {
  p3a::vector3<double> dx(0,0,0);
  p3a::vector3<double> const l = domain.extents();
  p3a::vector3<int> const n = nparts.extents();
  for (int axis = 0; axis < DIMS; ++axis) {
    if (n[axis] == 0) continue;
    dx[axis] = l[axis] / n[axis];
  }
  return dx;
}

[[nodiscard]] P3A_ALWAYS_INLINE P3A_HOST_DEVICE inline
p3a::box3<double> get_block_domain(
    Point const& base,
    Point const& node,
    p3a::box3<double> const& domain) {
  int const diff = node.depth - base.depth;
  p3a::vector3<int> const nbase = base.ijk;
  p3a::vector3<double> const base_dx = get_dx(domain, nbase);
  p3a::vector3<double> const dx = base_dx / std::pow(2, diff);
  p3a::vector3<double> const xmin = domain.lower();
  p3a::vector3<double> const start = xmin + p3a::hadamard_product(node.ijk, dx);
  p3a::vector3<double> const end = start + dx;
  return p3a::box3<double>(start, end);
}

[[nodiscard]] P3A_ALWAYS_INLINE P3A_HOST_DEVICE inline
p3a::vector3<double> get_cell_center(
    p3a::vector3<int> const& cell_ijk,
    p3a::vector3<double> const& origin,
    p3a::vector3<double> const& dx) {
  return p3a::vector3<double>(
      origin.x() + (cell_ijk.x() + 0.5) * dx.x(),
      origin.y() + (cell_ijk.y() + 0.5) * dx.y(),
      origin.z() + (cell_ijk.z() + 0.5) * dx.z());
}

[[nodiscard]] P3A_ALWAYS_INLINE P3A_HOST_DEVICE inline
p3a::vector3<double> get_x(
    p3a::vector3<int> const& cell_ijk,
    p3a::vector3<double> const& origin,
    p3a::vector3<double> const& dx,
    p3a::vector3<double> const& xi) {
  p3a::vector3<double> const center = get_cell_center(cell_ijk, origin, dx);
  return p3a::vector3<double>(
      center.x() + 0.5 * xi.x() * dx.x(),
      center.y() + 0.5 * xi.y() * dx.y(),
      center.z() + 0.5 * xi.z() * dx.z());
}

template <class BasisT>
[[nodiscard]] P3A_ALWAYS_INLINE P3A_HOST_DEVICE inline
p3a::vector3<double> get_intr_x(
    BasisT const& b,
    int pt,
    p3a::vector3<int> const& cell_ijk,
    p3a::vector3<double> const& origin,
    p3a::vector3<double> const& dx) {
  p3a::vector3<double> const xi = get_intr_pt(b, pt);
  return get_x(cell_ijk, origin, dx, xi);
}

template <class BasisT>
[[nodiscard]] P3A_ALWAYS_INLINE P3A_HOST_DEVICE inline
p3a::vector3<double> get_fine_x(
    BasisT const& b,
    int pt,
    p3a::vector3<int> const& cell_ijk,
    p3a::vector3<double> const& origin,
    p3a::vector3<double> const& dx) {
  p3a::vector3<double> const xi = get_fine_pt(b, pt);
  return get_x(cell_ijk, origin, dx, xi);
}

template <class BasisT>
[[nodiscard]] P3A_ALWAYS_INLINE P3A_HOST_DEVICE inline
p3a::vector3<double> get_side_x(
    BasisT const& b,
    int axis,
    int pt,
    p3a::vector3<int> const& side_ijk,
    p3a::vector3<double> const& origin,
    p3a::vector3<double> const& dx) {
  p3a::vector3<double> const xi = get_side_pt(b, axis, left, pt);
  return get_x(side_ijk, origin, dx, xi);
}

[[nodiscard]] P3A_ALWAYS_INLINE P3A_HOST_DEVICE inline
double get_volume(int dim, p3a::vector3<double> const& v) {
  double volume = v.x();
  if (dim > 1) volume *= v.y();
  if (dim > 2) volume *= v.z();
  return volume;
}

[[nodiscard]] P3A_ALWAYS_INLINE P3A_HOST_DEVICE inline
double get_cell_detJ(int dim, p3a::vector3<double> const& dx) {
  return std::pow(0.5, dim) * get_volume(dim, dx);
}

[[nodiscard]] P3A_ALWAYS_INLINE P3A_HOST_DEVICE inline
double get_side_detJ(int dim, int axis, p3a::vector3<double> const& dx) {
  return std::pow(0.5, dim-1) * get_volume(dim, dx) / dx[axis];
}

[[nodiscard]] P3A_ALWAYS_INLINE P3A_HOST_DEVICE inline
double get_amr_side_detJ(int dim, int axis, p3a::vector3<double> const& dx) {
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
