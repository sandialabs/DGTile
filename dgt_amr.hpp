#pragma once

#include "p3a_grid3.hpp"

#include "dgt_defines.hpp"
#include "dgt_point.hpp"
#include "dgt_views.hpp"

namespace dgt {

struct Basis;
class Mesh;
class Node;

[[nodiscard]] P3A_ALWAYS_INLINE P3A_HOST_DEVICE inline
vector3<int> get_local(int which_child) {
  vector3<int> local;
  int const a = which_child % 2;
  int const b = (which_child - a) / 2;
  local.x() = a;
  local.y() = b % 2;
  local.z() = b / 2;
  return local;
}

[[nodiscard]] P3A_ALWAYS_INLINE P3A_HOST_DEVICE inline
int get_which_child(vector3<int> const& local) {
  return local.x() + 2 * (local.y() + 2 * local.z());
}

[[nodiscard]] P3A_ALWAYS_INLINE P3A_HOST_DEVICE inline
vector3<int> get_local(int axis, int which_child) {
  vector3<int> local(0,0,0);
  int const ii = (axis == X) ? Y : X;
  int const jj = (axis == Z) ? Y : Z;
  local[ii] = which_child % 2;
  local[jj] = (which_child - local[ii]) / 2;
  return local;
}

[[nodiscard]] P3A_ALWAYS_INLINE P3A_HOST_DEVICE inline
int get_which_child(int axis, vector3<int> const& local) {
  int const ii = (axis == X) ? Y : X;
  int const jj = (axis == Z) ? Y : Z;
  return local[ii] + local[jj] * 2;
}

[[nodiscard]] P3A_ALWAYS_INLINE P3A_HOST_DEVICE inline
vector3<int> get_coarse_ijk(vector3<int> const& fine_ijk) {
  return vector3<int>(
      fine_ijk.x() >> 1,
      fine_ijk.y() >> 1,
      fine_ijk.z() >> 1);
}

[[nodiscard]] P3A_ALWAYS_INLINE P3A_HOST_DEVICE inline
vector3<int> get_fine_ijk(
    vector3<int> const& coarse_ijk,
    vector3<int> const& local) {
  return 2*coarse_ijk + local;
}

[[nodiscard]] P3A_ALWAYS_INLINE P3A_HOST_DEVICE inline
vector3<int> get_local_from_fine_ijk(vector3<int> const& fine_ijk) {
  vector3<int> const coarse_ijk = get_coarse_ijk(fine_ijk);
  vector3<int> const fine_corner_ijk = coarse_ijk * 2;
  vector3<int> const local = fine_ijk - fine_corner_ijk;
  return local;
}

[[nodiscard]] P3A_ALWAYS_INLINE P3A_HOST_DEVICE inline
vector3<int> map_to_fine(
    vector3<int> const& coarse_ijk,
    vector3<int> const& local,
    grid3 const& grid) {
  vector3<int> const n = grid.extents();
  return coarse_ijk + hadamard_product(local, n);
}

[[nodiscard]] P3A_ALWAYS_INLINE P3A_HOST_DEVICE inline
Point get_child_point(Point const& pt, vector3<int> const& local) {
  Point child_pt;
  child_pt.depth = pt.depth + 1;
  child_pt.ijk = get_fine_ijk(pt.ijk, local);
  return child_pt;
}

[[nodiscard]] P3A_ALWAYS_INLINE P3A_HOST_DEVICE inline
Point get_parent_point(Point const& pt) {
  Point parent_pt;
  parent_pt.depth = pt.depth - 1;
  parent_pt.ijk = get_coarse_ijk(pt.ijk);
  return parent_pt;
}

[[nodiscard]] P3A_ALWAYS_INLINE P3A_HOST_DEVICE inline
subgrid3 get_local_subgrid(grid3 const& grid, vector3<int> const& local) {
  vector3<int> const half = get_coarse_ijk(grid.extents());
  vector3<int> const start = hadamard_product(local, half);
  vector3<int> const end = start + half;
  return subgrid3(start, end);
}

void refine(int dim, Node* parent);
void coarsen(int dim, Node* parent);

void do_insertion(
    View<double***> from,
    View<double***> to,
    grid3 const& from_grid,
    grid3 const& to_grid,
    subgrid3 const& from_subgrid,
    subgrid3 const& to_subgrid);

void do_prolongation(
    Basis const& b,
    View<double***> from,
    View<double***> to,
    grid3 const& from_grid,
    grid3 const& to_grid,
    subgrid3 const& from_subgrid,
    subgrid3 const& to_subgrid);

void do_restriction(
    Basis const& b,
    View<double***> from,
    View<double***> to,
    grid3 const& from_grid,
    grid3 const& to_grid,
    subgrid3 const& from_subgrid,
    subgrid3 const& to_fine_subgrid);

void modify(Mesh& mesh, std::vector<int8_t> const& marks);

}
