#include "dgt_cartesian.hpp"
#include "dgt_tree.hpp"

namespace dgt {

Vec3<std::int8_t> get_adj_ijk_offset(Vec3<std::int8_t> const& ijk_offset)
{
  return std::int8_t(-1)*ijk_offset;
}

static int infer_dim(Grid3 const& cell_grid) {
  int const dim = infer_dimension(cell_grid);
  if (dim < 0) {
    throw std::runtime_error("dgt:cartesian -> infer dimension failed");
  }
  return dim;
}

int get_num_cells(Grid3 const& cell_grid)
{
  int const dim = infer_dim(cell_grid);
  return generalize(dim, cell_grid).size();
}

Grid3 get_face_grid(Grid3 const& cell_grid, int const axis)
{
  Vec3<int> const cell_extents = cell_grid.extents();
  Vec3<int> const face_extents = cell_extents + Vec3<int>::axis(axis);
  return Grid3(face_extents);
}

int get_num_faces(Grid3 const& cell_grid, int const axis)
{
  int const dim = infer_dim(cell_grid);
  Grid3 const face_grid = get_face_grid(cell_grid, axis);
  return generalize(dim, face_grid).size();
}

Subgrid3 get_cells(
    int const ownership,
    Grid3 const& cell_grid,
    Vec3<std::int8_t> const& ijk_offset)
{
  int const own_off = (ownership == OWNED) ? 1 : 0;
  int const dim = infer_dim(cell_grid);
  Vec3<int> const ncells = cell_grid.extents();
  Vec3<int> lower = Vec3<int>::zero();
  Vec3<int> upper = Vec3<int>::zero();
  for (int axis = 0; axis < dim; ++axis) {
    int const last = ncells[axis];
    int const lower_offsets[3] = {own_off, 1, last-1-own_off};
    int const upper_offsets[3] = {1+own_off, last-1, last-own_off};
    lower[axis] = lower_offsets[ijk_offset[axis] + 1];
    upper[axis] = upper_offsets[ijk_offset[axis] + 1];
  }
  return Subgrid3(lower, upper);
}

Subgrid3 get_fine_to_coarse_cells(
    int const ownership,
    Grid3 const& cell_grid,
    Vec3<std::int8_t> const& ijk_offset)
{
  if (ownership == GHOST) {
    return get_cells(GHOST, cell_grid, ijk_offset);
  }
  int const dim = infer_dim(cell_grid);
  Subgrid3 s = get_cells(OWNED, cell_grid, ijk_offset);
  for (int axis = 0; axis < dim; ++axis) {
    if (ijk_offset[axis] == -1) s.upper()[axis] += 1;
    if (ijk_offset[axis] ==  1) s.lower()[axis] -= 1;
  }
  return s;
}

Subgrid3 get_coarse_to_fine_cells(
    int const ownership,
    Grid3 const& cell_grid,
    Vec3<int> const& child_ijk,
    Vec3<std::int8_t> const& ijk_offset)
{
  int const dim = infer_dim(cell_grid);
  Vec3<int> const ncells = cell_grid.extents();
  Subgrid3 s = get_cells(ownership, cell_grid, ijk_offset);
  for (int axis = 0; axis < dim; ++axis) {
    if (ijk_offset[axis] == 0) {
      int const axis_half_cells = ncells[axis] / 2;
      if (child_ijk[axis] == 0) s.upper()[axis] = axis_half_cells;
      if (child_ijk[axis] == 1) s.lower()[axis] = axis_half_cells;
    }
  }
  return s;
}

}
