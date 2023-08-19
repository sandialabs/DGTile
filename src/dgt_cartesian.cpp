#include "dgt_cartesian.hpp"
#include "dgt_tree.hpp"

namespace dgt {

Vec3<std::int8_t> get_adj_ijk_offset(Vec3<std::int8_t> const& ijk_offset)
{
  return std::int8_t(-1)*ijk_offset;
}

Subgrid3 get_cells(
    int const ownership,
    Grid3 const& cell_grid,
    Vec3<std::int8_t> const& ijk_offset)
{
  int const own_off = (ownership == OWNED) ? 1 : 0;
  int const dim = cell_grid.dimension();
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
  int const dim = cell_grid.dimension();
  Subgrid3 s = get_cells(OWNED, cell_grid, ijk_offset);
  for (int axis = 0; axis < dim; ++axis) {
    if (ijk_offset[axis] == -1) s.upper()[axis] += 1;
    if (ijk_offset[axis] ==  1) s.lower()[axis] -= 1;
  }
  return s;
}

}
