#include <spdlog/spdlog.h>

#include "dgt_cartesian.hpp"
#include "dgt_tree.hpp"

namespace dgt {

static int get_offset(int const distribution)
{
  if (distribution == OWNED) return 1;
  else return 0;
}

Subgrid3 get_equal_cells(
    Grid3 const& cell_grid,
    Vec3<std::int8_t> const& meta_ijk,
    int const offset)
{
  Vec3<int> lower(0,0,0);
  Vec3<int> upper(0,0,0);
  int const dim = cell_grid.dimension();
  Vec3<int> const ncells = cell_grid.extents();
  for (int axis = 0; axis < dim; ++axis) {
    int const last = ncells[axis] - offset;
    int const lower_offsets[3] = {offset, offset+1, last-1};
    int const upper_offsets[3] = {offset+1, last-1, last};
    lower[axis] = lower_offsets[meta_ijk[axis]];
    upper[axis] = upper_offsets[meta_ijk[axis]];
  }
  return Subgrid3(lower, upper);
}

Subgrid3 get_cells(
    int const distribution,
    Grid3 const& cell_grid,
    Vec3<std::int8_t> const& meta_ijk,
    std::int8_t const level_difference)
{
  int const offset = get_offset(distribution);
  if (level_difference == tree::EQUAL) {
    return get_equal_cells(cell_grid, meta_ijk, offset);
  } else {
    spdlog::error("dgt:get_cells - invalid level difference");
    throw std::runtime_error("dgt:get_cells");
  }
}

}
