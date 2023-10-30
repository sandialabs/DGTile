#include <stdexcept>

#include "dgt_cartesian.hpp"
#include "dgt_tree.hpp"

namespace dgt {

static int infer_dim(Grid3 const& cell_grid) {
  int const dim = infer_dimension(cell_grid);
  if (dim < 0) {
    throw std::runtime_error("dgt:cartesian -> infer dimension failed");
  }
  return dim;
}

Grid3 get_face_grid(Grid3 const& cell_grid, int const axis)
{
  Vec3<int> const cell_extents = cell_grid.extents();
  Vec3<int> const face_extents = cell_extents + Vec3<int>::axis(axis);
  return Grid3(face_extents);
}

int get_num_cells(Grid3 const& cell_grid)
{
  int const dim = infer_dim(cell_grid);
  return generalize(dim, cell_grid).size();
}

int get_num_faces(Grid3 const& cell_grid, int const axis)
{
  int const dim = infer_dim(cell_grid);
  Grid3 const face_grid = get_face_grid(cell_grid, axis);
  return generalize(dim, face_grid).size();
}

Subgrid3 get_owned_cells(Grid3 const& cell_grid)
{
  int const dim = infer_dim(cell_grid);
  Vec3<int> lower = Vec3<int>::zero();
  Vec3<int> upper = Vec3<int>::zero();
  Vec3<int> const ncells = cell_grid.extents();
  for (int axis = 0; axis < dim; ++axis) {
    lower[axis] = 1;
    upper[axis] = ncells[axis] - 1;
  }
  return Subgrid3(lower, upper);
}

Subgrid3 get_owned_faces(Grid3 const& cell_grid, int const axis)
{
  Subgrid3 faces = get_owned_cells(cell_grid);
  faces.upper()[axis] += 1;
  return faces;
}

static Subgrid3 const get_cells(
    Grid3 const& cell_grid,
    int const ownership,
    int const axis,
    int const dir)
{
  int offset = (ownership == OWNED) ? 1 : 0;
  int const dim = infer_dim(cell_grid);
  Vec3<int> const ncells = cell_grid.extents();
  Vec3<int> lower = Vec3<int>::zero();
  Vec3<int> upper = Vec3<int>::zero();
  for (int d = 0; d < dim; ++d) {
    lower[d] = 1;
    upper[d] = ncells[d] - 1;
  }
  lower[axis] = (dir == LEFT) ? offset : ncells[axis]-offset-1;
  upper[axis] = (dir == LEFT) ? offset+1 : ncells[axis]-offset;
  return Subgrid3(lower, upper);
}

static Subgrid3 const get_ftc_cells(
    Grid3 const& cell_grid,
    int const ownership,
    int const axis,
    int const dir)
{
  Subgrid3 s = get_cells(cell_grid, ownership, axis, dir);
  if (ownership == GHOST) return s;
  if (dir == LEFT) s.upper()[axis] += 1;
  if (dir == RIGHT) s.lower()[axis] -= 1;
  return s;
}

static Subgrid3 const get_ctf_cells(
    Grid3 const& cell_grid,
    int const ownership,
    int const axis,
    int const dir,
    int const which_child)
{
  int const dim = infer_dim(cell_grid);
  Subgrid3 s = get_cells(cell_grid, ownership, axis, dir);
  Vec3<int> const child_ijk = get_child_ijk(which_child);
  Vec3<int> const ncells = cell_grid.extents();
  for (int d = 0; d < dim; ++d) {
    if (d == axis) continue;
    int const half_cells = ncells[d] / 2;
    if (child_ijk[axis] == 0) s.upper()[axis] = half_cells;
    if (child_ijk[axis] == 1) s.lower()[axis] = half_cells;
  }
  return s;
}

Subgrid3 get_cells(
    Grid3 const& cell_grid,
    int const ownership,
    int const adjacency_kind,
    int const axis,
    int const dir,
    int const which_child)
{
  if (adjacency_kind == tree::EQUAL) {
    return get_cells(cell_grid, ownership, axis, dir);
  } else if (adjacency_kind == tree::FINE_TO_COARSE) {
    return get_ftc_cells(cell_grid, ownership, axis, dir);
  } else if (adjacency_kind == tree::COARSE_TO_FINE) {
    return get_ctf_cells(cell_grid, ownership, axis, dir, which_child);
  } else {
    throw std::runtime_error("dgt:get_cells - invalid kind");
  }
}

}
