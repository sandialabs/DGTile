#include <stdexcept>

#include "dgt_defines.hpp"
#include "dgt_grid.hpp"
#include "dgt_point.hpp"

namespace dgt {

int get_dim(vector3<int> const& v) {
  if ((v.x() > 0) && (v.y() > 0) && (v.z() > 0)) return 3;
  if ((v.x() > 0) && (v.y() > 0) && (v.z() == 0)) return 2;
  if ((v.x() > 0) && (v.y() == 0) && (v.z() == 0)) return 1;
  throw std::runtime_error("get_dim");
}

int get_dim(grid3 const& g) {
  return get_dim(g.extents());
}

grid3 generalize(grid3 const& g) {
  int const dim = get_dim(g);
  vector3<int> n = g.extents();
  for (int axis = dim; axis < DIMS; ++axis) {
    n[axis] = 1;
  }
  return grid3(n);
}

subgrid3 generalize(subgrid3 const& s) {
  vector3<int> start = s.lower();
  vector3<int> end = s.upper();
  int const dim = get_dim(end);
  for (int axis = dim; axis < DIMS; ++axis) {
    start[axis] = 0;
    end[axis] = 1;
  }
  return subgrid3(start, end);
}

grid3 get_child_grid(int dim) {
  if (dim == 1) return grid3(2,0,0);
  if (dim == 2) return grid3(2,2,0);
  if (dim == 3) return grid3(2,2,2);
  throw std::runtime_error("get_child_grid");
}

grid3 get_side_grid(grid3 const& cells, int axis) {
  return grid3(cells.extents() + vector3<int>::axis(axis));
}

subgrid3 get_intr_sides(grid3 const& cells, int axis) {
  grid3 const sides = get_side_grid(cells, axis);
  vector3<int> const nsides = sides.extents();
  vector3<int> const start = vector3<int>::axis(axis);
  vector3<int> const end = nsides - vector3<int>::axis(axis);
  return subgrid3(start, end);
}

subgrid3 get_adj_sides(grid3 const& cells, int axis, int dir) {
  grid3 const sides = get_side_grid(cells, axis);
  vector3<int> const nsides = sides.extents();
  vector3<int> start = vector3<int>::zero();
  vector3<int> end = nsides;
  if (dir == left) end[axis] = 1;
  if (dir == right) start[axis] = nsides[axis] - 1;
  return subgrid3(start, end);
}

subgrid3 get_adj_cells(grid3 const& cells, int axis, int dir) {
  vector3<int> const ncells = cells.extents();
  vector3<int> start = vector3<int>::zero();
  vector3<int> end = ncells;
  if (dir == left) end[axis] = 1;
  if (dir == right) start[axis] = ncells[axis] - 1;
  return subgrid3(start, end);
}

grid3 get_viz_cell_grid(grid3 const& cells, int p) {
  vector3<int> const ncells = (p+1)*cells.extents();
  return grid3(ncells);
}

grid3 get_viz_point_grid(grid3 const& cells, int p) {
  int const dim = get_dim(cells);
  vector3<int> npoints = get_viz_cell_grid(cells, p).extents();
  for (int axis = 0; axis < dim; ++axis) {
    npoints[axis] += 1;
  }
  return grid3(npoints);
}

static vector3<int> map_to_fine(Point const& pt, int fine_depth) {
  int const diff = fine_depth - pt.depth;
  return vector3<int>(
      pt.ijk.x() << diff,
      pt.ijk.y() << diff,
      pt.ijk.z() << diff);
}

static bool contains(subgrid3 const& s, vector3<int> const& ijk) {
  bool in = true;
  int const dim = get_dim(s.upper());
  for (int axis = 0; axis < dim; ++axis) {
    if (ijk[axis] < s.lower()[axis]) in = false;
    if (ijk[axis] >= s.upper()[axis]) in = false;
  }
  for (int axis = dim; axis < DIMS; ++axis) {
    if (ijk[axis] != 0) in = false;
  }
  return in;
}

bool contains(grid3 const& g, subgrid3 const& s) {
  bool in = true;
  if (get_dim(g) != get_dim(s.upper())) in = false;
  if (!contains(g, s.lower())) in = false;
  if (!contains(g, s.upper())) in = false;
  return in;
}

bool contains(int depth, subgrid3 const& s, Point const& pt) {
  if (depth == pt.depth) {
    return contains(s, pt.ijk);
  } else if (depth > pt.depth) {
    vector3<int> const fine_pt = map_to_fine(pt, depth);
    return contains(s, fine_pt);
  } else {
    vector3<int> const start = map_to_fine({depth, s.lower()}, pt.depth);
    vector3<int> const end = map_to_fine({depth, s.upper()}, pt.depth);
    return contains(subgrid3(start, end), pt.ijk);
  }
}

int get_num_local(int ntotal, int nparts, int part) {
  int quotient = ntotal / nparts;
  int remainder = ntotal % nparts;
  return (part < remainder) ? (quotient + 1) : quotient;
}

int get_local_offset(int ntotal, int nparts, int part) {
  int quotient = ntotal / nparts;
  int remainder = ntotal % nparts;
  return (quotient * part) + min(remainder, part);
}

}
