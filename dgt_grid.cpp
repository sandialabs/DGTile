#include <stdexcept>

#include "dgt_defines.hpp"
#include "dgt_grid.hpp"
#include "dgt_point.hpp"

namespace dgt {

int get_dim(p3a::vector3<int> const& v) {
  if ((v.x() > 0) && (v.y() > 0) && (v.z() > 0)) return 3;
  if ((v.x() > 0) && (v.y() > 0) && (v.z() == 0)) return 2;
  if ((v.x() > 0) && (v.y() == 0) && (v.z() == 0)) return 1;
  throw std::runtime_error("get_dim");
}

int get_dim(p3a::grid3 const& g) {
  return get_dim(g.extents());
}

p3a::grid3 generalize(p3a::grid3 const& g) {
  int const dim = get_dim(g);
  p3a::vector3<int> n = g.extents();
  for (int axis = dim; axis < DIMS; ++axis) {
    n[axis] = 1;
  }
  return p3a::grid3(n);
}

p3a::subgrid3 generalize(p3a::subgrid3 const& s) {
  p3a::vector3<int> start = s.lower();
  p3a::vector3<int> end = s.upper();
  int const dim = get_dim(end);
  for (int axis = dim; axis < DIMS; ++axis) {
    start[axis] = 0;
    end[axis] = 1;
  }
  return p3a::subgrid3(start, end);
}

p3a::grid3 get_child_grid(int dim) {
  if (dim == 1) return p3a::grid3(2,0,0);
  if (dim == 2) return p3a::grid3(2,2,0);
  if (dim == 3) return p3a::grid3(2,2,2);
  throw std::runtime_error("get_child_grid");
}

p3a::grid3 get_node_grid(p3a::grid3 const& cells) {
  int const dim = get_dim(cells);
  p3a::vector3<int> nodes = cells.extents();
  for (int axis = 0; axis < dim; ++axis) {
    nodes[axis] += 1;
  }
  return p3a::grid3(nodes);
}

p3a::grid3 get_side_grid(p3a::grid3 const& cells, int axis) {
  return p3a::grid3(cells.extents() + p3a::vector3<int>::axis(axis));
}

p3a::subgrid3 get_intr_sides(p3a::grid3 const& cells, int axis) {
  p3a::grid3 const sides = get_side_grid(cells, axis);
  p3a::vector3<int> const nsides = sides.extents();
  p3a::vector3<int> const start = p3a::vector3<int>::axis(axis);
  p3a::vector3<int> const end = nsides - p3a::vector3<int>::axis(axis);
  return p3a::subgrid3(start, end);
}

p3a::subgrid3 get_adj_sides(p3a::grid3 const& cells, int axis, int dir) {
  p3a::grid3 const sides = get_side_grid(cells, axis);
  p3a::vector3<int> const nsides = sides.extents();
  p3a::vector3<int> start = p3a::vector3<int>::zero();
  p3a::vector3<int> end = nsides;
  if (dir == left) end[axis] = 1;
  if (dir == right) start[axis] = nsides[axis] - 1;
  return p3a::subgrid3(start, end);
}

p3a::subgrid3 get_adj_cells(p3a::grid3 const& cells, int axis, int dir) {
  p3a::vector3<int> const ncells = cells.extents();
  p3a::vector3<int> start = p3a::vector3<int>::zero();
  p3a::vector3<int> end = ncells;
  if (dir == left) end[axis] = 1;
  if (dir == right) start[axis] = ncells[axis] - 1;
  return p3a::subgrid3(start, end);
}

static p3a::vector3<int> map_to_fine(Point const& pt, int fine_depth) {
  int const diff = fine_depth - pt.depth;
  return p3a::vector3<int>(
      pt.ijk.x() << diff,
      pt.ijk.y() << diff,
      pt.ijk.z() << diff);
}

static bool contains(p3a::subgrid3 const& s, p3a::vector3<int> const& ijk) {
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

bool contains(p3a::grid3 const& g, p3a::subgrid3 const& s) {
  bool in = true;
  if (get_dim(g) != get_dim(s.upper())) in = false;
  if (!contains(g, s.lower())) in = false;
  if (!contains(g, s.upper())) in = false;
  return in;
}

bool contains(int depth, p3a::subgrid3 const& s, Point const& pt) {
  if (depth == pt.depth) {
    return contains(s, pt.ijk);
  } else if (depth > pt.depth) {
    p3a::vector3<int> const fine_pt = map_to_fine(pt, depth);
    return contains(s, fine_pt);
  } else {
    p3a::vector3<int> const start = map_to_fine({depth, s.lower()}, pt.depth);
    p3a::vector3<int> const end = map_to_fine({depth, s.upper()}, pt.depth);
    return contains(p3a::subgrid3(start, end), pt.ijk);
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
  return (quotient * part) + p3a::min(remainder, part);
}

}
