#include <list>

#include "dgt_cartesian.hpp"
#include "dgt_for_each.hpp"
#include "dgt_tree.hpp"

namespace dgt {
namespace tree {

ID get_level_offset(int const dim, int const level)
{
  ID const num = (ID(1) << (dim * level)) - 1;
  ID const den = (ID(1) << dim) - 1;
  return num / den;
}

ID get_level_id(int const dim, Point const& pt)
{
  (void)dim;
  ID const L = pt.level;
  ID const i = pt.ijk.x();
  ID const j = pt.ijk.y();
  ID const k = pt.ijk.z();
  ID const idx1 = j << L;
  ID const idx2 = k << (2*L);
  return i + idx1 + idx2;
}

ID get_global_id(int const dim, Point const& pt)
{
  ID const level_offset = get_level_offset(dim, pt.level);
  ID const level_id = get_level_id(dim, pt);
  return level_offset + level_id;
}

std::int8_t get_level(int const dim, ID const global_id)
{
  int L = 0;
  while (global_id >= get_level_offset(dim, L+1)) {
    L++;
  }
  return std::int8_t(L);
}

Point get_point(int const dim, ID const global_id)
{
  int const L = get_level(dim, global_id);
  int const grid_extent = ID(1) << L;
  ID const level_offset = get_level_offset(dim, L);
  ID idx = global_id - level_offset;
  int const i = idx % grid_extent;
  idx = idx >> L;
  int const j = idx % grid_extent;
  idx = idx >> L;
  int const k = idx;
  Vec3<int> const ijk(i,j,k);
  return Point(L, dimensionalize(dim, ijk));
}

Point get_coarse_point(int const dim, Point const& pt)
{
  int const level = pt.level - 1;
  Vec3<int> const ijk = pt.ijk / 2;
  return Point(level, dimensionalize(dim,ijk));
}

Point get_fine_point(int const dim, Point const& pt, Vec3<int> const& child_ijk)
{
  int const level = pt.level + 1;
  Vec3<int> const ijk = pt.ijk*2 + child_ijk;
  return Point(level, dimensionalize(dim, ijk));
}

static int get_level(Grid3 const& g)
{
  int level = 0;
  int const gmax = max(g.extents());
  while ((1 << level) < gmax) {
    level++;
  }
  return level;
}

Leaves create(int const dim, Grid3 const& grid)
{
  Leaves leaves;
  int const level = get_level(grid);
  auto functor = [&] (Vec3<int> const& ijk) {
    Point const pt(level, ijk);
    ID const global_id = get_global_id(dim, pt);
    leaves.insert(global_id);
  };
  seq_for_each(dimensionalize(dim, grid), functor);
  return leaves;
}

static bool is_leaf(ID const id, Leaves const& leaves)
{
  return leaves.count(id);
}

static void recursively_order(
    ZLeaves& z_leaves,
    int const dim,
    Leaves const& leaves,
    int const max_level,
    ID const global_id)
{
  if (is_leaf(global_id, leaves)) {
    z_leaves.push_back(global_id);
    return;
  }
  int const level = get_level(dim, global_id);
  if (level == max_level) return;
  Point const pt = get_point(dim, global_id);
  auto functor = [&] (Vec3<int> const& child_ijk) {
    Point const fine_pt = get_fine_point(dim, pt, child_ijk);
    ID const fine_global_id = get_global_id(dim, fine_pt);
    recursively_order(z_leaves, dim, leaves, max_level, fine_global_id);
  };
  seq_for_each(dimensionalize(dim, child_grid), functor);
}

ZLeaves order(int const dim, Leaves const& leaves)
{
  int const max_level = get_max_level(dim, leaves);
  ZLeaves z_leaves;
  recursively_order(z_leaves, dim, leaves, max_level, 0);
  return z_leaves;
}

int min_op(int const a, int const b) { return (a < b) ? a : b; }
int max_op(int const a, int const b) { return (a > b) ? a : b; }

template <class LeavesT, class Op>
int reduce_level(
    int const dim,
    LeavesT const& leaves,
    Op const& op,
    int const start)
{
  int level = start;
  for (ID const global_id : leaves) {
    int const leaf_level = get_level(dim, global_id);
    level = op(level, leaf_level);
  }
  return level;
}

template <class LeavesT>
int get_min_level(int const dim, LeavesT const& leaves)
{
  return reduce_level(dim, leaves, min_op, INT_MAX);
}

template int get_min_level(int const, Leaves const&);
template int get_min_level(int const, ZLeaves const&);

template <class Leaves>
int get_max_level(int const dim, Leaves const& leaves)
{
  return reduce_level(dim, leaves, max_op, INT_MIN);
}

template int get_max_level(int const, Leaves const&);
template int get_max_level(int const, ZLeaves const&);

template <class LeavesT>
Point get_base_point(int const dim, LeavesT const& leaves)
{
  Vec3<int> ijk(0,0,0);
  int const min_level = get_min_level(dim, leaves);
  for (ID const global_id : leaves) {
    int const leaf_level = get_level(dim, global_id);
    int const diff = leaf_level - min_level;
    Point const leaf_pt = get_point(dim, global_id);
    ijk.x() = std::max(ijk.x(), (leaf_pt.ijk.x() >> diff) + 1);
    ijk.y() = std::max(ijk.y(), (leaf_pt.ijk.y() >> diff) + 1);
    ijk.z() = std::max(ijk.z(), (leaf_pt.ijk.z() >> diff) + 1);
  }
  return Point(min_level, dimensionalize(dim, ijk));
}

template Point get_base_point(int const, Leaves const&);
template Point get_base_point(int const, ZLeaves const&);

struct AdjImpl
{
  AdjacentToLeaf adj;
  bool should_refine = false;
};

static Box3<int> get_grid_bounds(
    int const level,
    Point const& base_pt)
{
  int const dim = infer_dimension(base_pt.ijk);
  int const diff = level - base_pt.level;
  Vec3<int> min = Vec3<int>::zero();
  Vec3<int> max = Vec3<int>::zero();
  for (int axis = 0; axis < dim; ++axis) {
    int const naxis = base_pt.ijk[axis];
    max[axis] = (1 << diff) * naxis - 1;
  }
  return Box3<int>(min, max);
}

static bool is_in(
    int const dim,
    Point const& pt,
    Box3<int> const& bounds)
{
  for (int axis = 0; axis < dim; ++axis) {
    if (pt.ijk[axis] < bounds.lower()[axis]) return false;
    if (pt.ijk[axis] > bounds.upper()[axis]) return false;
  }
  return true;
}

static Point make_periodic(
    Point const& pt,
    Box3<int> const& bounds,
    int const axis,
    int const dir)
{
  Point periodic = pt;
  Vec3<int> const min_ijk = bounds.lower();
  Vec3<int> const max_ijk = bounds.upper();
  if (dir == LEFT) periodic.ijk[axis] = max_ijk[axis];
  if (dir == RIGHT) periodic.ijk[axis] = min_ijk[axis];
  return periodic;
}

static std::vector<Vec3<int>> get_adj_children(
    int const dim,
    int const axis,
    int const dir)
{
  std::vector<Vec3<int>> children;
  auto functor = [&] (Vec3<int> const& child_ijk) {
    if (child_ijk[axis] == dir) return;
    children.push_back(child_ijk);
  };
  seq_for_each(generalize(dim, child_grid), functor);
  return children;
}

static std::vector<ID> get_fine_ids(
    int const dim,
    int const axis,
    int const dir,
    Point const& pt)
{
  std::vector<Vec3<int>> children = get_adj_children(dim, axis, dir);
  std::vector<ID> ids(children.size());
  for (std::size_t i = 0; i < children.size(); ++i) {
    Vec3<int> const child_ijk = children[i];
    Point const fine_pt = get_fine_point(dim, pt, child_ijk);
    ids[i] = get_global_id(dim, fine_pt);
  }
  return ids;
}

static std::vector<ID> get_finer_ids(
    int const dim,
    int const axis,
    int const dir,
    std::vector<ID> const adj_ids)
{
  std::vector<ID> finer_ids;
  for (ID const global_id : adj_ids) {
    Point const pt = get_point(dim, global_id);
    std::vector<ID> const ids = get_fine_ids(dim, axis, dir, pt);
    finer_ids.insert(finer_ids.end(), ids.begin(), ids.end());
  }
  return finer_ids;
}

static bool are_leaves(std::vector<ID> const& ids, Leaves const& leaves)
{
  bool are = true;
  for (ID const id : ids) {
    if (!is_leaf(id, leaves)) are = false;
  }
  return are;
}

static bool is_leaf_in(std::vector<ID> const& ids, Leaves const& leaves)
{
  bool is_in = false;
  for (ID const id : ids) {
    if (is_leaf(id, leaves)) is_in = true;
  }
  return is_in;
}

static AdjImpl get_adj(
    int const dim,
    ID const global_id,
    Leaves const& leaves,
    Point const& base_pt,
    Periodic const& periodic)
{
  AdjImpl result;
  Point const pt = get_point(dim, global_id);
  Box3<int> const bounds = get_grid_bounds(pt.level, base_pt);
  for (int axis = 0; axis < dim; ++axis) {
    for (int dir = 0; dir < DIRECTIONS; ++dir) {
      Vec3<int> const offset = get_dir_sign(dir) * Vec3<int>::axis(axis);
      Point adj_pt(pt.level, pt.ijk + offset);
      if (!is_in(dim, adj_pt, bounds)) {
        if (periodic[axis]) {
          adj_pt = make_periodic(adj_pt, bounds, axis, dir);
        } else {
          continue;
        }
      }
      Adjacent adj;
      adj.axis = static_cast<std::int8_t>(axis);
      adj.dir = static_cast<std::int8_t>(dir);
      ID const adj_id = get_global_id(dim, adj_pt);
      if (is_leaf(adj_id, leaves)) {
        adj.neighbor = adj_id;
        result.adj.push_back(adj);
      } else {
        Point const coarse_adj_pt = get_coarse_point(dim, adj_pt);
        ID const coarse_adj_id = get_global_id(dim, coarse_adj_pt);
        std::vector<ID> const fine_adj_ids = get_fine_ids(dim, axis, dir, adj_pt);
        std::vector<ID> const finer_adj_ids = get_finer_ids(dim, axis, dir, fine_adj_ids);
        bool const is_fine_to_coarse = is_leaf(coarse_adj_id, leaves);
        bool const is_coarse_to_fine = are_leaves(fine_adj_ids, leaves);
        bool const needs_refinement = is_leaf_in(finer_adj_ids, leaves);
        if (is_fine_to_coarse) {
          adj.kind = FINE_TO_COARSE;
          adj.neighbor = coarse_adj_id;
          result.adj.push_back(adj);
        }
        if (is_coarse_to_fine) {
          int const idir = invert_dir(dir);
          auto const my_child_ijks = get_adj_children(dim, axis, idir);
          for (std::size_t i = 0; i < fine_adj_ids.size(); ++i) {
            adj.kind = COARSE_TO_FINE;
            adj.neighbor = fine_adj_ids[i];
            adj.my_child_ijk = my_child_ijks[i];
            result.adj.push_back(adj);
          }
        }
        if (needs_refinement) {
          result.should_refine = true;
        }
        int const invalid_flag =
          int(is_fine_to_coarse) +
          int(is_coarse_to_fine) +
          int(needs_refinement);
        if (invalid_flag > 1) {
          throw std::runtime_error("dgt::tree.get_adj - invalid tree");
        }
      }
    }
  }
  return result;
}

Adjacencies get_adjacencies(
    int const dim,
    ZLeaves const& z_leaves,
    Leaves const& leaves,
    Point const& base_pt,
    Periodic const& periodic)
{
  Adjacencies result;
  result.resize(z_leaves.size());
  for (std::size_t i = 0; i < z_leaves.size(); ++i) {
    ID const global_id = z_leaves[i];
    AdjImpl const impl = get_adj(dim, global_id, leaves, base_pt, periodic);
    result[i] = impl.adj;
    if (impl.should_refine) {
      throw std::runtime_error("dgt::tree::get_adjacencies - invalid tree");
    }
  }
  return result;
}

Box3<real> get_domain(
    int const dim,
    ID const global_id,
    Point const& base_pt,
    Box3<real> const& domain)
{
  Point const leaf_pt = get_point(dim, global_id);
  int const diff = leaf_pt.level - base_pt.level;
  Vec3<real> const length = domain.extents();
  Vec3<real> const base_ijk(base_pt.ijk.x(), base_pt.ijk.y(), base_pt.ijk.z());
  Vec3<real> const leaf_ijk(leaf_pt.ijk.x(), leaf_pt.ijk.y(), leaf_pt.ijk.z());
  Vec3<real> const g_base_ijk = generalize(dim, base_ijk);
  Vec3<real> const dx = comp_division(length, g_base_ijk) / std::pow(2., diff);
  Vec3<real> const min = domain.lower() + comp_product(leaf_ijk, dx);
  Vec3<real> const max = min + dx;
  Vec3<real> const d_min = dimensionalize(dim, min);
  Vec3<real> const d_max = dimensionalize(dim, max);
  return Box3<real>(d_min, d_max);
}

}
}
