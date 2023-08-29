#include <list>

#include "dgt_cartesian.hpp"
#include "dgt_for_each.hpp"
#include "dgt_tree.hpp"

namespace dgt {
namespace tree {

static bool operator==(Point const& a, Point const& b)
{
  return
    (a.level == b.level) &&
    (a.ijk == b.ijk);
}

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

int get_level(int const dim, ID const global_id)
{
  int L = 0;
  while (global_id >= get_level_offset(dim, L+1)) {
    L++;
  }
  return L;
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

static Box3<int> get_grid_bounds(int const level, Point const& base_pt)
{
  int const dim = base_pt.ijk.dimension();
  int const diff = level - base_pt.level;
  Vec3<int> min = Vec3<int>::zero();
  Vec3<int> max = Vec3<int>::zero();
  for (int axis = 0; axis < dim; ++axis) {
    int const naxis = base_pt.ijk[axis];
    max[axis] = (1 << diff)*naxis - 1;
  }
  return Box3<int>(min, max);
}

static bool is_in(
    int const dim,
    Point const& adj_pt,
    Box3<int> const& bounds)
{
  for (int axis = 0; axis < dim; ++axis) {
    if (adj_pt.ijk[axis] < bounds.lower()[axis]) return false;
    if (adj_pt.ijk[axis] > bounds.upper()[axis]) return false;
  }
  return true;
}

static Point make_periodic(
    Point& pt,
    int const dim,
    Periodic const& periodic,
    Vec3<int> const& offset,
    Box3<int> const& bounds)
{
  if (periodic == Vec3<bool>::zero()) return pt;
  Point result = pt;
  for (int axis = 0; axis < dim; ++axis) {
    if (!periodic[axis]) {
      if (pt.ijk[axis] < bounds.lower()[axis]) return pt;
      if (pt.ijk[axis] > bounds.upper()[axis]) return pt;
      continue;
    }
    if (offset[axis] == -1) result.ijk[axis] = bounds.upper()[axis];
    if (offset[axis] ==  1) result.ijk[axis] = bounds.lower()[axis];
  }
  return result;
}

static Vec3<std::int8_t> toi8(Vec3<int> const& ijk)
{
  return Vec3<std::int8_t>(
      std::int8_t(ijk.x()),
      std::int8_t(ijk.y()),
      std::int8_t(ijk.z()));
}

static std::vector<Vec3<int>> get_adj_children(
    int const dim,
    Vec3<int> const& offset)
{
  std::vector<Vec3<int>> children;
  children.push_back(offset);
  for (int axis = 0; axis < dim; ++axis) {
    std::size_t num_children = children.size();
    for (std::size_t i = 0; i < num_children; ++i) {
      Vec3<int>& child = children[i];
      if      (child[axis] == -1) child[axis] = 1;
      else if (child[axis] ==  1) child[axis] = 0;
      else children.push_back(child + Vec3<int>::axis(axis));
    }
  }
  return children;
}

static std::vector<ID> get_fine_ids(
    int const dim,
    Point const& adj_pt,
    Vec3<int> const& offset)
{
  std::vector<Vec3<int>> children = get_adj_children(dim, offset);
  std::vector<ID> ids(children.size());
  for (std::size_t i = 0; i < children.size(); ++i) {
    Vec3<int> const child_ijk = children[i];
    Point const fine_pt = get_fine_point(dim, adj_pt, child_ijk);
    ids[i] = get_global_id(dim, fine_pt);
  }
  return ids;
}

static std::vector<ID> get_finer_ids(
    int const dim,
    std::vector<ID> const& adj_ids,
    Vec3<int> const& offset)
{
  std::vector<ID> finer_ids;
  for (ID const global_id : adj_ids) {
    Point const pt = get_point(dim, global_id);
    std::vector<ID> const ids = get_fine_ids(dim, pt, offset);
    finer_ids.insert(finer_ids.end(), ids.begin(), ids.end());
  }
  return finer_ids;
}

struct AdjImpl
{
  Adjacent adjacent;
  bool should_refine = false;
};

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
  Subgrid3 const grid = dimensionalize(dim, offset_grid);
  auto functor = [&] (Vec3<int> const& offset) {
    if (offset == Vec3<int>::zero()) return;
    Point adj_pt(pt.level, pt.ijk + offset);
    if (!is_in(dim, adj_pt, bounds)) {
      Point const pt = make_periodic(adj_pt, dim, periodic, offset, bounds);
      bool const not_periodic = (pt == adj_pt);
      if (not_periodic) return;
      else adj_pt = pt;
    }
    ID const adj_id = get_global_id(dim, adj_pt);
    if (is_leaf(adj_id, leaves)) {
      result.adjacent.push_back({adj_id, 0, toi8(offset)});
    } else {
      Point const coarse_adj_pt = get_coarse_point(dim, adj_pt);
      ID const coarse_adj_id = get_global_id(dim, coarse_adj_pt);
      std::vector<ID> const fine_adj_ids = get_fine_ids(dim, adj_pt, offset);
      std::vector<ID> const finer_adj_ids = get_finer_ids(dim, fine_adj_ids, offset);
      bool const is_fine_to_coarse = is_leaf(coarse_adj_id, leaves);
      bool const is_coarse_to_fine = are_leaves(fine_adj_ids, leaves);
      bool const needs_refinement = is_leaf_in(finer_adj_ids, leaves);
      if (is_fine_to_coarse) {
        result.adjacent.push_back({coarse_adj_id, -1, toi8(offset)});
      }
      if (is_coarse_to_fine) {
        for (ID const fine_adj_id : fine_adj_ids) {
          result.adjacent.push_back({fine_adj_id, 1, toi8(offset)});
        }
      }
      if (needs_refinement) {
        result.should_refine = true;
      }
      int const sum =
        int(is_fine_to_coarse) +
        int(is_coarse_to_fine) +
        int(needs_refinement);
      if (sum > 1) {
        throw std::runtime_error("dgt::tree.get_adj - invalid tree");
      }
    }
  };
  seq_for_each(grid, functor);
  return result;
}

Adjacencies get_adjacencies(
    int const dim,
    Leaves const& leaves,
    Point const& base_pt,
    Periodic const& periodic)
{
  Adjacencies result;
  for (ID const global_id : leaves) {
    AdjImpl const impl = get_adj(dim, global_id, leaves, base_pt, periodic);
    result[global_id] = impl.adjacent;
    if (impl.should_refine) {
      throw std::runtime_error("dgt::tree:get_adjacencies - invalid tree");
    }
  }
  return result;
}

static void insert_children(
    Leaves& leaves,
    int const dim,
    ID const global_id)
{
  Point const pt = get_point(dim, global_id);
  auto functor = [&] (Vec3<int> const& child_ijk) {
    Point const fine_pt = get_fine_point(dim, pt, child_ijk);
    ID const fine_global_id = get_global_id(dim, fine_pt);
    leaves.insert(fine_global_id);
  };
  seq_for_each(dimensionalize(dim, child_grid), functor);
}

static void insert_parent(
    Leaves& leaves,
    int const dim,
    ID const global_id)
{
  Point const pt = get_point(dim, global_id);
  Point const coarse_pt = get_coarse_point(dim, pt);
  ID const coarse_global_id = get_global_id(dim, coarse_pt);
  leaves.insert(coarse_global_id);
}

static std::vector<ID> get_refines(
    int const dim,
    Leaves const& leaves,
    Point const& base_pt,
    Periodic const& periodic)
{
  std::vector<ID> result;
  for (ID const global_id : leaves) {
    AdjImpl const impl = get_adj(dim, global_id, leaves, base_pt, periodic);
    if (impl.should_refine) {
      result.push_back(global_id);
    }
  }
  return result;
}

Leaves balance(
    int const dim,
    Leaves const& leaves,
    Point const& base_pt,
    Periodic const& periodic)
{
  Leaves result = leaves;
  bool needs_modification = true;
  while (needs_modification) {
    ZLeaves const z_leaves = order(dim, result);
    std::vector<ID> const ids = get_refines(dim, result, base_pt, periodic);
    std::size_t const num_refines = ids.size();
    if (num_refines > 0) {
      for (std::size_t i = 0; i < num_refines; ++i) {
        ID const leaf = ids[i];
        insert_children(result, dim, leaf);
        result.erase(leaf);
      }
    } else {
      needs_modification = false;
    }
  }
  return result;
}

static std::unordered_set<ID> get_valid_parents(
    int const dim,
    ZLeaves const& z_leaves,
    Marks const& marks)
{
  std::list<ID> sorted_parents;
  std::list<ID> unique_parents;
  std::unordered_set<ID> result;
  for (std::size_t i = 0; i < z_leaves.size(); ++i) {
    if (marks[i] == DEREFINE) {
      ID const leaf = z_leaves[i];
      Point const pt = get_point(dim, leaf);
      Point const coarse_pt = get_coarse_point(dim, pt);
      ID const parent = get_global_id(dim, coarse_pt);
      unique_parents.push_back(parent);
    }
  }
  unique_parents.sort();
  sorted_parents = unique_parents;
  unique_parents.unique();
  int const num_child = dimensionalize(dim, child_grid).size();
  for (ID const parent : unique_parents) {
    int const count = std::count(
        sorted_parents.begin(), sorted_parents.end(), parent);
    if (count == num_child) result.insert(parent);
  }
  return result;
}

static bool can_derefine(
    int const dim,
    ID const leaf,
    std::unordered_set<ID> const& parents)
{
  Point const pt = get_point(dim, leaf);
  Point const coarse_pt = get_coarse_point(dim, pt);
  ID const parent = get_global_id(dim, coarse_pt);
  if (parents.count(parent)) return true;
  else return false;
}

Leaves modify(
    int const dim,
    ZLeaves const& z_leaves,
    Marks const& marks)
{
  Leaves leaves;
  auto const parents = get_valid_parents(dim, z_leaves, marks);
  for (std::size_t i = 0; i < z_leaves.size(); ++i) {
    ID const leaf = z_leaves[i];
    if (marks[i] == REMAIN) leaves.insert(leaf);
    if (marks[i] == REFINE) insert_children(leaves, dim, leaf);
    if (marks[i] == DEREFINE) {
      if (can_derefine(dim, leaf, parents)) {
        insert_parent(leaves, dim, leaf);
      } else {
        leaves.insert(leaf);
      }
    }
  }
  return leaves;
}

Leaves modify(
    int const dim,
    ZLeaves const& z_leaves,
    Levels const& levels,
    Point const& base_pt,
    Periodic const& periodic,
    int const min_level,
    int const max_level)
{
  Marks marks(z_leaves.size(), REMAIN);
  for (std::size_t i = 0; i < z_leaves.size(); ++i) {
    ID const leaf = z_leaves[i];
    int const current_level = get_level(dim, leaf);
    int const desired_level = int(levels[i]);
    bool const can_refine = (current_level < max_level);
    bool const can_derefine = (current_level > min_level);
    bool const should_refine = (desired_level > current_level);
    bool const should_derefine = (desired_level < current_level);
    if (should_refine && can_refine) marks[i] = REFINE;
    if (should_derefine && can_derefine) marks[i] = DEREFINE;
  }
  Leaves leaves = modify(dim, z_leaves, marks);
  leaves = balance(dim, leaves, base_pt, periodic);
  return leaves;
}

}
}
