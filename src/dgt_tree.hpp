#pragma once

#include <cstdint>
#include <unordered_map>
#include <unordered_set>

#include "dgt_grid3.hpp"

namespace dgt {
namespace tree {

using ID = std::uint64_t;

enum {DEREFINE = -1, REMAIN = 0, REFINE = 1};

struct Point
{
  int level = -1;
  Vec3<int> ijk = Vec3<int>::zero();
  DGT_METHOD constexpr Point(int l, Vec3<int> const& loc)
    :level(l)
    ,ijk(loc)
  {
  }
};

struct Adjacency {
  ID neighbor = -1;
  std::int8_t level_difference = -1;
  Vec3<std::int8_t> ijk = Vec3<std::int8_t>::zero();
};

using Leaves = std::unordered_set<ID>;
using ZLeaves = std::vector<ID>;
using Adjacent = std::vector<Adjacency>;
using Adjacencies = std::unordered_map<ID, Adjacent>;
using Marks = std::vector<std::int8_t>;

static constexpr Grid3 child_grid = {2,2,2};
static constexpr Grid3 adj_grid = {3,3,3};
static constexpr Grid3 fine_adj_grid = {4,4,4};

[[nodiscard]] ID get_level_offset(int const dim, int const level);
[[nodiscard]] ID get_level_id(int const dim, Point const& pt);
[[nodiscard]] ID get_global_id(int const dim, Point const& pt);
[[nodiscard]] int get_level(int const dim, ID const global_id);
[[nodiscard]] Point get_point(int const dim, ID const global_id);
[[nodiscard]] Point get_coarse_point(int const dim, ID const global_id);
[[nodiscard]] Point get_fine_point(int const dim, Point const& pt, Vec3<int> const& child_ijk);
[[nodiscard]] Leaves create(int const dim, Grid3 const& grid);
[[nodiscard]] ZLeaves order(int const dim, Leaves const& leaves);
[[nodiscard]] Leaves modify(int const dim, ZLeaves const& z_leaves, Marks  const& marks);
template <class LeavesT> [[nodiscard]] int get_max_level(int const dim, LeavesT const& leaves);
template <class LeavesT> [[nodiscard]] int get_min_level(int const dim, LeavesT const& leaves);
template <class LeavesT> [[nodiscard]] Point get_base_point(int const dim, LeavesT const& leaves);

// get_domain
// get_adjacencies
// check_marks
// balance_tree

}
}