#pragma once

#include <cstdint>
#include <string>
#include <unordered_map>
#include <unordered_set>

#include "dgt_box3.hpp"
#include "dgt_grid3.hpp"

namespace dgt {
namespace tree {

using ID = std::uint64_t;

enum {DEREFINE=-1, REMAIN=0, REFINE=1};
enum {FINE_TO_COARSE=-1, EQUAL=0, COARSE_TO_FINE=1};

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

struct Adjacent
{
  ID neighbor = 0;
  std::int8_t axis = -1;
  std::int8_t dir = -1;
  std::int8_t kind = EQUAL;
  std::int8_t which_child = 0;
  bool boundary = false;
};

using AdjacentToLeaf = std::vector<Adjacent>;
using Adjacencies = std::vector<AdjacentToLeaf>;
using GlobalToZ = std::unordered_map<ID, int>;
using Leaves = std::unordered_set<ID>;
using Levels = std::vector<std::int8_t>;
using Marks = std::vector<std::int8_t>;
using Periodic = Vec3<bool>;
using ZLeaves = std::vector<ID>;

[[nodiscard]] ID get_level_offset(int const dim, int const level);

[[nodiscard]] ID get_level_id(int const dim, Point const& pt);

[[nodiscard]] ID get_global_id(int const dim, Point const& pt);

[[nodiscard]] std::int8_t get_level(int const dim, ID const global_id);

[[nodiscard]] Point get_point(int const dim, ID const global_id);

[[nodiscard]] Point get_coarse_point(int const dim, Point const& pt);

[[nodiscard]] Point get_fine_point(
    int const dim, Point const& pt, Vec3<int> const& child_ijk);

[[nodiscard]] Leaves create(int const dim, Grid3 const& grid);

[[nodiscard]] ZLeaves order(int const dim, Leaves const& leaves);

[[nodiscard]] GlobalToZ invert(int const dim, ZLeaves const& z_leaves);

template <class LeavesT>
[[nodiscard]] int get_max_level(int const dim, LeavesT const& leaves);

template <class LeavesT>
[[nodiscard]] int get_min_level(int const dim, LeavesT const& leaves);

template <class LeavesT>
[[nodiscard]] Point get_base_point(int const dim, LeavesT const& leaves);

[[nodiscard]] Adjacencies get_adjacencies(
    int const dim, ZLeaves const& z_leaves, Leaves const& leaves,
    Point const& base_pt, Periodic const& periodic);

[[nodiscard]] Leaves balance(
    int const dim, Leaves const& leaves,
    Point const& base_pt, Periodic const& periodic);

[[nodiscard]] Leaves modify(
    int const dim, ZLeaves const& z_leaves, Marks const& marks);

[[nodiscard]] Leaves modify(
    int const dim, ZLeaves const& z_leaves, Levels const& levels,
    Point const& base_pt, Periodic const& periodic,
    int const min_level = 0, int const max_level = 16);

[[nodiscard]] Box3<real> get_domain(
    int const dim, ID const global_id,
    Point const& base_pt, Box3<real> const& d);

void write_vtu(
    int const dim, std::string const& prefix,
    ZLeaves const& zl, Box3<real> const& d);

}
}
