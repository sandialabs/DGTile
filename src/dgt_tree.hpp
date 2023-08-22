#pragma once

#include <cstdint>
#include <unordered_map>
#include <unordered_set>

#include "dgt_box3.hpp"
#include "dgt_grid3.hpp"

namespace dgt {
namespace tree {

using ID = std::uint64_t;

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

struct Adjacency
{
  ID neighbor = -1;
  std::int8_t level_offset = -1;
  Vec3<std::int8_t> ijk_offset = Vec3<std::int8_t>::zero();
};

using Adjacent = std::vector<Adjacency>;
using Adjacencies = std::unordered_map<ID, Adjacent>;
using Leaves = std::unordered_set<ID>;
using Levels = std::vector<std::int8_t>;
using Periodic = Vec3<bool>;
using ZLeaves = std::vector<ID>;

[[nodiscard]] ID get_level_offset(int const dim, int const level);

[[nodiscard]] ID get_level_id(int const dim, Point const& pt);

[[nodiscard]] ID get_global_id(int const dim, Point const& pt);

[[nodiscard]] int get_level(int const dim, ID const global_id);

[[nodiscard]] Point get_point(int const dim, ID const global_id);

[[nodiscard]] Point get_coarse_point(int const dim, ID const global_id);

[[nodiscard]] Point get_fine_point(
    int const dim, Point const& pt, Vec3<int> const& child_ijk);

[[nodiscard]] Leaves create(int const dim, Grid3 const& grid);

[[nodiscard]] ZLeaves order(int const dim, Leaves const& leaves);

template <class LeavesT>
[[nodiscard]] int get_max_level(int const dim, LeavesT const& leaves);

template <class LeavesT>
[[nodiscard]] int get_min_level(int const dim, LeavesT const& leaves);

template <class LeavesT>
[[nodiscard]] Point get_base_point(int const dim, LeavesT const& leaves);

[[nodiscard]] Box3<real> get_domain(
    int const dim, ID const global_id,
    Point const& base, Box3<real> const& d);

[[nodiscard]] Adjacencies get_adjacencies(
    int const dim, Leaves const& leaves,
    Point const& base_pt, Periodic const& periodic);

void write_vtu(
    int const dim, std::string const& prefix,
    ZLeaves const& zl, Box3<real> const& d);

}
}
