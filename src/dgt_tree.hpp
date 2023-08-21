#pragma once

#include <cstdint>
#include <unordered_map>
#include <unordered_set>

#include "dgt_box3.hpp"
#include "dgt_grid3.hpp"

namespace dgt {
namespace tree {

using ID = std::uint64_t;

enum {DEREFINE = -1, REMAIN = 0, REFINE = 1};
enum {COARSE_TO_FINE = -1, EQUAL = 0, FINE_TO_COARSE = 1};

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
  std::int8_t level_offset = -1;
  Vec3<std::int8_t> ijk_offset = Vec3<std::int8_t>::zero();
};

using Leaves = std::unordered_set<ID>;
using ZLeaves = std::vector<ID>;
using Boundaries = std::vector<std::vector<std::vector<ID>>>;
using Adjacent = std::vector<Adjacency>;
using Adjacencies = std::unordered_map<ID, Adjacent>;
using Marks = std::vector<std::int8_t>;

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
[[nodiscard]] Box3<real> get_domain(int const dim, ID const gid, Point const& base, Box3<real> const& d);
[[nodiscard]] Boundaries get_boundaries(int const dim, Leaves const& leaves, Point const& base);
[[nodiscard]] Adjacencies get_adjacencies(int const dim, Leaves const& leaves);
void write_vtu(int const dim, std::string const& prefix, ZLeaves const& zl, Box3<real> const& d);

// TODO:
// get_boundaries
// add periodic to adjacencies
// check_marks
// balance_tree

}
}
