#pragma once

#include <cstdint>
#include <unordered_map>
#include <unordered_set>

#include "dgt_grid3.hpp"

namespace dgt {
namespace tree {

using ID = std::uint64_t;

struct Point
{
  int level = -1;
  Vec3<int> ijk = Vec3<int>::zero();
  DGT_ALWAYS_INLINE constexpr Point(int l, Vec3<int> const& loc)
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

DGT_METHOD ID get_level_offset(int const dim, int const level);
DGT_METHOD ID get_level_id(int const dim, Point const& pt);

}
}
