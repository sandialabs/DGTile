#include <algorithm>

#include "dgt_partitioning.hpp"

namespace dgt {
namespace linear_partitioning {

int get_num_local(
    int const num_total,
    int const num_parts,
    int const which_part)
{
  int const quotient = num_total / num_parts;
  int const remainder = num_total % num_parts;
  return (which_part < remainder) ? (quotient + 1) : quotient;
}

int get_local_offset(
    int const num_total,
    int const num_parts,
    int const which_part)
{
  int const quotient = num_total / num_parts;
  int const remainder = num_total % num_parts;
  return (quotient * which_part) + std::min(remainder, which_part);
}

tree::ZLeaves get_owned_leaves(
    int const rank,
    int const num_ranks,
    tree::ZLeaves const& z_leaves)
{
  int const nleaves = int(z_leaves.size());
  int const nlocal = get_num_local(nleaves, num_ranks, rank);
  int const offset = get_local_offset(nleaves, num_ranks, rank);
  auto begin = z_leaves.begin() + offset;
  auto end = z_leaves.begin() + offset + nlocal;
  return std::vector<tree::ID>(begin, end);
}

}
}
