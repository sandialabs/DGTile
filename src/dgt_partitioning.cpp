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

PartInfo get_part_info(
    int const num_parts,
    tree::ID const global_id,
    tree::GlobalToZ const& inv_z_leaves)
{
  PartInfo result;
  int const num_total = int(inv_z_leaves.size());
  int const zid = inv_z_leaves.at(global_id);
  int const quotient = num_total / num_parts;
  int const remainder = num_total % num_parts;
  int const transition_zid = get_local_offset(num_total, num_parts, remainder);
  if (zid < transition_zid) {
    result.rank = zid / (quotient + 1);
    result.block = zid - result.rank * (quotient + 1);
  } else {
    result.rank = remainder + (zid - transition_zid) / quotient;
    result.block = (zid - transition_zid) - (result.rank - remainder) * quotient;
  }
  return result;
}

}
}
