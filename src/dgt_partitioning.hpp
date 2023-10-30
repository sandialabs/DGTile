#pragma once

#include "dgt_tree.hpp"

namespace dgt {

struct PartInfo
{
  int rank = -1;
  int block = -1;
};

namespace linear_partitioning {

[[nodiscard]] int get_num_local(
    int const num_total,
    int const num_parts,
    int const which_part);

[[nodiscard]] int get_local_offset(
    int const num_total,
    int const num_parts,
    int const which_part);

[[nodiscard]] tree::ZLeaves get_owned_leaves(
    int const rank,
    int const num_ranks,
    tree::ZLeaves const& z_leaves);

[[nodiscard]] PartInfo get_part_info(
    int const num_ranks,
    tree::ID const global_id,
    tree::GlobalToZ const& inv_z_leaves);

}
}
