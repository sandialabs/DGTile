#pragma once

namespace dgt {
namespace linear_partitioning {

[[nodiscard]] int get_num_local(
    int const num_total,
    int const num_parts,
    int const which_part);

[[nodiscard]] int get_local_offset(
    int const num_total,
    int const num_parts,
    int const which_part);

}
}
