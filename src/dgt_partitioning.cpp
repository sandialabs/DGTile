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

}
}
