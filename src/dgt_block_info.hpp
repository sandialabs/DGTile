#pragma once

#include "dgt_tree.hpp"
#include "dgt_view.hpp"

namespace dgt {

struct BlockInfo
{
  View<tree::ID*> global_ids;
  View<std::int8_t*> levels;
  View<Box3<real>*> domains;
  View<Vec3<real>*> dxs;
  View<real*> cell_detJs;
  View<real*> face_detJs[DIMENSIONS];
};

BlockInfo create_block_info(
    int const dim,
    Box3<real> const& domain,
    tree::Leaves const& ids,
    tree::Point const& base_pt);

}
