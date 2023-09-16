#pragma once

#include "dgt_tree.hpp"
#include "dgt_view.hpp"

namespace dgt {

template <template <class> class ViewT>
struct BlockInfo
{
  ViewT<tree::ID*> global_ids;
  ViewT<std::int8_t*> levels;
  ViewT<Box3<real>*> domains;
  ViewT<Vec3<real>*> dxs;
  ViewT<real*> cell_detJs;
  ViewT<real*> face_detJs[DIMENSIONS];
};

template <template <class> class ViewT>
BlockInfo<ViewT> build_block_info(
    int const dim,
    Box3<real> const& domain,
    tree::OwnedLeaves const& ids,
    tree::Point const& base_pt);

}
