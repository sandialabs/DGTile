#include "dgt_tree.hpp"

namespace dgt {
namespace tree {

ID get_level_offset(int const dim, int const level)
{
  ID const num = (ID(1) << (dim * level)) - 1;
  ID const den = (ID(1) << dim) - 1;
  return num / den;
}

ID get_level_id(int const dim, Point const& pt)
{
  (void)dim;
  ID const L = pt.level;
  ID const i = pt.ijk.x();
  ID const j = pt.ijk.y();
  ID const k = pt.ijk.z();
  ID const idx1 = j << L;
  ID const idx2 = k << (2*L);
  return i + idx1 + idx2;
}

}
}
