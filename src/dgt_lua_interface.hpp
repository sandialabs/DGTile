#pragma once

#include "dgt_dg.hpp"
#include "dgt_view.hpp"
#include "dgt_when.hpp"

namespace dgt {

namespace lua {
class table;
}

WhenPtr make_when(lua::table const& in);
Basis<View> make_basis(lua::table const& in);

}
