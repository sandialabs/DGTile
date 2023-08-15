#pragma once

#include "dgt_when.hpp"

namespace dgt {

namespace lua {
class table;
}

WhenPtr make_when(lua::table const& in);

}
