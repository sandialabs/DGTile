#pragma once

#include <filesystem>
#include <vector>

#include "dgt_defines.hpp"
#include "dgt_tree.hpp"
#include "dgt_vec.hpp"

namespace dgt {
namespace stl {

using Triangle = Vec<Vec3<real>, 4>;
using Triangles = std::vector<Triangle>;

Triangles read(std::filesystem::path const& path);

}
}
