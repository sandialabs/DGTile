#pragma once

#include <filesystem>
#include <vector>

#include "dgt_field.hpp"
#include "dgt_vec.hpp"
#include "dgt_vec3.hpp"
#include "dgt_view.hpp"

namespace dgt {

class Mesh;

namespace stl {

using Triangle = Vec<Vec3<real>, 4>;

HostView<Triangle*> read(std::filesystem::path const& path);
View<Triangle*> to_device(HostView<Triangle*> const triangles);
Field<real**> compute_vfs(Mesh const& mesh, View<Triangle*> const triangles);

}
}
