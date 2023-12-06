#include <fstream>
#include <stdexcept>

#include "dgt_cartesian.hpp"
#include "dgt_for_each.hpp"
#include "dgt_mesh.hpp"
#include "dgt_vec3.hpp"
#include "dgt_stl.hpp"

namespace dgt {
namespace stl {

static float parse_float(std::ifstream& file)
{
  char float_buffer[4];
  file.read(float_buffer, 4);
  float* float_ptr = (float*)float_buffer;
  return *float_ptr;
}

static void parse_dummy(std::ifstream& file)
{
  char dummy[2];
  file.read(dummy, 2);
}

HostView<Triangle*> read(std::filesystem::path const& path)
{
  std::ifstream file(path.string().c_str(), std::ios::in | std::ios::binary);
  if (!file.is_open()) {
    throw std::runtime_error("dgt::stl::read - could not open " + path.string());
  }
  char header_char[80] = "";
  char num_triangles_char[4];
  file.read(header_char, 80);
  file.read(num_triangles_char, 4);
  int* num_triangles_ptr = (int*)num_triangles_char;
  int const num_triangles = *num_triangles_ptr;
  HostView<Triangle*> triangles("stl_triangles", num_triangles);
  for (int i = 0; i < num_triangles; ++i) {
    for (int pt = 0; pt < 4; ++pt) {
      triangles[i][pt][X] = parse_float(file);
      triangles[i][pt][Y] = parse_float(file);
      triangles[i][pt][Z] = parse_float(file);
    }
    parse_dummy(file);
  }
  file.close();
  return triangles;
}

DGT_METHOD inline bool ray_intesects_triangle(
    Vec3<real> const& ray_origin,
    Vec3<real> const& ray_vector,
    Triangle const& triangle)
{
  const double EPSILON = 1.e-8;
  Vec3<real> const& v0 = triangle[0];
  Vec3<real> const& v1 = triangle[1];
  Vec3<real> const& v2 = triangle[2];
  Vec3<real> const e1 = v1-v0;
  Vec3<real> const e2 = v2-v0;
  Vec3<real> const rayXe2 = cross(ray_vector, e2);
  real const det = dot(e1, rayXe2);
  if (std::abs(det) < EPSILON) return false;
  real const inv_det = 1./det;
  Vec3<real> const s = ray_origin - v0;
  real const u = inv_det * dot(s, rayXe2);
  if ((u < 0.) || (u > 1.)) return false;
  Vec3<real> const sXe1 = cross(s, e1);
  real const v = inv_det * dot(ray_vector, sXe1);
  if ((v < 0.) || (u+v > 1.)) return false;
  real const t = inv_det * dot(e2, sXe1);
  if (t > EPSILON) return true;
  else return false;
}

View<Triangle*> to_device(HostView<Triangle*> const triangles)
{
  size_t num_triangles = triangles.size();
  View<Triangle*> result("stl_triangles", num_triangles);
  Kokkos::deep_copy(result, triangles);
  return result;
}

Field<real**> compute_vfs(Mesh const& mesh, View<Triangle*> const triangles)
{
  Vec3<real> const ray(0.1, 0.1, 0.1);
  static constexpr int CELL = basis_locations::CELL;
  int const ntriangles = triangles.size();
  Grid3 const cell_grid = mesh.cell_grid();
  int const nblocks = mesh.num_owned_blocks();
  int const ncells = generalize(mesh.dim(), cell_grid).size();
  Basis<View> const& B = mesh.basis();
  BlockInfo<View> const& info = mesh.block_info();
  Field<real**> vfs_field;
  vfs_field.create("stl_vfs", nblocks, ncells, B.num_cell_pts);
  auto vfs = vfs_field.get();
  auto functor = [=] DGT_DEVICE (
      int const block,
      Vec3<int> const& cell_ijk) DGT_ALWAYS_INLINE
  {
    int const cell = cell_grid.index(cell_ijk);
    Vec3<real> const origin = info.domains[block].lower();
    Vec3<real> const dx = info.cell_dxs[block];
    for (int pt = 0; pt < B.num_cell_pts; ++pt) {
      Vec3<real> const xi = get_point(B, CELL, pt);
      Vec3<real> const x = map_to_physical(cell_ijk, origin, dx, xi);
      int num_intersections = 0;
      for (int tri = 0; tri < ntriangles; ++tri) {
        Triangle const& triangle = triangles[tri];
        bool const intersects = ray_intesects_triangle(x, ray, triangle);
        if (intersects) num_intersections++;
      }
      if (num_intersections % 2) vfs[block](cell, pt);
    }
  };
  for_each("stl::compute_vfs", nblocks, cell_grid, functor);
  return vfs_field;
}

}
}
