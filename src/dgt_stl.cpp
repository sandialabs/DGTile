#include <fstream>
#include <stdexcept>

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

Triangles read(std::filesystem::path const& path)
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
  Triangles triangles(num_triangles);
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

}
}
