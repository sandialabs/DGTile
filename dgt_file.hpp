#pragma once

#include <filesystem>
#include <sstream>
#include <string>
#include <stdexcept>

#include "dgt_defines.hpp"
#include "dgt_views.hpp"

namespace dgt {

class Block;
class Mesh;
class Tree;

template <class T>
T string_to_type(std::string const& s);

template <>
inline int string_to_type(std::string const& s) {
  return std::stoi(s);
}

template <>
inline double string_to_type(std::string const& s) {
  return std::stod(s);
}

template <>
inline bool string_to_type(std::string const& s) {
  if (s == "true") return true;
  if (s == "false") return false;
  throw std::runtime_error("to_type<bool>");
}

void write_stream(
    std::filesystem::path const& path,
    std::stringstream const& stream);

namespace base64 {
std::size_t encoded_size(std::size_t size);
std::string encode(void const* data, std::size_t size);
void decode(std::string const& text, void* data, std::size_t size);
}

namespace ascii {
void write_mesh(std::filesystem::path const& path, Mesh const& mesh);
void read_mesh(std::filesystem::path const& path, Mesh& mesh);
}

namespace binary {
void write_mesh(std::filesystem::path const& path, Mesh const& mesh);
void read_mesh(std::filesystem::path const& path, Mesh& mesh);
}

namespace vtk {

template <class T>
using VizView = Kokkos::DualView<T**, Kokkos::LayoutRight>;

void write_vtr_start(
    std::stringstream& stream, Block const& block, double time, int step);
template <class T>
void write_field(std::stringstream& stream, std::string const& name, VizView<T> f);
void write_vtr_end(std::stringstream& stream);
void write_vtm(std::stringstream& stream, std::string const& prefix, int nblocks);

void write_tree(
    std::filesystem::path const& path,
    Tree& tree,
    p3a::box3<double> const& domain);

}

}
