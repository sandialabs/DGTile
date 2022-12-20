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

void write_tree(std::filesystem::path const& path, Mesh const& mesh);

void write_pvtu_start(std::stringstream& stream, int nblocks);
void write_pvtu_point_data_start(std::stringstream& stream);
void write_pvtu_point_data_end(std::stringstream& stream);
void write_pvtu_cell_data_start(std::stringstream& stream);
void write_pvtu_cell_data_end(std::stringstream& stream);
void write_pvtu_end(std::stringstream& stream);

void write_vtu_start(std::stringstream& stream, Block const& block);
void write_vtu_point_data_start(std::stringstream& stream);
void write_vtu_point_data_end(std::stringstream& stream);
void write_vtu_cell_data_start(std::stringstream& stream);
void write_vtu_cell_data_end(std::stringstream& stream);
void write_vtu_end(std::stringstream& stream);

template <class T>
void write_vtu_cell_field(
    std::stringstream& stream, std::string const& name, VizView<T> f);

template <class T>
void write_pvtu_cell_field(
    std::stringstream& stream, std::string const& name, int ncomps);

}

}
