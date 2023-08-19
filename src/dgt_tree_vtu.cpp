#include <fstream>
#include <filesystem>
#include <stdexcept>
#include <sstream>

#include "dgt_tree.hpp"

namespace dgt {
namespace tree {

static void write_vtu_header(std::stringstream& out)
{
  out << "<VTKFile type=\"UnstructuredGrid\" ";
  out << "header_type=\"UInt64\">\n";
}

static void write_data_start(
    std::stringstream& stream,
    std::string const& type,
    std::string const& name,
    int ncomps)
{
  stream << "<";
  stream << "DataArray type=\"" << type << "\" ";
  stream << "Name=\"" << name << "\" ";
  stream << "NumberOfComponents=\"" << ncomps << "\" ";
  stream << "format=\"ascii\"";
  stream << ">\n";
}

static void write_data_end(std::stringstream& stream) {
  stream << "</DataArray>\n";
}

static void write_types(
    int const dim,
    std::stringstream& stream,
    std::size_t const nleaves)
{
  static constexpr int vtk_types[] = {1,3,9,12};
  int const type = vtk_types[dim];
  write_data_start(stream, "Int8", "types", 1);
  for (std::size_t i = 0; i < nleaves; ++i) {
    stream << type << "\n";
  }
  write_data_end(stream);
}

static void write_offsets(
    std::stringstream& stream,
    std::size_t const nleaves,
    std::size_t const ncorners)
{
  write_data_start(stream, "Int32", "offsets", 1);
  std::size_t offset = ncorners;
  for (std::size_t i = 0; i < nleaves; ++i) {
    stream << offset << "\n";
    offset += ncorners;
  }
  write_data_end(stream);
}

static void write_connectivity(
    std::stringstream& stream,
    std::size_t const npoints)
{
  write_data_start(stream, "Int32", "connectivity", 1);
  for (std::size_t i = 0; i < npoints; ++i) {
    stream << i << "\n";
  }
  write_data_end(stream);
}

static void write_coordinates(
    std::stringstream& stream,
    int const dim,
    int const ncorners,
    ZLeaves const& zlvs,
    Box3<real> const& domain)
{
  int constexpr vtk_corners[8][3] = {
    {0,0,0},
    {1,0,0},
    {1,1,0},
    {0,1,0},
    {0,0,1},
    {1,0,1},
    {1,1,1},
    {0,1,1},
  };
  write_data_start(stream, "Float64", "coordinates", 3);
  Point const base_pt = get_base_point(dim, zlvs);
  for (ID const global_id : zlvs) {
    Box3<real> const leaf_domain = get_domain(dim, global_id, base_pt, domain);
    Vec3<double> const dx = leaf_domain.extents();
    Vec3<double> const o = leaf_domain.lower();
    for (int c = 0; c < ncorners; ++c) {
      double const x = o.x() + dx.x() * vtk_corners[c][0];
      double const y = o.y() + dx.y() * vtk_corners[c][1];
      double const z = o.z() + dx.z() * vtk_corners[c][2];
      stream << x << " " << y << " " << z << "\n";
    }
  }
  write_data_end(stream);
}

static void write_zcurve_ids(
    std::stringstream& stream,
    ZLeaves const& zlvs)
{
  write_data_start(stream, "Int32", "zcurve_id", 1);
  for (std::size_t i = 0; i < zlvs.size(); ++i) {
    stream << i << "\n";
  }
  write_data_end(stream);
}

static void write_global_ids(
    std::stringstream& stream,
    ZLeaves const& zlvs)
{
  write_data_start(stream, "Int64", "global_id", 1);
  for (ID const global_id : zlvs) {
    stream << global_id << "\n";
  }
  write_data_end(stream);
}

static void write_levels(
    std::stringstream& stream,
    int const dim,
    ZLeaves const& zlvs)
{
  write_data_start(stream, "Int32", "level", 1);
  for (ID const global_id : zlvs) {
    Point const pt = get_point(dim, global_id);
    stream << pt.level << "\n";
  }
  write_data_end(stream);
}

static void write_ijks(
    std::stringstream& stream,
    int const dim,
    ZLeaves const& zlvs)
{
  write_data_start(stream, "Int32", "ijk", 3);
  for (ID const global_id : zlvs) {
    Point const pt = get_point(dim, global_id);
    stream << pt.ijk.x() << " " << pt.ijk.y() << " " << pt.ijk.z() << "\n";
  }
  write_data_end(stream);
}

static void write_stream(
    std::filesystem::path const& path,
    std::stringstream const& stream)
{
  std::ofstream file(path.c_str());
  if (!file.is_open()) {
    throw std::runtime_error("intree: could not open: " + path.string());
  }
  file << stream.rdbuf();
  file.close();
}

void write_vtu(
    int const dim,
    std::string const& path,
    ZLeaves const& zlvs,
    Box3<real> const& dom)
{
  std::stringstream stream;
  write_vtu_header(stream);
  std::size_t nleaves = zlvs.size();
  std::size_t ncorners = 1 << dim;
  std::size_t npoints = ncorners*nleaves;
  stream << "<UnstructuredGrid>\n";
  stream << "<Piece NumberOfPoints=\"" << npoints << "\" ";
  stream << "NumberOfCells=\"" << nleaves << "\">\n";
  stream << "<Cells>\n";
  write_types(dim, stream, nleaves);
  write_offsets(stream, nleaves, ncorners);
  write_connectivity(stream, npoints);
  stream << "</Cells>\n";
  stream << "<Points>\n";
  write_coordinates(stream, dim, ncorners, zlvs, dom);
  stream << "</Points>\n";
  stream << "<CellData>\n";
  write_zcurve_ids(stream, zlvs);
  write_global_ids(stream, zlvs);
  write_levels(stream, dim, zlvs);
  write_ijks(stream, dim, zlvs);
  stream << "</CellData>\n";
  stream << "<PointData>\n";
  stream << "</PointData>\n";
  stream << "</Piece>\n";
  stream << "</UnstructuredGrid>\n";
  stream << "</VTKFile>\n";
  std::filesystem::path const file_path(path + ".vtu");
  write_stream(file_path, stream);
}

}
}
