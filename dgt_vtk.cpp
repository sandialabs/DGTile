#include <fstream>

#include "zlib.h"

#include "caliper/cali.h"

#include "dgt_basis.hpp"
#include "dgt_file.hpp"
#include "dgt_mesh.hpp"
#include "dgt_grid.hpp"
#include "dgt_spatial.hpp"

namespace dgt {

namespace vtk {

template <class T> std::string vtk_type_name();
template <> std::string vtk_type_name<int>() { return "Int32"; }
template <> std::string vtk_type_name<float>() { return "Float32"; }
template <> std::string vtk_type_name<double>() { return "Float64"; }

static void write_vtr_header(std::stringstream& stream) {
  stream << "<VTKFile type=\"RectilinearGrid\" ";
  stream << "version=\"1.0\" ";
  stream << "compressor=\"vtkZLibDataCompressor\" ";
  stream << "header_type=\"UInt64\">\n";
}

static void write_vtr_rectilinear_start(
    std::stringstream& stream,
    Block const& block) {
  int const p = block.basis().p;
  p3a::grid3 const cells = block.cell_grid();
  p3a::vector3<int> const n = (p+1)*cells.extents();
  stream << "<RectilinearGrid WholeExtent=\"";
  stream << "0 " << n.x() << " ";
  stream << "0 " << n.y() << " ";
  stream << "0 " << n.z() << "\">\n";
}

static void write_fdata_start(
    std::stringstream& stream,
    std::string const& type,
    std::string const& name,
    int ntuples) {
  stream << "<DataArray ";
  stream << "type=\"" << type << "\"";
  stream << " Name=\"" << name << "\"";
  stream << " NumberOfTuples=\"" << ntuples << "\"";
  stream << " format=\"ascii\"";
  stream << ">\n";
}

static void write_data_start(
    std::stringstream& stream,
    std::string const& type,
    std::string const& name,
    int ncomps) {
  stream << "<";
  stream << "DataArray type=\"" << type << "\" ";
  stream << "Name=\"" << name << "\" ";
  stream << "NumberOfComponents=\"" << ncomps << "\" ";
  stream << "format=\"binary\"";
  stream << ">\n";
}

static void write_data_end(std::stringstream& stream) {
  stream << "</DataArray>\n";
}

static void write_vtr_time(std::stringstream& stream, double time) {
  write_fdata_start(stream, "Float64", "TIME", 1);
  stream << time << "\n";
  write_data_end(stream);
}

static void write_vtr_step(std::stringstream& stream, int step) {
  write_fdata_start(stream, "Int32", "STEP", 1);
  stream << step << "\n";
  write_data_end(stream);
}

static void write_vtr_block_depth(std::stringstream& stream, Block const& block) {
  Node const* node = block.node();
  Point const pt = node->pt();
  int const depth = pt.depth;
  write_fdata_start(stream, "Int32", "block_depth", 1);
  stream << depth << "\n";
  write_data_end(stream);
}

static void write_vtr_block_ijk(std::stringstream& stream, Block const& block) {
  Node const* node = block.node();
  Point const pt = node->pt();
  p3a::vector3<int> const ijk = pt.ijk;
  write_fdata_start(stream, "Int32", "block_ijk", DIMS);
  stream << ijk.x() << " " << ijk.y() << " " << ijk.z() << "\n";
  write_data_end(stream);
}

static void write_vtr_block_id(std::stringstream& stream, Block const& block) {
  int const id = block.id();
  write_fdata_start(stream, "Int32", "block_id", 1);
  stream << id << "\n";
  write_data_end(stream);
}

static void write_vtr_block_owner(std::stringstream& stream, Block const& block) {
  int const owner = block.owner();
  write_fdata_start(stream, "Int32", "block_owner", 1);
  stream << owner << "\n";
  write_data_end(stream);
}

static void write_vtr_field_data(
    std::stringstream& stream,
    Block const& block,
    double time,
    int step) {
  stream << std::scientific << std::setprecision(12);
  stream << "<FieldData>\n";
  write_vtr_time(stream, time);
  write_vtr_step(stream, step);
  write_vtr_block_depth(stream, block);
  write_vtr_block_ijk(stream, block);
  write_vtr_block_id(stream, block);
  write_vtr_block_owner(stream, block);
  stream << "</FieldData>\n";
}

static void write_piece_start(std::stringstream& stream, Block const& block) {
  int const p = block.basis().p;
  p3a::grid3 const cells = block.cell_grid();
  p3a::vector3<int> const n = (p+1)*cells.extents();
  stream << "<Piece Extent=\"";
  stream << "0 " << n.x() << " ";
  stream << "0 " << n.y() << " ";
  stream << "0 " << n.z() << "\">\n";
}

template <class T>
void write_data(
    std::stringstream& stream,
    VizView<T> dual,
    bool copy = true) {
  if (copy) dual.template sync<typename VizView<T>::host_mirror_space>();
  auto field = dual.h_view;
  std::uint64_t uncompressed_bytes = sizeof(T) * static_cast<uint64_t>(field.size());
  uLong source_bytes = uncompressed_bytes;
  uLong dest_bytes = ::compressBound(source_bytes);
  auto compressed = new ::Bytef[dest_bytes];
  int ret = ::compress2(compressed, &dest_bytes,
      reinterpret_cast<const ::Bytef*>(field.data()),
      source_bytes, Z_BEST_SPEED);
  if (ret != Z_OK) throw std::runtime_error("vtk - zlib error");
  std::string const encoded = base64::encode(compressed, dest_bytes);
  delete[] compressed;
  std::uint64_t header[4] = {1, uncompressed_bytes, uncompressed_bytes, dest_bytes};
  std::string const enc_header = base64::encode(header, sizeof(header));
  stream.write(enc_header.data(), std::streamsize(enc_header.length()));
  stream.write(encoded.data(), std::streamsize(encoded.length()));
  stream.write("\n", 1);
}

static void write_coordinate(std::stringstream& stream, Block const& block, int axis) {
  std::string const axis_name[DIMS] = {"x", "y", "z"};
  int const p = block.basis().p;
  int const num_pts = (p+1)*block.cell_grid().extents()[axis] + 1;
  double const o = block.domain().lower()[axis];
  double const dx = block.dx()[axis] / (p+1);
  VizView<float> coord;
  Kokkos::resize(coord, num_pts, 1);
  double const offset[max_p+1][max_p+1] = {
    {0.,  0.,    0.},
    {0.,  0.,    0.},
    {0., -2./9., 2./9.}
  };
  for (int i = 0; i < num_pts; ++i) {
    int const mod = i%(p+1);
    coord.h_view(i, 0) = o + i*dx + offset[p][mod]*dx;
  }
  write_data_start(stream, "Float32", axis_name[axis], 1);
  write_data(stream, coord, false);
  write_data_end(stream);
}

static void write_coordinates(std::stringstream& stream, Block const& block) {
  stream << "<Coordinates>\n";
  write_coordinate(stream, block, X);
  write_coordinate(stream, block, Y);
  write_coordinate(stream, block, Z);
  stream << "</Coordinates>\n";
}

void write_vtr_start(
    std::stringstream& stream,
    Block const& block,
    double time,
    int step) {
  write_vtr_header(stream);
  write_vtr_rectilinear_start(stream, block);
  write_vtr_field_data(stream, block, time, step);
  write_piece_start(stream, block);
  write_coordinates(stream, block);
  stream << "<CellData>\n";
}

template <class T>
void write_field(
    std::stringstream& stream,
    std::string const& name,
    VizView<T> f) {
  CALI_CXX_MARK_FUNCTION;
  int const ncomps = f.d_view.extent(1);
  write_data_start(stream, vtk_type_name<T>(), name, ncomps);
  write_data(stream, f, true);
  write_data_end(stream);
}

template void write_field<int>(std::stringstream&, std::string const&, VizView<int> f);
template void write_field<float>(std::stringstream&, std::string const&, VizView<float> f);
template void write_field<double>(std::stringstream&, std::string const&, VizView<double> f);

void write_vtr_end(std::stringstream& stream) {
  stream << "</CellData>\n";
  stream << "</Piece>\n";
  stream << "</RectilinearGrid>\n";
  stream << "</VTKFile>\n";
}

static void write_vtm_header(std::stringstream& stream) {
  stream << "<VTKFile type=\"vtkMultiBlockDataSet\" ";
  stream << "version=\"1.0\">\n";
  stream << "<vtkMultiBlockDataSet>\n";
}

static void write_vtm_end(std::stringstream& stream) {
  stream << "</vtkMultiBlockDataSet>\n";
  stream << "</VTKFile>";
}

static void write_vtm_source_file(
    std::stringstream& stream,
    int i,
    std::string const& file) {
  stream << "<DataSet index=\"" << i << "\" ";
  stream << "file=\"" << file << "\"/>\n";
}

void write_vtm(std::stringstream& stream, std::string const& prefix, int nblocks) {
  write_vtm_header(stream);
  for (int i = 0; i < nblocks; ++i) {
    std::string const file = prefix + std::to_string(i) + ".vtr";
    write_vtm_source_file(stream, i, file);
  }
  write_vtm_end(stream);
}

}

}
