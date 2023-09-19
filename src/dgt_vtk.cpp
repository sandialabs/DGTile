#include <iomanip>

#include <zlib.h>

#include "dgt_base64.hpp"
#include "dgt_cartesian.hpp"
#include "dgt_mesh.hpp"
#include "dgt_vtk.hpp"

namespace dgt {
namespace vtk {

static void check_q(int const q)
{
  if (q >= max_1D_quadrature_points) {
    std::string const msg = fmt::format(
        "dgt::vtk-> invalid quadrature rule `{}`", q);
    throw std::runtime_error(msg);
  }
}

static void write_vtr_header(std::stringstream& stream)
{
  stream << "<VTKFile type=\"RectilinearGrid\" ";
  stream << "version=\"1.0\" ";
  stream << "compressor=\"vtkZLibDataCompressor\" ";
  stream << "header_type=\"UInt64\">\n";
}

static Vec3<int> get_nviz_cells(Grid3 const& cell_grid, int const q)
{
  Subgrid3 const owned_cells = get_owned_cells(cell_grid);
  Vec3<int> const ncells = owned_cells.extents();
  Vec3<int> const nviz_cells = q * ncells;
  return nviz_cells;
}

void write_vtr_rectilinear_start(std::stringstream& stream, Mesh const& mesh)
{
  int const q = mesh.basis().q;
  Grid3 const cell_grid = mesh.cell_grid();
  Vec3<int> const nviz_cells = get_nviz_cells(cell_grid, q);
  stream << "<RectilinearGrid WholeExtent=\"";
  stream << "0 " << nviz_cells.x() << " ";
  stream << "0 " << nviz_cells.y() << " ";
  stream << "0 " << nviz_cells.z() << "\">\n";
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

static void write_vtr_time(std::stringstream& stream, real const time)
{
  write_fdata_start(stream, "Float64", "TIME", 1);
  stream << time << "\n";
  write_data_end(stream);
}

static void write_vtr_step(std::stringstream& stream, int const step) {
  write_fdata_start(stream, "Int32", "STEP", 1);
  stream << step << "\n";
  write_data_end(stream);
}

static void write_vtr_block_level(
    std::stringstream& stream,
    int const block,
    Mesh const& mesh)
{
  int const dim = mesh.dim();
  tree::ID const global_id = mesh.owned_leaves()[block];
  std::int8_t level = tree::get_level(dim, global_id);
  write_fdata_start(stream, "Int32", "block_level", 1);
  stream << int(level) << "\n";
  write_data_end(stream);
}

static void write_vtr_block_ijk(
    std::stringstream& stream,
    int const block,
    Mesh const& mesh)
{
  int const dim = mesh.dim();
  tree::ID const global_id = mesh.owned_leaves()[block];
  tree::Point const point = tree::get_point(dim, global_id);
  Vec3<int> const ijk = point.ijk;
  write_fdata_start(stream, "Int32", "block_ijk", DIMENSIONS);
  stream << ijk.x() << " " << ijk.y() << " " << ijk.z() << "\n";
  write_data_end(stream);
}

static void write_vtr_zcurve_id(
    std::stringstream& stream,
    int const block,
    Mesh const& mesh)
{
  tree::ZLeaves const& z_leaves = mesh.z_leaves();
  tree::ID const global_id = mesh.owned_leaves()[block];
  int zcurve_id = -1;
  for (std::size_t i = 0; i < z_leaves.size(); ++i) {
    if (global_id == z_leaves[i]) {
      zcurve_id = int(i);
      break;
    }
  }
  write_fdata_start(stream, "Int32", "block_zcurve_id", 1);
  stream << zcurve_id << "\n";
  write_data_end(stream);
}

static void write_vtr_block_owner(
    std::stringstream& stream,
    int const owner)
{
  write_fdata_start(stream, "Int32", "block_owner", 1);
  stream << owner << "\n";
  write_data_end(stream);
}

static void write_vtr_field_data(
    std::stringstream& stream,
    int const block,
    Mesh const& mesh,
    real const time,
    int const step)
{
  int const rank = mesh.comm()->rank();
  stream << std::scientific << std::setprecision(12);
  stream << "<FieldData>\n";
  write_vtr_time(stream, time);
  write_vtr_step(stream, step);
  write_vtr_block_level(stream, block, mesh);
  write_vtr_block_ijk(stream, block, mesh);
  write_vtr_zcurve_id(stream, block, mesh);
  write_vtr_block_owner(stream, rank);
  stream << "</FieldData>\n";
}

static void write_vtr_piece_start(
    std::stringstream& stream,
    Mesh const& mesh)
{
  int const q = mesh.basis().q;
  Grid3 const cell_grid = mesh.cell_grid();
  Vec3<int> const nviz_cells = get_nviz_cells(cell_grid, q);
  stream << "<Piece Extent=\"";
  stream << "0 " << nviz_cells.x() << " ";
  stream << "0 " << nviz_cells.y() << " ";
  stream << "0 " << nviz_cells.z() << "\">\n";
}

template <class T>
void write_data(
    std::stringstream& stream,
    VtkView<T> dual,
    bool copy = true) {
  if (copy) dual.template sync<typename VtkView<T>::host_mirror_space>();
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

static void write_coordinate(
    std::stringstream& stream,
    int const block,
    Mesh const& mesh,
    int const axis)
{
  std::string const axis_name = get_axis_name(axis);
  int const q = mesh.basis().q;
  check_q(q);
  int const num_pts = get_nviz_cells(mesh.cell_grid(), q)[axis] + 1;
  Box3<real> const block_domain = mesh.block_info_h().domains[block];
  Vec3<real> const cell_dx = mesh.block_info_h().cell_dxs[block];
  real const o = block_domain.lower()[axis] + cell_dx[axis];
  real const dx = cell_dx[axis] / q;
  real const q4 = std::sqrt(30.)/36.;
  VtkView<float> coord;
  Kokkos::resize(coord, num_pts, 1);
  real const offset[max_1D_quadrature_points-1][max_1D_quadrature_points-1] = {
    {0.,  0.,    0.,    0.},
    {0.,  0.,    0.,    0.},
    {0., -2./9., 2./9., 0.},
    {0., -q4,    0.,    q4}
  };
  for (int i = 0; i < num_pts; ++i) {
    int const mod = i % q;
    coord.h_view(i, 0) = o + i*dx + offset[q-1][mod]*dx;
  }
  write_data_start(stream, "Float32", axis_name, 1);
  write_data(stream, coord, false);
  write_data_end(stream);
}

static void write_vtr_coordinates(
    std::stringstream& stream,
    int const block,
    Mesh const& mesh)
{
  stream << "<Coordinates>\n";
  write_coordinate(stream, block, mesh, X);
  write_coordinate(stream, block, mesh, Y);
  write_coordinate(stream, block, mesh, Z);
  stream << "</Coordinates>\n";
}

void write_vtr_start(
    std::stringstream& stream,
    int const block,
    Mesh const& mesh,
    real const time,
    int const step)
{
  write_vtr_header(stream);
  write_vtr_rectilinear_start(stream, mesh);
  write_vtr_field_data(stream, block, mesh, time, step);
  write_vtr_piece_start(stream, mesh);
  write_vtr_coordinates(stream, block, mesh);
  stream << "<CellData>\n";
}

void write_vtr_end(std::stringstream& stream)
{
  stream << "</CellData>\n";
  stream << "</Piece>\n";
  stream << "</RectilinearGrid>\n";
  stream << "</VTKFile>\n";
}

template <class T> std::string vtk_type_name();
template <> std::string vtk_type_name<int>() { return "Int32"; }
template <> std::string vtk_type_name<float>() { return "Float32"; }
template <> std::string vtk_type_name<double>() { return "Float64"; }

template <class T>
void write_vtr_field(
    std::stringstream& stream,
    std::string const& name,
    VtkView<T> f)
{
  int const ncomps = f.d_view.extent(1);
  write_data_start(stream, vtk_type_name<T>(), name, ncomps);
  write_data(stream, f, true);
  write_data_end(stream);
}

template void write_vtr_field<int>(std::stringstream&, std::string const&, VtkView<int>);
template void write_vtr_field<float>(std::stringstream&, std::string const&, VtkView<float>);
template void write_vtr_field<real>(std::stringstream&, std::string const&, VtkView<real>);

static void write_vtm_header(std::stringstream& stream) {
  stream << "<VTKFile type=\"vtkMultiBlockDataSet\" ";
  stream << "version=\"1.0\">\n";
  stream << "<vtkMultiBlockDataSet>\n";
}

static void write_vtm_source_file(
    std::stringstream& stream,
    int const i,
    std::string const& file) {
  stream << "<DataSet index=\"" << i << "\" ";
  stream << "file=\"" << file << "\"/>\n";
}

static void write_vtm_end(std::stringstream& stream) {
  stream << "</vtkMultiBlockDataSet>\n";
  stream << "</VTKFile>";
}

void write_vtm(
    std::stringstream& stream,
    std::string const& prefix,
    int const num_blocks)
{
  write_vtm_header(stream);
  for (int block = 0; block < num_blocks; ++block) {
    std::string const file = fmt::format("{}{}.vtr", prefix, block);
    write_vtm_source_file(stream, block, file);
  }
  write_vtm_end(stream);
}

Grid3 get_viz_cell_grid(Grid3 const& cell_grid, int const q)
{
  return Grid3(get_nviz_cells(cell_grid, q));
}

}
}
