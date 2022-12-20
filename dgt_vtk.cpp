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
template <> std::string vtk_type_name<double>() { return "Float64"; }

static void write_vtu_header(std::stringstream& out) {
  out << "<VTKFile type=\"UnstructuredGrid\" ";
  out << "header_type=\"UInt64\" ";
  out << "compressor=\"vtkZLibDataCompressor\">\n";
}

static void write_data_start(
    std::stringstream& stream,
    std::string const& type,
    std::string const& name,
    int ncomps,
    bool is_pdata = false) {
  stream << "<";
  if (is_pdata) stream << "P";
  stream << "DataArray type=\"" << type << "\" ";
  stream << "Name=\"" << name << "\" ";
  stream << "NumberOfComponents=\"" << ncomps << "\" ";
  stream << "format=\"binary\"";
  if (is_pdata) stream << "/";
  stream << ">\n";
}

static void write_data_end(std::stringstream& stream) {
  stream << "</DataArray>\n";
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

static void write_types(std::stringstream& stream, int dim, int num) {
  static constexpr std::int8_t vtk_types[] = {1,3,9,12};
  std::int8_t const type = vtk_types[dim];
  VizView<std::int8_t> types;
  Kokkos::resize(types, num, 1);
  for (int i = 0; i < num; ++i) {
    types.h_view(i, 0) = type;
  }
  write_data_start(stream, "Int8", "types", 1);
  write_data(stream, types, false);
  write_data_end(stream);
}

static void write_offsets(std::stringstream& stream, int num, int nents) {
  int offset = nents;
  VizView<int> offsets;
  Kokkos::resize(offsets, num, 1);
  for (int i = 0; i < num; ++i) {
    offsets.h_view(i, 0) = offset;
    offset += nents;
  }
  write_data_start(stream, "Int32", "offsets", 1);
  write_data(stream, offsets, false);
  write_data_end(stream);
}

static void write_tree_connectivity(std::stringstream& stream, int npoints) {
  VizView<int> connectivity;
  Kokkos::resize(connectivity, npoints, 1);
  for (int i = 0; i < npoints; ++i) {
    connectivity.h_view(i, 0) = i;
  }
  write_data_start(stream, "Int32", "connectivity", 1);
  write_data(stream, connectivity, false);
  write_data_end(stream);
}

static constexpr vector3<int> vtk_corners[] = {
  {0,0,0},
  {1,0,0},
  {1,1,0},
  {0,1,0},
  {0,0,1},
  {1,0,1},
  {1,1,1},
  {0,1,1}
};

static void write_tree_coords(
    std::stringstream& stream,
    Mesh const& mesh,
    int npoints,
    int ncorners) {
  int idx = 0;
  Point const base = mesh.tree().base();
  box3<double> const domain = mesh.domain();
  std::vector<Node*> const& leaves = mesh.leaves();
  VizView<double> coords;
  Kokkos::resize(coords, npoints, DIMS);
  for (Node* leaf : leaves) {
    Point const pt = leaf->pt();
    box3<double> const box = get_block_domain(base, pt, domain);
    vector3<double> const o = box.lower();
    vector3<double> const dx = box.extents();
    for (int c = 0; c < ncorners; ++c) {
      vector3<double> const x = o + hadamard_product(dx, vtk_corners[c]);
      for (int axis = 0; axis < DIMS; ++axis) {
        coords.h_view(idx, axis) = x[axis];
      }
      idx++;
    }
  }
  write_data_start(stream, "Float64", "coordinates", DIMS);
  write_data(stream, coords, false);
  write_data_end(stream);
}

static void write_leaf_ids(
    std::stringstream& stream,
    std::vector<Node*> const& leaves) {
  int const nleaves = leaves.size();
  VizView<int> ids;
  Kokkos::resize(ids, nleaves, 1);
  for (int i = 0; i < nleaves; ++i) {
    ids.h_view(i, 0) = leaves[i]->block.id();
  }
  write_data_start(stream, "Int32", "ID", 1);
  write_data(stream, ids, false);
  write_data_end(stream);
}

static void write_leaf_depths(
    std::stringstream& stream,
    std::vector<Node*> const& leaves) {
  int const nleaves = leaves.size();
  VizView<int> depths;
  Kokkos::resize(depths, nleaves, 1);
  for (int i = 0; i < nleaves; ++i) {
    depths.h_view(i, 0) = leaves[i]->pt().depth;
  }
  write_data_start(stream, "Int32", "depth", 1);
  write_data(stream, depths, false);
  write_data_end(stream);
}

static void write_leaf_ijks(
    std::stringstream& stream,
    std::vector<Node*> const& leaves) {
  int const nleaves = leaves.size();
  VizView<int> ijks;
  Kokkos::resize(ijks, nleaves, DIMS);
  for (int i = 0; i < nleaves; ++i) {
    vector3<int> const ijk = leaves[i]->pt().ijk;
    for (int axis = 0; axis < DIMS; ++axis) {
      ijks.h_view(i, axis) = ijk[axis];
    }
  }
  write_data_start(stream, "Int32", "ijk", DIMS);
  write_data(stream, ijks, false);
  write_data_end(stream);
}

static void write_leaf_owners(
    std::stringstream& stream,
    std::vector<Node*> const& leaves) {
  int const nleaves = leaves.size();
  VizView<int> owners;
  Kokkos::resize(owners, nleaves, 1);
  for (int i = 0; i < nleaves; ++i) {
    owners.h_view(i, 0) = leaves[i]->block.owner();
  }
  write_data_start(stream, "Int32", "owner", 1);
  write_data(stream, owners, false);
  write_data_end(stream);
}

static void write_block_connectivity(
    std::stringstream& stream,
    Block const& block,
    int ncorners) {
  VizView<int> connectivity;
  int const p = block.basis().p;
  grid3 const cell_grid = block.cell_grid();
  grid3 const viz_cell_grid = generalize(get_viz_cell_grid(cell_grid, p));
  grid3 const viz_point_grid = generalize(get_viz_point_grid(cell_grid, p));
  Kokkos::resize(connectivity, viz_cell_grid.size(), ncorners);
  auto f = [&] (vector3<int> const& viz_cell_ijk) {
    int const viz_cell = viz_cell_grid.index(viz_cell_ijk);
    for (int c = 0; c < ncorners; ++c) {
      vector3<int> const offset = vtk_corners[c];
      vector3<int> const viz_point_ijk = viz_cell_ijk + offset;
      int const point = viz_point_grid.index(viz_point_ijk);
      connectivity.h_view(viz_cell, c) = point;
    }
  };
  for_each(execution::seq, viz_cell_grid, f);
  write_data_start(stream, "Int32", "connectivity", 1);
  write_data(stream, connectivity, false);
  write_data_end(stream);
}

static void write_block_coords(
    std::stringstream& stream,
    Block const& block) {
  VizView<double> coords;
  int const p = block.basis().p;
  vector3<double> const dx = block.dx() / (p+1);
  vector3<double> const origin = block.domain().lower();
  grid3 const cell_grid = block.cell_grid();
  grid3 const viz_point_grid = generalize(get_viz_point_grid(cell_grid, p));
  Kokkos::resize(coords, viz_point_grid.size(), DIMS);
  auto f = [&] (vector3<int> const& viz_point_ijk) {
    int const point = viz_point_grid.index(viz_point_ijk);
    coords.h_view(point, X) = origin.x() + viz_point_ijk.x() * dx.x();
    coords.h_view(point, Y) = origin.y() + viz_point_ijk.y() * dx.y();
    coords.h_view(point, Z) = origin.z() + viz_point_ijk.z() * dx.z();
  };
  for_each(execution::seq, viz_point_grid, f);
  write_data_start(stream, "Float64", "coordinates", DIMS);
  write_data(stream, coords, false);
  write_data_end(stream);
}

void write_tree(std::filesystem::path const& path, Mesh const& mesh) {
  CALI_CXX_MARK_FUNCTION;
  std::stringstream stream;
  write_vtu_header(stream);
  int const dim = mesh.dim();
  int const nleaves = mesh.leaves().size();
  int const ncorners = ipow(2, dim);
  int const npoints = nleaves * ncorners;
  stream << "<UnstructuredGrid>\n";
  stream << "<Piece NumberOfPoints=\"" << npoints << "\" ";
  stream << "NumberOfCells=\"" << nleaves << "\">\n";
  stream << "<Cells>\n";
  write_types(stream, dim, nleaves);
  write_offsets(stream, nleaves, ncorners);
  write_tree_connectivity(stream, npoints);
  stream << "</Cells>\n";
  stream << "<Points>\n";
  write_tree_coords(stream, mesh, npoints, ncorners);
  stream << "</Points>\n";
  stream << "<CellData>\n";
  write_leaf_ids(stream, mesh.leaves());
  write_leaf_depths(stream, mesh.leaves());
  write_leaf_ijks(stream, mesh.leaves());
  write_leaf_owners(stream, mesh.leaves());
  stream << "</CellData>\n";
  stream << "<PointData>\n";
  stream << "</PointData>\n";
  stream << "</Piece>\n";
  stream << "</UnstructuredGrid>\n";
  stream << "</VTKFile>\n";
  std::filesystem::path const file_path(path.string() + ".vtu");
  write_stream(file_path, stream);
}

void write_pvtu_start(std::stringstream& stream, int nblocks) {
  CALI_CXX_MARK_FUNCTION;
  stream << "<VTKFile type=\"PUnstructuredGrid\">\n";
  stream << "<PUnstructuredGrid>\n";
  for (int block = 0; block < nblocks; ++block) {
    std::filesystem::path const block_path(std::to_string(block) + ".vtu");
    stream << "<Piece Source=\"" << block_path.string() << "\"/>\n";
  }
  stream << "<PPoints>\n";
  write_data_start(stream, "Float64", "coordinates", DIMS, true);
  stream << "</PPoints>\n";
  stream << "<PCellData>\n";
}

void write_pvtu_end(std::stringstream& stream) {
  CALI_CXX_MARK_FUNCTION;
  stream << "</PCellData>\n";
  stream << "<PPointData>\n";
  stream << "</PPointData>\n";
  stream << "</PUnstructuredGrid>\n";
  stream << "</VTKFile>\n";
}

void write_vtu_start(std::stringstream& stream, Block const& block) {
  CALI_CXX_MARK_FUNCTION;
  int const dim = block.dim();
  int const p = block.basis().p;
  grid3 const cell_grid = block.cell_grid();
  int const ncorners = ipow(2, dim);
  int const ncells = generalize(get_viz_cell_grid(cell_grid, p)).size();
  int const npoints = generalize(get_viz_point_grid(cell_grid, p)).size();
  write_vtu_header(stream);
  stream << "<UnstructuredGrid>\n";
  stream << "<Piece NumberOfPoints=\"" << npoints << "\" ";
  stream << "NumberOfCells=\"" << ncells << "\">\n";
  stream << "<Cells>\n";
  write_types(stream, dim, ncells);
  write_offsets(stream, ncells, ncorners);
  write_block_connectivity(stream, block, ncorners);
  stream << "</Cells>\n";
  stream << "<Points>\n";
  write_block_coords(stream, block);
  stream << "</Points>\n";
  stream << "<CellData>\n";
}

void write_vtu_end(std::stringstream& stream) {
  CALI_CXX_MARK_FUNCTION;
  stream << "</CellData>\n";
  stream << "<PointData>\n";
  stream << "</PointData>\n";
  stream << "</Piece>\n";
  stream << "</UnstructuredGrid>\n";
  stream << "</VTKFile>\n";
}

template <class T>
void write_pvtu_cell_field(
    std::stringstream& stream,
    std::string const& name,
    int ncomps) {
  CALI_CXX_MARK_FUNCTION;
  write_data_start(stream, vtk_type_name<T>(), name, ncomps, true);
}

template <class T>
void write_vtu_cell_field(
    std::stringstream& stream,
    std::string const& name,
    VizView<T> f) {
  CALI_CXX_MARK_FUNCTION;
  int const ncomps = f.d_view.extent(1);
  write_data_start(stream, vtk_type_name<T>(), name, ncomps);
  write_data(stream, f, true);
  write_data_end(stream);
}

template void write_pvtu_cell_field<int>(std::stringstream&, std::string const&, int);
template void write_pvtu_cell_field<double>(std::stringstream&, std::string const&, int);

template void write_vtu_cell_field<int>(std::stringstream&, std::string const&, VizView<int> f);
template void write_vtu_cell_field<double>(std::stringstream&, std::string const&, VizView<double> f);

}

}
