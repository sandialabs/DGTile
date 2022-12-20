#include <fstream>
#include <iomanip>

#include "caliper/cali.h"

#include "dgt_file.hpp"
#include "dgt_grid.hpp"
#include "dgt_mesh.hpp"

namespace dgt {

namespace ascii {

static void verify_extension(std::filesystem::path const& path) {
  if (path.extension() != ".dga") {
    throw std::runtime_error("ascii::read_mesh - extension != .dga");
  }
}

static void verify_file(
    std::fstream const& f,
    std::filesystem::path const& path) {
  if (!f.is_open()) {
    throw std::runtime_error("ascii - could not open: " + path.string());
  }
}

template <class T>
void write_vec3(std::stringstream& stream, vector3<T> const& v) {
  stream << v.x() << " " << v.y() << " " << v.z();
}

static void write_point(std::stringstream& stream, Point const& pt) {
  stream << pt.depth << " ";
  write_vec3(stream, pt.ijk);
}

static void write_meta(std::filesystem::path const& base, Mesh const& mesh) {
  CALI_CXX_MARK_FUNCTION;
  std::stringstream stream;
  stream << std::scientific;
  stream << std::setprecision(17);
  stream << mesh.dim() << " ";
  stream << mesh.basis().p << "\n";
  stream << mesh.basis().tensor << "\n";
  stream << mesh.nsoln() << " ";
  stream << mesh.neq() << "\n";
  write_vec3(stream, mesh.domain().lower()); stream << "\n";
  write_vec3(stream, mesh.domain().upper()); stream << "\n";
  write_vec3(stream, mesh.periodic()); stream << "\n";
  write_vec3(stream, mesh.cell_grid().extents()); stream << "\n";
  write_point(stream, mesh.tree().base()); stream << "\n";
  stream << mesh.leaves().size() << "\n";
  for (Node* leaf : mesh.leaves()) {
    write_point(stream, leaf->pt()); stream << "\n";
  }
  std::filesystem::path const file_path = base / "mesh.dga";
  write_stream(file_path, stream);
}

static void read_meta(std::filesystem::path const& path, Mesh& mesh) {
  CALI_CXX_MARK_FUNCTION;
  std::filesystem::path const file_path = path / "mesh.dga";
  std::fstream file;
  file.open(file_path, std::ios::in);
  verify_file(file, file_path);
  bool tensor;
  int dim, p, nsoln, neq, nblocks;
  vector3<double> xmin, xmax;
  vector3<bool> periodic;
  vector3<int> cells;
  Point base;
  std::vector<Point> leaf_pts;
  file >> dim;
  file >> p;
  file >> tensor;
  file >> nsoln;
  file >> neq;
  file >> xmin.x() >> xmin.y() >> xmin.z();
  file >> xmax.x() >> xmax.y() >> xmax.z();
  file >> periodic.x() >> periodic.y() >> periodic.z();
  file >> cells.x() >> cells.y() >> cells.z();
  file >> base.depth >> base.ijk.x() >> base.ijk.y() >> base.ijk.z();
  file >> nblocks;
  leaf_pts.resize(nblocks);
  for (Point& pt : leaf_pts) {
    file >> pt.depth >> pt.ijk.x() >> pt.ijk.y() >> pt.ijk.z();
  }
  file.close();
  mesh.set_domain({xmin, xmax});
  mesh.set_periodic(periodic);
  mesh.set_cell_grid(cells);
  mesh.set_nsoln(nsoln);
  mesh.set_neq(neq);
  mesh.init(grid3(base.ijk), p, tensor);
  for (Point& pt : leaf_pts) {
    mesh.tree().insert(pt);
  }
  mesh.rebuild();
  mesh.allocate();
}

static void write_block(std::filesystem::path const& base, Block const& block) {
  CALI_CXX_MARK_FUNCTION;
  std::stringstream stream;
  stream << std::scientific;
  stream << std::setprecision(17);
  grid3 const cell_grid = generalize(block.cell_grid());
  View<double***> U = block.soln(0);
  auto U_host = Kokkos::create_mirror_view(U);
  Kokkos::deep_copy(U_host, U);
  Basis const b = block.basis();
  int const neq = U.extent(1);
  auto f = [&] (vector3<int> const& cell_ijk) {
    int const cell = cell_grid.index(cell_ijk);
    for (int eq = 0; eq < neq; ++eq) {
      for (int m = 0; m < b.nmodes; ++m) {
        stream << std::setw(24) << U_host(cell, eq, m) << " ";
      }
      stream << "\n";
    }
  };
  for_each(execution::seq, cell_grid, f);
  int const id = block.id();
  std::string const block_name = std::to_string(id) + ".dga";
  std::filesystem::path const file_path = base / block_name;
  write_stream(file_path, stream);
}

static void read_block(std::filesystem::path const& base, Block& block) {
  CALI_CXX_MARK_FUNCTION;
  int const id = block.id();
  std::string const block_name = std::to_string(id) + ".dga";
  std::filesystem::path const file_path = base / block_name;
  std::fstream file;
  file.open(file_path, std::ios::in);
  verify_file(file, file_path);
  grid3 const cell_grid = generalize(block.cell_grid());
  View<double***> U = block.soln(0);
  auto U_host = Kokkos::create_mirror_view(U);
  Basis const b = block.basis();
  int const neq = U.extent(1);
  auto f = [&] (vector3<int> const& cell_ijk) {
    int const cell = cell_grid.index(cell_ijk);
    for (int eq = 0; eq < neq; ++eq) {
      for (int m = 0; m < b.nmodes; ++m) {
        file >> U_host(cell, eq, m);
      }
    }
  };
  for_each(execution::seq, cell_grid, f);
  file.close();
  Kokkos::deep_copy(U, U_host);
}

void write_mesh(std::filesystem::path const& path, Mesh const& mesh) {
  CALI_CXX_MARK_FUNCTION;
  std::filesystem::path const base(path.string() + ".dga");
  std::filesystem::create_directory(base);
  if (mesh.comm()->rank() == 0) {
    write_meta(base, mesh);
  }
  for (Node* leaf : mesh.owned_leaves()) {
    write_block(base, leaf->block);
  }
}

void read_mesh(std::filesystem::path const& path, Mesh& mesh) {
  CALI_CXX_MARK_FUNCTION;
  verify_extension(path);
  read_meta(path, mesh);
  for (Node* leaf : mesh.owned_leaves()) {
    read_block(path, leaf->block);
  }
}

}

}
