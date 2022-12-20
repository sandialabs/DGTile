#include "gtest/gtest.h"

#include "dgt_amr.hpp"
#include "dgt_defines.hpp"
#include "dgt_file.hpp"
#include "dgt_mesh.hpp"

using namespace dgt;

TEST(file, base64) {
  double in_data[3];
  double out_data[3];
  in_data[0] = 100.;
  in_data[1] = 2.123;
  in_data[2] = 1.23;
  std::uint64_t const bytes = sizeof(double) * static_cast<uint64_t>(3);
  std::string const encoded = base64::encode(in_data, bytes);
  base64::decode(encoded, out_data, bytes);
  for (int i = 0; i < 3; ++i) {
    ASSERT_EQ(in_data[i], out_data[i]);
  }
}

static void init_test_mesh(int dim, Mesh& m) {
  vector3<int> block_grid(0,0,0);
  vector3<int> cell_grid(0,0,0);
  box3<double> domain(vector3<double>(0,0,0), vector3<double>(0,0,0));
  for (int axis = 0; axis < dim; ++axis) {
    block_grid[axis] = 1;
    cell_grid[axis] = 2;
    domain.upper()[axis] = 1.;
  }
  m.set_domain(domain);
  m.set_cell_grid(cell_grid);
  m.set_nsoln(1);
  m.set_neq(5);
  m.init(block_grid, 1, true);
  m.rebuild();
  refine(m.dim(), m.leaves().back());
  m.rebuild();
  refine(m.dim(), m.leaves().back());
  m.rebuild();
  refine(m.dim(), m.leaves().back());
  m.rebuild();
}

TEST(file, vtk_write_tree_1D) {
  Mesh mesh;
  mpicpp::comm world = mpicpp::comm::world();
  mesh.set_comm(&world);
  init_test_mesh(1, mesh);
  vtk::write_tree("tree1D", mesh);
}

TEST(file, vtk_write_tree_2D) {
  Mesh mesh;
  mpicpp::comm world = mpicpp::comm::world();
  mesh.set_comm(&world);
  init_test_mesh(2, mesh);
  vtk::write_tree("tree2D", mesh);
}

TEST(file, vtk_write_tree_3D) {
  Mesh mesh;
  mpicpp::comm world = mpicpp::comm::world();
  mesh.set_comm(&world);
  init_test_mesh(3, mesh);
  vtk::write_tree("tree3D", mesh);
}

TEST(file, vtk_write_ascii_2D) {
  Mesh mesh;
  mpicpp::comm world = mpicpp::comm::world();
  mesh.set_comm(&world);
  init_test_mesh(2, mesh);
  mesh.allocate();
  ascii::write_mesh("test", mesh);
}

TEST(file, vtk_read_ascii_2D) {
  Mesh mesh;
  mpicpp::comm world = mpicpp::comm::world();
  mesh.set_comm(&world);
  ascii::read_mesh("test.dga", mesh);
  ASSERT_EQ(mesh.dim(), 2);
  ASSERT_EQ(mesh.basis().p, 1);
  ASSERT_EQ(mesh.owned_leaves().size(), 10);
}
