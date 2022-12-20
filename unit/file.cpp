#include "gtest/gtest.h"

#include "dgt_amr.hpp"
#include "dgt_defines.hpp"
#include "dgt_file.hpp"
#include "dgt_mesh.hpp"

TEST(file, base64) {
  double in_data[3];
  double out_data[3];
  in_data[0] = 100.;
  in_data[1] = 2.123;
  in_data[2] = 1.23;
  std::uint64_t const bytes = sizeof(double) * static_cast<uint64_t>(3);
  std::string const encoded = dgt::base64::encode(in_data, bytes);
  dgt::base64::decode(encoded, out_data, bytes);
  for (int i = 0; i < 3; ++i) {
    ASSERT_EQ(in_data[i], out_data[i]);
  }
}

static void init_test_mesh(int dim, dgt::Mesh& m) {
  p3a::vector3<int> block_grid(0,0,0);
  p3a::vector3<int> cell_grid(0,0,0);
  p3a::box3<double> domain(p3a::vector3<double>(0,0,0), p3a::vector3<double>(0,0,0));
  for (int axis = 0; axis < dim; ++axis) {
    block_grid[axis] = 1;
    cell_grid[axis] = 2;
    domain.upper()[axis] = 1.;
  }
  m.set_domain(domain);
  m.set_cell_grid(cell_grid);
  m.set_nsoln(1);
  m.set_nmodal_eq(5);
  m.set_nflux_eq(5);
  m.init(block_grid, 1, true);
  m.add_field("test_field", dim-1, 3);
  m.rebuild();
  refine(m.dim(), m.leaves().back());
  m.rebuild();
  refine(m.dim(), m.leaves().back());
  m.rebuild();
  refine(m.dim(), m.leaves().back());
  m.rebuild();
}

TEST(file, vtk_write_tree_1D) {
  dgt::Mesh mesh;
  mpicpp::comm world = mpicpp::comm::world();
  mesh.set_comm(&world);
  init_test_mesh(1, mesh);
  dgt::vtk::write_tree("tree1D", mesh);
}

TEST(file, vtk_write_tree_2D) {
  dgt::Mesh mesh;
  mpicpp::comm world = mpicpp::comm::world();
  mesh.set_comm(&world);
  init_test_mesh(2, mesh);
  dgt::vtk::write_tree("tree2D", mesh);
}

TEST(file, vtk_write_tree_3D) {
  dgt::Mesh mesh;
  mpicpp::comm world = mpicpp::comm::world();
  mesh.set_comm(&world);
  init_test_mesh(3, mesh);
  dgt::vtk::write_tree("tree3D", mesh);
}

TEST(file, vtk_write_ascii_2D) {
  dgt::Mesh mesh;
  mpicpp::comm world = mpicpp::comm::world();
  mesh.set_comm(&world);
  init_test_mesh(2, mesh);
  mesh.allocate();
  dgt::ascii::write_mesh("test", mesh);
}

TEST(file, vtk_read_ascii_2D) {
  dgt::Mesh mesh;
  mpicpp::comm world = mpicpp::comm::world();
  mesh.set_comm(&world);
  dgt::ascii::read_mesh("test.dga", mesh);
  ASSERT_EQ(mesh.dim(), 2);
  ASSERT_EQ(mesh.basis().p, 1);
  ASSERT_EQ(mesh.owned_leaves().size(), 10);
  ASSERT_EQ(mesh.fields().size(), 1);
  ASSERT_EQ(mesh.fields()[0].name, "test_field");
  ASSERT_EQ(mesh.fields()[0].ent_dim, 1);
  ASSERT_EQ(mesh.fields()[0].ncomps, 3);
}
