#include "gtest/gtest.h"

#include "p3a_for_each.hpp"

#include "dgt_mesh.hpp"

TEST(mesh, init_1D) {
  dgt::Mesh mesh;
  mpicpp::comm world = mpicpp::comm::world();
  mesh.set_comm(&world);
  mesh.set_domain({p3a::vector3<double>(0,0,0), p3a::vector3<double>(1,0,0)});
  mesh.set_cell_grid({1,0,0});
  mesh.init({2,0,0}, 1, true);
  mesh.rebuild();
  for (dgt::Node* leaf : mesh.leaves()) {
    ASSERT_EQ(leaf->block.dx(), p3a::vector3<double>(.5,0,0));
  }
  dgt::Node* node = mesh.tree().find({1, {1,0,0}});
  ASSERT_TRUE(node->is_leaf());
  ASSERT_EQ(&mesh, node->block.mesh());
  ASSERT_EQ(node, node->block.node());
  ASSERT_EQ(node->block.border(dgt::X,dgt::left).type(), dgt::STANDARD);
  ASSERT_EQ(node->block.border(dgt::X,dgt::left).adj(), mesh.tree().find({1, {0,0,0}}));
  ASSERT_EQ(node->block.border(dgt::X,dgt::right).type(), dgt::BOUNDARY);
  ASSERT_EQ(node->block.border(dgt::X,dgt::right).adj(), nullptr);
}

TEST(mesh, init_2D) {
  dgt::Mesh mesh;
  mpicpp::comm world = mpicpp::comm::world();
  mesh.set_comm(&world);
  mesh.set_domain({p3a::vector3<double>(0,0,0), p3a::vector3<double>(1,1,0)});
  mesh.set_cell_grid({1,1,0});
  mesh.init({2,2,0}, 1, true);
  mesh.rebuild();
  for (dgt::Node* leaf : mesh.leaves()) {
    ASSERT_EQ(leaf->block.dx(), p3a::vector3<double>(.5,.5,0));
  }
  dgt::Node* node = mesh.tree().find({1, {1,0,0}});
  ASSERT_TRUE(node->is_leaf());
  ASSERT_EQ(&mesh, node->block.mesh());
  ASSERT_EQ(node, node->block.node());
  ASSERT_EQ(node->block.border(dgt::X,dgt::left).type(), dgt::STANDARD);
  ASSERT_EQ(node->block.border(dgt::X,dgt::left).adj(), mesh.tree().find({1, {0,0,0}}));
  ASSERT_EQ(node->block.border(dgt::X,dgt::right).type(), dgt::BOUNDARY);
  ASSERT_EQ(node->block.border(dgt::X,dgt::right).adj(), nullptr);
}

TEST(mesh, init_3D) {
  dgt::Mesh mesh;
  mpicpp::comm world = mpicpp::comm::world();
  mesh.set_comm(&world);
  mesh.set_domain({p3a::vector3<double>(0,0,0), p3a::vector3<double>(1,1,1)});
  mesh.set_cell_grid({1,1,1});
  mesh.init({2,2,2}, 1, true);
  mesh.rebuild();
  for (dgt::Node* leaf : mesh.leaves()) {
    ASSERT_EQ(leaf->block.dx(), p3a::vector3<double>(.5,.5,.5));
  }
  dgt::Node* node = mesh.tree().find({1, {1,0,0}});
  ASSERT_TRUE(node->is_leaf());
  ASSERT_EQ(&mesh, node->block.mesh());
  ASSERT_EQ(node, node->block.node());
  ASSERT_EQ(node->block.border(dgt::X,dgt::left).type(), dgt::STANDARD);
  ASSERT_EQ(node->block.border(dgt::X,dgt::left).adj(), mesh.tree().find({1, {0,0,0}}));
  ASSERT_EQ(node->block.border(dgt::X,dgt::right).type(), dgt::BOUNDARY);
  ASSERT_EQ(node->block.border(dgt::X,dgt::right).adj(), nullptr);
}

TEST(mesh, scale) {
  dgt::Mesh mesh;
  mpicpp::comm world = mpicpp::comm::world();
  mesh.set_comm(&world);
  mesh.set_domain({p3a::vector3<double>(0,0,0), p3a::vector3<double>(1,1,1)});
  mesh.set_cell_grid({1,1,1});
  mesh.init({2,2,2}, 1, true);
  mesh.rebuild();
  mesh.scale(0.01);
  ASSERT_EQ(mesh.domain().lower(), p3a::vector3<double>(0., 0., 0.));
  ASSERT_EQ(mesh.domain().upper(), p3a::vector3<double>(0.01, 0.01, 0.01));
}

TEST(mesh, allocate) {
  dgt::Mesh mesh;
  mpicpp::comm world = mpicpp::comm::world();
  mesh.set_comm(&world);
  mesh.set_domain({p3a::vector3<double>(0,0,0), p3a::vector3<double>(1,1,1)});
  mesh.set_cell_grid({1,1,1});
  mesh.set_nsoln(2);
  mesh.set_nmodal_eq(5);
  mesh.set_nflux_eq(5);
  mesh.init({2,2,2}, 1, true);
  mesh.rebuild();
  mesh.allocate();
}
