#include "gtest/gtest.h"

#include "p3a_for_each.hpp"

#include "dgt_mesh.hpp"

using namespace dgt;

TEST(mesh, init_1D) {
  Mesh mesh;
  mpicpp::comm world = mpicpp::comm::world();
  mesh.set_comm(&world);
  mesh.set_domain({vector3<double>(0,0,0), vector3<double>(1,0,0)});
  mesh.set_cell_grid({1,0,0});
  mesh.init({2,0,0}, 1, true);
  mesh.rebuild();

  for (Node* leaf : mesh.leaves()) {
    ASSERT_EQ(leaf->block.dx(), vector3<double>(.5,0,0));
  }
  Node* node = mesh.tree().find({1, {1,0,0}});
  ASSERT_TRUE(node->is_leaf());
  ASSERT_EQ(&mesh, node->block.mesh());
  ASSERT_EQ(node, node->block.node());
  ASSERT_EQ(node->block.border(X,left).type(), STANDARD);
  ASSERT_EQ(node->block.border(X,left).adj(), mesh.tree().find({1, {0,0,0}}));
  ASSERT_EQ(node->block.border(X,right).type(), BOUNDARY);
  ASSERT_EQ(node->block.border(X,right).adj(), nullptr);
}

TEST(mesh, init_2D) {
  Mesh mesh;
  mpicpp::comm world = mpicpp::comm::world();
  mesh.set_comm(&world);
  mesh.set_domain({vector3<double>(0,0,0), vector3<double>(1,1,0)});
  mesh.set_cell_grid({1,1,0});
  mesh.init({2,2,0}, 1, true);
  mesh.rebuild();
  for (Node* leaf : mesh.leaves()) {
    ASSERT_EQ(leaf->block.dx(), vector3<double>(.5,.5,0));
  }
  Node* node = mesh.tree().find({1, {1,0,0}});
  ASSERT_TRUE(node->is_leaf());
  ASSERT_EQ(&mesh, node->block.mesh());
  ASSERT_EQ(node, node->block.node());
  ASSERT_EQ(node->block.border(X,left).type(), STANDARD);
  ASSERT_EQ(node->block.border(X,left).adj(), mesh.tree().find({1, {0,0,0}}));
  ASSERT_EQ(node->block.border(X,right).type(), BOUNDARY);
  ASSERT_EQ(node->block.border(X,right).adj(), nullptr);
}

TEST(mesh, init_3D) {
  Mesh mesh;
  mpicpp::comm world = mpicpp::comm::world();
  mesh.set_comm(&world);
  mesh.set_domain({vector3<double>(0,0,0), vector3<double>(1,1,1)});
  mesh.set_cell_grid({1,1,1});
  mesh.init({2,2,2}, 1, true);
  mesh.rebuild();
  for (Node* leaf : mesh.leaves()) {
    ASSERT_EQ(leaf->block.dx(), vector3<double>(.5,.5,.5));
  }
  Node* node = mesh.tree().find({1, {1,0,0}});
  ASSERT_TRUE(node->is_leaf());
  ASSERT_EQ(&mesh, node->block.mesh());
  ASSERT_EQ(node, node->block.node());
  ASSERT_EQ(node->block.border(X,left).type(), STANDARD);
  ASSERT_EQ(node->block.border(X,left).adj(), mesh.tree().find({1, {0,0,0}}));
  ASSERT_EQ(node->block.border(X,right).type(), BOUNDARY);
  ASSERT_EQ(node->block.border(X,right).adj(), nullptr);
}

TEST(mesh, scale) {
  Mesh mesh;
  mpicpp::comm world = mpicpp::comm::world();
  mesh.set_comm(&world);
  mesh.set_domain({vector3<double>(0,0,0), vector3<double>(1,1,1)});
  mesh.set_cell_grid({1,1,1});
  mesh.init({2,2,2}, 1, true);
  mesh.rebuild();
  mesh.scale(0.01);
  ASSERT_EQ(mesh.domain().lower(), vector3<double>(0., 0., 0.));
  ASSERT_EQ(mesh.domain().upper(), vector3<double>(0.01, 0.01, 0.01));
}

TEST(mesh, allocate) {
  Mesh mesh;
  mpicpp::comm world = mpicpp::comm::world();
  mesh.set_comm(&world);
  mesh.set_domain({vector3<double>(0,0,0), vector3<double>(1,1,1)});
  mesh.set_cell_grid({1,1,1});
  mesh.set_nsoln(2);
  mesh.set_neq(5);
  mesh.init({2,2,2}, 1, true);
  mesh.rebuild();
  mesh.allocate();
}
