#include "gtest/gtest.h"

#include "p3a_for_each.hpp"

#include "dgt_grid.hpp"
#include "dgt_tree.hpp"

using namespace dgt;

static void recurse(int dim, Node* node) {
 if (node->is_leaf()) {
    int const depth = node->pt().depth;
    vector3<int> const ijk = node->pt().ijk;
    ASSERT_EQ(depth, 1);
    for (int axis = 0; axis < dim; ++axis) {
      ASSERT_TRUE(ijk[axis] == 0 || ijk[axis] == 1);
    }
    for (int axis = dim; axis < DIMS; ++axis) {
      ASSERT_TRUE(ijk[axis] == 0);
    }
  }
  auto f = [&] (vector3<int> const& c) {
    Node* child = node->child(c);
    if (child) recurse(dim, child);
  };
  for_each(execution::seq, generalize(get_child_grid(dim)), f);
}

TEST(tree, init_1D) {
  Tree tree;
  tree.init(grid3(2,0,0));
  ASSERT_EQ(tree.dim(), 1);
  recurse(tree.dim(), tree.root());
}

TEST(tree, init_2D) {
  Tree tree;
  tree.init(grid3(2,2,0));
  ASSERT_EQ(tree.dim(), 2);
  recurse(tree.dim(), tree.root());
}

TEST(tree, init_3D) {
  Tree tree;
  tree.init(grid3(2,2,2));
  ASSERT_EQ(tree.dim(), 3);
  recurse(tree.dim(), tree.root());
}

TEST(tree, find_1D) {
  Tree tree;
  tree.init(grid3(2,0,0));
  Node* in = tree.find({1, {1,0,0}});
  Node* out1 = tree.find({1, {2,0,0}});
  Node* out2 = tree.find({1, {1,0,1}});
  ASSERT_EQ(in->pt().depth, 1);
  ASSERT_EQ(in->pt().ijk, vector3<int>(1,0,0));
  ASSERT_EQ(out1, nullptr);
  ASSERT_EQ(out2, nullptr);
}

TEST(tree, find_2D) {
  Tree tree;
  tree.init(grid3(2,2,0));
  Node* in = tree.find({1, {1,1,0}});
  Node* out1 = tree.find({1, {2,0,0}});
  Node* out2 = tree.find({1, {1,1,1}});
  ASSERT_EQ(in->pt().depth, 1);
  ASSERT_EQ(in->pt().ijk, vector3<int>(1,1,0));
  ASSERT_EQ(out1, nullptr);
  ASSERT_EQ(out2, nullptr);
}

TEST(tree, find_3D) {
  Tree tree;
  tree.init(grid3(2,2,2));
  Node* in = tree.find({1, {1,1,1}});
  Node* out1 = tree.find({1, {2,0,0}});
  Node* out2 = tree.find({1, {1,1,2}});
  ASSERT_EQ(in->pt().depth, 1);
  ASSERT_EQ(in->pt().ijk, vector3<int>(1,1,1));
  ASSERT_EQ(out1, nullptr);
  ASSERT_EQ(out2, nullptr);
}

TEST(tree, move_2D) {
  Tree orig;
  orig.init(grid3(2,2,0));
  Tree moved = std::move(orig);
  Node* in = moved.find({1,{1,1,0}});
  Node* out1 = moved.find({1, {2,0,0}});
  Node* out2 = moved.find({1, {1,1,1}});
  ASSERT_EQ(in->pt().depth, 1);
  ASSERT_EQ(in->pt().ijk, vector3<int>(1,1,0));
  ASSERT_EQ(out1, nullptr);
  ASSERT_EQ(out2, nullptr);
}

TEST(tree, move_3D) {
  Tree orig;
  orig.init(grid3(2,2,2));
  Tree moved = std::move(orig);
  Node* in = moved.find({1,{1,1,1}});
  Node* out1 = moved.find({1, {2,0,0}});
  Node* out2 = moved.find({1, {1,1,2}});
  ASSERT_EQ(in->pt().depth, 1);
  ASSERT_EQ(in->pt().ijk, vector3<int>(1,1,1));
  ASSERT_EQ(out1, nullptr);
  ASSERT_EQ(out2, nullptr);
}
