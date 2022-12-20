#include "gtest/gtest.h"

#include "p3a_for_each.hpp"

#include "dgt_grid.hpp"
#include "dgt_tree.hpp"

static void recurse(int dim, dgt::Node* node) {
 if (node->is_leaf()) {
    int const depth = node->pt().depth;
    p3a::vector3<int> const ijk = node->pt().ijk;
    ASSERT_EQ(depth, 1);
    for (int axis = 0; axis < dim; ++axis) {
      ASSERT_TRUE(ijk[axis] == 0 || ijk[axis] == 1);
    }
    for (int axis = dim; axis < dgt::DIMS; ++axis) {
      ASSERT_TRUE(ijk[axis] == 0);
    }
  }
  auto f = [&] (p3a::vector3<int> const& c) {
    dgt::Node* child = node->child(c);
    if (child) recurse(dim, child);
  };
  p3a::for_each(p3a::execution::seq, dgt::generalize(dgt::get_child_grid(dim)), f);
}

TEST(tree, init_1D) {
  dgt::Tree tree;
  tree.init(p3a::grid3(2,0,0));
  ASSERT_EQ(tree.dim(), 1);
  recurse(tree.dim(), tree.root());
}

TEST(tree, init_2D) {
  dgt::Tree tree;
  tree.init(p3a::grid3(2,2,0));
  ASSERT_EQ(tree.dim(), 2);
  recurse(tree.dim(), tree.root());
}

TEST(tree, init_3D) {
  dgt::Tree tree;
  tree.init(p3a::grid3(2,2,2));
  ASSERT_EQ(tree.dim(), 3);
  recurse(tree.dim(), tree.root());
}

TEST(tree, find_1D) {
  dgt::Tree tree;
  tree.init(p3a::grid3(2,0,0));
  dgt::Node* in = tree.find({1, {1,0,0}});
  dgt::Node* out1 = tree.find({1, {2,0,0}});
  dgt::Node* out2 = tree.find({1, {1,0,1}});
  ASSERT_EQ(in->pt().depth, 1);
  ASSERT_EQ(in->pt().ijk, p3a::vector3<int>(1,0,0));
  ASSERT_EQ(out1, nullptr);
  ASSERT_EQ(out2, nullptr);
}

TEST(tree, find_2D) {
  dgt::Tree tree;
  tree.init(p3a::grid3(2,2,0));
  dgt::Node* in = tree.find({1, {1,1,0}});
  dgt::Node* out1 = tree.find({1, {2,0,0}});
  dgt::Node* out2 = tree.find({1, {1,1,1}});
  ASSERT_EQ(in->pt().depth, 1);
  ASSERT_EQ(in->pt().ijk, p3a::vector3<int>(1,1,0));
  ASSERT_EQ(out1, nullptr);
  ASSERT_EQ(out2, nullptr);
}

TEST(tree, find_3D) {
  dgt::Tree tree;
  tree.init(p3a::grid3(2,2,2));
  dgt::Node* in = tree.find({1, {1,1,1}});
  dgt::Node* out1 = tree.find({1, {2,0,0}});
  dgt::Node* out2 = tree.find({1, {1,1,2}});
  ASSERT_EQ(in->pt().depth, 1);
  ASSERT_EQ(in->pt().ijk, p3a::vector3<int>(1,1,1));
  ASSERT_EQ(out1, nullptr);
  ASSERT_EQ(out2, nullptr);
}

TEST(tree, move_2D) {
  dgt::Tree orig;
  orig.init(p3a::grid3(2,2,0));
  dgt::Tree moved = std::move(orig);
  dgt::Node* in = moved.find({1,{1,1,0}});
  dgt::Node* out1 = moved.find({1, {2,0,0}});
  dgt::Node* out2 = moved.find({1, {1,1,1}});
  ASSERT_EQ(in->pt().depth, 1);
  ASSERT_EQ(in->pt().ijk, p3a::vector3<int>(1,1,0));
  ASSERT_EQ(out1, nullptr);
  ASSERT_EQ(out2, nullptr);
}

TEST(tree, move_3D) {
  dgt::Tree orig;
  orig.init(p3a::grid3(2,2,2));
  dgt::Tree moved = std::move(orig);
  dgt::Node* in = moved.find({1,{1,1,1}});
  dgt::Node* out1 = moved.find({1, {2,0,0}});
  dgt::Node* out2 = moved.find({1, {1,1,2}});
  ASSERT_EQ(in->pt().depth, 1);
  ASSERT_EQ(in->pt().ijk, p3a::vector3<int>(1,1,1));
  ASSERT_EQ(out1, nullptr);
  ASSERT_EQ(out2, nullptr);
}
