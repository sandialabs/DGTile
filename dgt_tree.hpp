#pragma once

#include "mpicpp.hpp"

#include "p3a_grid3.hpp"

#include "dgt_block.hpp"
#include "dgt_defines.hpp"
#include "dgt_point.hpp"

namespace dgt {

using namespace p3a;

class Node {
  private:
    friend class Tree;
  private:
    Point m_pt = {0, {0,0,0}};
    Node* m_parent = nullptr;
    std::unique_ptr<Node> m_child[2][2][2] = {{{nullptr}}};
  public:
    Block block;
  public:
    Node() = default;
    Node(Node const& other) = delete;
    Node operator=(Node const& other) = delete;
    Node(Node&& other) = default;
    Node& operator=(Node&& other) = default;
    [[nodiscard]] Point pt() const;
    [[nodiscard]] Node* parent();
    [[nodiscard]] Node const* parent() const;
    [[nodiscard]] Node* child(vector3<int> const& local);
    [[nodiscard]] Node const* child(vector3<int> const& local) const;
    [[nodiscard]] bool is_leaf() const;
    void add_child(vector3<int> const& local);
    void rm_child(vector3<int> const& local);
  private:
    Node(Node* parent, vector3<int> const& local);
    void create(int dim, Point const& base);
    void insert(int dim, Point const& pt);
};

class Tree {
  private:
    int m_dim = 0;
    Point m_base_pt = {0, {0,0,0}};
    std::unique_ptr<Node> m_root;
  public:
    Tree();
    Tree(Tree const& other) = delete;
    Tree operator=(Tree const& other) = delete;
    Tree(Tree&& other) = default;
    Tree& operator=(Tree&& other) = default;
    [[nodiscard]] int dim() const;
    [[nodiscard]] Point base() const;
    [[nodiscard]] Node* root();
    [[nodiscard]] Node const* root() const;
    [[nodiscard]] Node* find(Point const& pt);
    void set_dim(int dim);
    void set_base(Point const& pt);
    void insert(Point const& pt);
    void init(grid3 const& base);
};

std::vector<Node*> collect_leaves(Tree& tree);

void partition_leaves(
    mpicpp::comm* comm,
    std::vector<Node*> const& leaves);

std::vector<Node*> collect_owned_leaves(
    mpicpp::comm* comm,
    std::vector<Node*> const& leaves);

}
