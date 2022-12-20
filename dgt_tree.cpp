#include "mpicpp.hpp"

#include "p3a_for_each.hpp"

#include "caliper/cali.h"

#include "dgt_amr.hpp"
#include "dgt_grid.hpp"
#include "dgt_tree.hpp"

namespace dgt {

Point Node::pt() const {
  return m_pt;
}

Node* Node::parent() {
  return m_parent;
}

Node const* Node::parent() const {
  return m_parent;
}

Node* Node::child(p3a::vector3<int> const& local) {
  return m_child[local.x()][local.y()][local.z()].get();
}

Node const* Node::child(p3a::vector3<int> const& local) const {
  return m_child[local.x()][local.y()][local.z()].get();
}

bool Node::is_leaf() const {
  return !(m_child[0][0][0]);
}

void Node::add_child(p3a::vector3<int> const& local) {
  m_child[local.x()][local.y()][local.z()] =
    std::unique_ptr<Node>(new Node(this, local));
}

void Node::rm_child(p3a::vector3<int> const& local) {
  m_child[local.x()][local.y()][local.z()].reset();
}

Node::Node(Node* parent, p3a::vector3<int> const& local) {
  m_parent = parent;
  m_pt = get_child_point(m_parent->pt(), local);
}

void Node::create(int dim, Point const& base) {
  if (base.depth == m_pt.depth) return;
  auto f = [&] (p3a::vector3<int> const& local) {
    Point const child_pt = get_child_point(m_pt, local);
    if (contains(base.depth, p3a::subgrid3(base.ijk), child_pt)) {
      m_child[local.x()][local.y()][local.z()] =
        std::unique_ptr<Node>(new Node(this, local));
      m_child[local.x()][local.y()][local.z()]->create(dim, base);
    }
  };
  p3a::for_each(p3a::execution::seq, generalize(get_child_grid(dim)), f);
}

void Node::insert(int dim, Point const& pt) {
  if (pt.depth == m_pt.depth) return;
  int const child_depth = m_pt.depth + 1;
  auto f = [&] (p3a::vector3<int> const& local) {
    p3a::vector3<int> const child_ijk = 2*m_pt.ijk + local;
    p3a::subgrid3 const s(child_ijk, child_ijk + p3a::vector3<int>::ones());
    if (contains(child_depth, s, pt)) {
      if (!m_child[local.x()][local.y()][local.z()]) {
        m_child[local.x()][local.y()][local.z()] =
          std::unique_ptr<Node>(new Node(this, local));
      }
      m_child[local.x()][local.y()][local.z()]->insert(dim, pt);
    }
  };
  p3a::for_each(p3a::execution::seq, generalize(get_child_grid(dim)), f);
}

static int get_tree_depth(p3a::grid3 const& grid) {
  int depth = 0;
  p3a::vector3<int> const ex = grid.extents();
  int const limit = p3a::max(ex.x(), p3a::max(ex.y(), ex.z()));
  while ((1 << depth) < limit) {
    depth++;
  }
  return depth;
}

int Tree::dim() const {
  return m_dim;
}

Point Tree::base() const {
  return m_base_pt;
}

Node* Tree::root() {
  return m_root.get();
}

Node const* Tree::root() const {
  return m_root.get();
}

Tree::Tree() {
  m_root = std::make_unique<Node>();
}

static Node* find_node(Node* node, Point const& pt) {
  int const node_depth = node->pt().depth;
  p3a::vector3<int> const node_begin = node->pt().ijk;
  p3a::vector3<int> const node_end = node_begin + p3a::vector3<int>::ones();
  p3a::subgrid3 const subgrid(node_begin, node_end);
  if (!contains(node_depth, subgrid, pt)) return nullptr;
  if (node_depth == pt.depth) return node;
  int const shift = pt.depth - node_depth - 1;
  p3a::vector3<int> const local = {
    static_cast<int>((pt.ijk.x() >> shift) & 1L),
    static_cast<int>((pt.ijk.y() >> shift) & 1L),
    static_cast<int>((pt.ijk.z() >> shift) & 1L) };
  Node* child = node->child(local);
  if (!child) return nullptr;
  return find_node(child, pt);
}

Node* Tree::find(Point const& pt) {
  return find_node(m_root.get(), pt);
}

void Tree::set_dim(int dim) {
  m_dim = dim;
}

void Tree::set_base(Point const& pt) {
  m_base_pt = pt;
}

void Tree::insert(Point const& pt) {
  m_root->insert(m_dim, pt);
}

void Tree::init(p3a::grid3 const& base) {
  CALI_CXX_MARK_FUNCTION;
  m_dim = get_dim(base);
  m_base_pt.depth = get_tree_depth(base);
  m_base_pt.ijk = base.extents();
  m_root->create(m_dim, m_base_pt);
}

static void collect_leaves(int dim, Node* node, std::vector<Node*>& leaves) {
  if (node->is_leaf()) {
    node->block.reset();
    leaves.push_back(node);
  }
  auto f = [&] (p3a::vector3<int> const& local) {
    Node* child = node->child(local);
    if (child) collect_leaves(dim, child, leaves);
  };
  p3a::for_each(p3a::execution::seq, generalize(get_child_grid(dim)), f);
}

std::vector<Node*> collect_leaves(Tree& tree) {
  CALI_CXX_MARK_FUNCTION;
  std::vector<Node*> leaves;
  int const dim = tree.dim();
  collect_leaves(dim, tree.root(), leaves);
  return leaves;
}

void partition_leaves(mpicpp::comm* comm, std::vector<Node*> const& leaves) {
  CALI_CXX_MARK_FUNCTION;
  int const nleaves = leaves.size();
  int const nranks = comm->size();
  for (int rank = 0; rank < nranks; ++rank) {
    int const nowned = get_num_local(nleaves, nranks, rank);
    int const offset = get_local_offset(nleaves, nranks, rank);
    for (int i = 0; i < nowned; ++i) {
      int const id = offset + i;
      leaves[id]->block.set_owner(rank);
      leaves[id]->block.set_id(id);
    }
  }
}

std::vector<Node*> collect_owned_leaves(
    mpicpp::comm* comm,
    std::vector<Node*> const& leaves) {
  CALI_CXX_MARK_FUNCTION;
  std::vector<Node*> owned_leaves;
  int const rank = comm->rank();
  for (Node* leaf : leaves) {
    if (leaf->block.owner() == rank) {
      owned_leaves.push_back(leaf);
    }
  }
  return owned_leaves;
}

}
