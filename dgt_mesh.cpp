#include <stdexcept>

#include "p3a_for_each.hpp"

#include "caliper/cali.h"

#include "dgt_amr.hpp"
#include "dgt_grid.hpp"
#include "dgt_mesh.hpp"
#include "dgt_spatial.hpp"

namespace dgt {

static void verify_comm(mpicpp::comm* comm) {
  if (!comm) {
    throw std::runtime_error("Mesh - unset comm");
  }
}

static void verify_cell_grid(p3a::grid3 const& cell_grid) {
  p3a::grid3 const gg = generalize(cell_grid);
  if (gg.extents().volume() == 0) {
    throw std::runtime_error("Mesh - unset cell grid");
  }
}

static void verify_dims(p3a::grid3 const& cell_grid, p3a::grid3 const& block_grid) {
  int const cell_dim = get_dim(cell_grid);
  int const block_dim = get_dim(block_grid);
  if (cell_dim != block_dim) {
    throw std::runtime_error("Mesh - invalid grids");
  }
}

static void verify_domain(int dim, p3a::box3<double> const& domain) {
  for (int axis = 0; axis < dim; ++axis) {
    if (domain.lower()[axis] >= domain.upper()[axis]) {
      throw std::runtime_error("Mesh - xmin >= xmax");
    }
  }
  for (int axis = dim; axis < DIMS; ++axis) {
    if ((domain.lower()[axis] != 0) || domain.upper()[axis] != 0) {
      throw std::runtime_error("Mesh - domain dim != mesh dim");
    }
  }
}

static void verify_leaves(mpicpp::comm* comm, std::vector<Node*> leaves) {
  int const nranks = comm->size();
  int const nleaves = leaves.size();
  if (nranks > nleaves) {
    throw std::runtime_error("Mesh - too few leaves");
  }
}

static void verify_initial_tree(Tree& tree) {
  if (!tree.root()->is_leaf()) {
    throw std::runtime_error("Mesh - tree exists");
  }
}

static void verify_solution(int nsoln, int nmodal_eq, int nflux_eq) {
  if (nsoln < 1) {
    throw std::runtime_error("Mesh - invalid nsoln");
  }
  if (nmodal_eq < 1) {
    throw std::runtime_error("Mesh - invalid nmodal_eq");
  }
  if (nflux_eq < 0) {
    throw std::runtime_error("Mesh - invalid nflux_eq");
  }
}

static void verify_fine_to_coarse_node(Node* n) {
  if (!n) {
    throw std::runtime_error("Mesh - not 2to1 in FTC");
  }
}

static void verify_no_field(
    std::string const& name,
    std::vector<FieldInfo> const& fields) {
  for (FieldInfo const& f : fields) {
    if (f.name == name) {
      throw std::runtime_error("Mesh - field " + name + " exists");
    }
  }
}

static int get_max_ijk(Point const& pt, Point const& base, int axis) {
  int const nblocks = base.ijk[axis];
  int const diff = pt.depth - base.depth;
  int const max_ijk = int(std::pow(2, diff)*nblocks) - 1;
  return max_ijk;
}

static void verify_boundary_ijk(
    Point const& pt, Point const& base, int axis, int dir) {
  int const min_ijk = 0;
  int const max_ijk = get_max_ijk(pt, base, axis);
  bool is_left = (pt.ijk[axis] == min_ijk) && (dir == left);
  bool is_right = (pt.ijk[axis] == max_ijk) && (dir == right);
  if (!(is_left || is_right)) {
    throw std::runtime_error("Mesh - invalid tree");
  }
}

mpicpp::comm* Mesh::comm() const {
  return m_comm;
}

Tree& Mesh::tree() {
  return m_tree;
}

Tree const& Mesh::tree() const {
  return m_tree;
}

int Mesh::dim() const {
  return m_tree.dim();
}

p3a::box3<double> Mesh::domain() const {
  return m_domain;
}

p3a::vector3<bool> Mesh::periodic() const {
  return m_periodic;
}

p3a::grid3 Mesh::cell_grid() const {
  return m_cell_grid;
}

Basis& Mesh::basis() {
  return m_basis;
}

Basis const& Mesh::basis() const {
  return m_basis;
}

int Mesh::nmodal_eq() const {
  return m_nmodal_eq;
}

int Mesh::nflux_eq() const {
  return m_nflux_eq;
}

int Mesh::nsoln() const {
  return m_nsoln;
}

std::vector<Node*> const& Mesh::leaves() const {
  return m_leaves;
}

std::vector<Node*> const& Mesh::owned_leaves() const {
  return m_owned_leaves;
}

std::vector<FieldInfo> const& Mesh::fields() const {
  return m_fields;
}

void Mesh::set_comm(mpicpp::comm* comm) {
  m_comm = comm;
}

void Mesh::set_domain(p3a::box3<double> const& domain) {
  m_domain = domain;
}

void Mesh::set_periodic(p3a::vector3<bool> const& periodic) {
  m_periodic = periodic;
}

void Mesh::set_cell_grid(p3a::grid3 const& cell_grid) {
  m_cell_grid = cell_grid;
}

void Mesh::set_nsoln(int nsoln) {
  m_nsoln = nsoln;
}

void Mesh::set_nmodal_eq(int neq) {
  m_nmodal_eq = neq;
}

void Mesh::set_nflux_eq(int neq) {
  m_nflux_eq = neq;
}

void Mesh::set_tree(Tree& tree) {
  m_tree = std::move(tree);
}

void Mesh::add_field(std::string name, int ent_dim, int ncomps) {
  verify_no_field(name, m_fields);
  FieldInfo info;
  info.name = name;
  info.ent_dim = ent_dim;
  info.ncomps = ncomps;
  m_fields.push_back(info);
}

void Mesh::init(p3a::grid3 const& block_grid, int p, bool tensor) {
  CALI_CXX_MARK_FUNCTION;
  verify_initial_tree(tree());
  verify_dims(m_cell_grid, block_grid);
  m_tree.init(block_grid);
  m_basis.init(dim(), p, tensor);
}

static Point get_adj_pt(Point const& pt, int axis, int dir) {
  Point adj_pt = pt;
  std::array<int, 2> offset = {-1,1};
  adj_pt.ijk[axis] += offset[dir];
  return adj_pt;
}

static Point make_adj_periodic(
    Point const& pt, Point const& base, int axis, int dir) {
  Point periodic = pt;
  int const min_ijk = 0;
  int const max_ijk = get_max_ijk(pt, base, axis);
  if (dir == left) periodic.ijk[axis] = max_ijk;
  if (dir == right) periodic.ijk[axis] = min_ijk;
  return periodic;
}

static Point get_parent_pt(Point const& pt) {
  Point parent_pt;
  parent_pt.depth = pt.depth - 1;
  parent_pt.ijk = get_coarse_ijk(pt.ijk);
  return parent_pt;
}

static void init_borders(
    Tree& tree,
    Node* leaf,
    p3a::vector3<bool> const& periodic) {
  CALI_CXX_MARK_FUNCTION;
  Block& block = leaf->block;
  int const dim = tree.dim();
  Point const base = tree.base();
  Point const pt = leaf->pt();
  for (int axis = 0; axis < dim; ++axis) {
    for (int dir = 0; dir < ndirs; ++dir) {
      Border& border = block.border(axis, dir);
      border.reset();
      border.set_node(leaf);
      border.set_axis(axis);
      border.set_dir(dir);
      Point adj_pt = get_adj_pt(pt, axis, dir);
      bool const in = contains(base.depth, p3a::subgrid3(base.ijk), adj_pt);
      if (!in) {
        verify_boundary_ijk(pt, base, axis, dir);
        if (periodic[axis]) {
          adj_pt = make_adj_periodic(adj_pt, base, axis, dir);
        } else {
          border.set_type(BOUNDARY);
          continue;
        }
      }
      Node* adj = tree.find(adj_pt);
      if (adj) {
        if (adj->is_leaf()) {
          border.set_type(STANDARD);
          border.set_adj(adj);
        } else {
          border.set_type(COARSE_TO_FINE);
          border.set_adj(adj);
        }
      } else {
        Point const adj_parent_pt = get_parent_pt(adj_pt);
        Node* adj_parent = tree.find(adj_parent_pt);
        verify_fine_to_coarse_node(adj_parent);
        border.set_type(FINE_TO_COARSE);
        border.set_adj(adj_parent);
      }
    }
  }
}

static void init_fields(
    Node* leaf,
    std::vector<FieldInfo> const& fields) {
  Block& block = leaf->block;
  for (FieldInfo const& info : fields) {
    int const idx = block.field_idx(info.name);
    if (idx == -1) {
      leaf->block.add_field(info);
    }
  }
}

void init_leaves(
    Mesh* mesh,
    Tree& tree,
    p3a::vector3<bool> const& periodic,
    std::vector<Node*> const& leaves) {
  CALI_CXX_MARK_FUNCTION;
  for (Node* leaf : leaves) {
    leaf->block.set_mesh(mesh);
    leaf->block.set_node(leaf);
    init_borders(tree, leaf, periodic);
    init_fields(leaf, mesh->fields());
  }
}

void Mesh::rebuild() {
  CALI_CXX_MARK_FUNCTION;
  m_leaves.resize(0);
  m_owned_leaves.resize(0);
  verify_comm(m_comm);
  verify_cell_grid(m_cell_grid);
  verify_domain(get_dim(m_cell_grid), m_domain);
  m_leaves = collect_leaves(tree());
  partition_leaves(m_comm, m_leaves);
  init_leaves(this, m_tree, m_periodic, m_leaves);
  m_owned_leaves = collect_owned_leaves(m_comm, m_leaves);
}

void Mesh::scale(double l) {
  m_domain.lower() *= l;
  m_domain.upper() *= l;
}

void Mesh::verify() {
  verify_comm(m_comm);
  verify_cell_grid(m_cell_grid);
  verify_domain(get_dim(m_cell_grid), m_domain);
  verify_leaves(m_comm, m_leaves);
}

void Mesh::allocate() {
  CALI_CXX_MARK_FUNCTION;
  verify_solution(m_nsoln, m_nmodal_eq, m_nflux_eq);
  for (Node* leaf : m_owned_leaves) {
    leaf->block.allocate(m_nsoln, m_nmodal_eq, m_nflux_eq);
  }
}

static void free_branch_node(int dim, Node* node) {
  CALI_CXX_MARK_FUNCTION;
  if (!node->is_leaf()) {
    node->block.deallocate();
    node->block.reset();
  }
  auto f = [&] (p3a::vector3<int> const& local) {
    Node* child = node->child(local);
    if (child) free_branch_node(dim, child);
  };
  p3a::for_each(p3a::execution::seq, generalize(get_child_grid(dim)), f);
}

void Mesh::clean() {
  CALI_CXX_MARK_FUNCTION;
  free_branch_node(dim(), m_tree.root());
}

}
