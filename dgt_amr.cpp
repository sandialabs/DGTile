#include <algorithm>
#include <list>
#include <stdexcept>

#include "caliper/cali.h"

#include "p3a_for_each.hpp"

#include "dgt_amr.hpp"
#include "dgt_grid.hpp"
#include "dgt_interp.hpp"
#include "dgt_marks.hpp"
#include "dgt_mesh.hpp"
#include "dgt_spatial.hpp"

namespace dgt {

static void verify_refine(Node* parent) {
  if (!parent->is_leaf()) {
    throw std::runtime_error("refine - invalid parent");
  }
  auto f = [&] (p3a::vector3<int> const& local) {
    Node* child = parent->child(local);
    if (child) {
      throw std::runtime_error("refine - invalid parent");
    }
  };
  p3a::for_each(p3a::execution::seq, get_child_grid(DIMS), f);
}

static void verify_coarsen(Node* parent) {
  if (parent->is_leaf()) {
    throw std::runtime_error("coarsen - invalid parent, not leaf");
  }
}

static void verify_marks(Mesh const& mesh, std::vector<int8_t> const& marks) {
  if (marks.size() != mesh.leaves().size()) {
    throw std::runtime_error("modify - invalid marks");
  }
}

static void verify_mesh(Mesh const& mesh) {
  int const dim = mesh.dim();
  p3a::vector3<int> const ncells = mesh.cell_grid().extents();
  for (int axis = 0; axis < dim; ++axis) {
    if ((ncells[axis] % 2) != 0) {
      throw std::runtime_error("modify - invalid cell grid");
    }
  }
}

static void verify_leaf(Node const* leaf) {
  if (!leaf) {
    throw std::runtime_error("modify - leaf doesn't exist");
  }
  if (!leaf->is_leaf()) {
    throw std::runtime_error("modify - not a leaf");
  }
}

static void verify_2to1_parent(Node const* leaf) {
  if (!leaf) {
    throw std::runtime_error("ensure2to1 - leaf doesn't exist");
  }
  if (!leaf->is_leaf()) {
    throw std::runtime_error("ensure2to1 - something crazy happened");
  }
}

static void verify_marks(
    std::vector<int8_t> const& marks,
    std::vector<Node*> const& leaves) {
  if (marks.size() != leaves.size()) {
    throw std::runtime_error("modify - invalid marks");
  }
}

static void verify_insertion(
    View<double***> from,
    View<double***> to,
    p3a::grid3 const& from_grid,
    p3a::grid3 const& to_grid,
    p3a::subgrid3 const& from_subgrid,
    p3a::subgrid3 const& to_subgrid) {
  bool valid = true;
  p3a::grid3 const fg = generalize(from_grid);
  p3a::grid3 const tg = generalize(to_grid);
  p3a::subgrid3 const fsg = generalize(from_subgrid);
  p3a::subgrid3 const tsg = generalize(to_subgrid);
  if (int(from.extent(0)) != fg.size()) valid = false;
  if (int(to.extent(0)) != tg.size()) valid = false;
  if (from.extent(1) != to.extent(1)) valid = false;
  if (from.extent(2) != to.extent(2)) valid = false;
  if (fsg.size() != tsg.size()) valid = false;
  if (!valid) {
    throw std::runtime_error("amr insertion - invalid inputs");
  }
}

static void verify_transfer(
    Basis const& b,
    View<double***> from,
    View<double***> to,
    p3a::grid3 const& from_grid,
    p3a::grid3 const& to_grid,
    p3a::subgrid3 const& coarse_subgrid,
    p3a::subgrid3 const& fine_subgrid) {
  bool valid = true;
  p3a::subgrid3 const csg = generalize(coarse_subgrid);
  p3a::subgrid3 const fsg = generalize(fine_subgrid);
  if (int(from.extent(0)) != generalize(from_grid).size()) valid = false;
  if (int(to.extent(0)) != generalize(to_grid).size()) valid = false;
  if (int(from.extent(1)) != int(to.extent(1))) valid = false;
  if (int(from.extent(2)) != b.nmodes) valid = false;
  if (int(to.extent(2)) != b.nmodes) valid = false;
  if (fsg.size() != num_child(b.dim)*csg.size()) valid = false;
  if (!valid) {
    throw std::runtime_error("amr transfer - invalid inputs");
  }
}

void refine(int dim, Node* parent) {
  verify_refine(parent);
  auto f = [&] (p3a::vector3<int> const& local) {
    parent->add_child(local);
  };
  p3a::for_each(p3a::execution::seq, generalize(get_child_grid(dim)), f);
}

void coarsen(int dim, Node* parent) {
  verify_coarsen(parent);
  auto f = [&] (p3a::vector3<int> const& local) {
    parent->rm_child(local);
  };
  p3a::for_each(p3a::execution::seq, generalize(get_child_grid(dim)), f);
}

void do_insertion(
    View<double***> from,
    View<double***> to,
    p3a::grid3 const& from_grid,
    p3a::grid3 const& to_grid,
    p3a::subgrid3 const& from_subgrid,
    p3a::subgrid3 const& to_subgrid) {
  CALI_CXX_MARK_FUNCTION;
  verify_insertion(from, to, from_grid, to_grid, from_subgrid, to_subgrid);
  int const extent1 = from.extent(1);
  int const extent2 = from.extent(2);
  p3a::grid3 const general_from_grid = generalize(from_grid);
  p3a::grid3 const general_to_grid = generalize(to_grid);
  auto f = [=] P3A_DEVICE (p3a::vector3<int> const& from_cell_ijk) {
    p3a::vector3<int> const offset = from_cell_ijk - from_subgrid.lower();
    p3a::vector3<int> const to_cell_ijk = offset + to_subgrid.lower();
    int const from_cell = general_from_grid.index(from_cell_ijk);
    int const to_cell = general_to_grid.index(to_cell_ijk);
    for (int i = 0; i < extent1; ++i) {
      for (int j = 0; j < extent2; ++j) {
        to(to_cell, i, j) = from(from_cell, i, j);
      }
    }
  };
  p3a::for_each(p3a::execution::par, generalize(from_subgrid), f);
}

void do_prolongation(
    Basis const& b,
    View<double***> from,
    View<double***> to,
    p3a::grid3 const& from_grid,
    p3a::grid3 const& to_grid,
    p3a::subgrid3 const& from_subgrid,
    p3a::subgrid3 const& to_subgrid) {
  CALI_CXX_MARK_FUNCTION;
  verify_transfer(b, from, to, from_grid, to_grid, from_subgrid, to_subgrid);
  int const neq = from.extent(1);
  int const nchild = num_child(b.dim);
  int const nintr_pts = num_pts(b.dim, b.p);
  p3a::grid3 const general_from_grid = generalize(from_grid);
  p3a::grid3 const general_to_grid = generalize(to_grid);
  auto f = [=] P3A_DEVICE (p3a::vector3<int> const& from_cell_ijk) {
    int const from_cell = general_from_grid.index(from_cell_ijk);
    p3a::vector3<int> const coarse_offset = from_cell_ijk - from_subgrid.lower();
    for (int child = 0; child < nchild; ++child) {
      p3a::vector3<int> const local = get_local(child);
      p3a::vector3<int> const fine_offset = get_fine_ijk(coarse_offset, local);
      p3a::vector3<int> const to_cell_ijk = to_subgrid.lower() + fine_offset;
      int const to_cell = general_to_grid.index(to_cell_ijk);
      for (int pt = 0; pt < nintr_pts; ++pt) {
        double const wt = b.wt_intr(pt);
        for (int m = 0; m < b.nmodes; ++m) {
          double const phi = b.phi_intr(pt, m);
          double const mass = b.mass(m);
          for (int eq = 0; eq < neq; ++eq) {
            double const from_eq = interp_scalar_child_intr(
                from, b, from_cell, child, pt, eq);
            to(to_cell, eq, m) += from_eq * phi * wt / mass;
          }
        }
      }
    }
  };
  p3a::for_each(p3a::execution::par, generalize(from_subgrid), f);
}

void do_restriction(
    Basis const& b,
    View<double***> from,
    View<double***> to,
    p3a::grid3 const& from_grid,
    p3a::grid3 const& to_grid,
    p3a::subgrid3 const& from_subgrid,
    p3a::subgrid3 const& to_subgrid) {
  CALI_CXX_MARK_FUNCTION;
  verify_transfer(b, from, to, from_grid, to_grid, to_subgrid, from_subgrid);
  int const neq = from.extent(1);
  int const nchild = num_child(b.dim);
  int const nintr_pts = num_pts(b.dim, b.p);
  p3a::grid3 const general_from_grid = generalize(from_grid);
  p3a::grid3 const general_to_grid = generalize(to_grid);
  double const factor = std::pow(0.5, b.dim);
  auto f = [=] P3A_DEVICE (p3a::vector3<int> const& to_cell_ijk) {
    int const to_cell = general_to_grid.index(to_cell_ijk);
    p3a::vector3<int> const coarse_offset = to_cell_ijk - to_subgrid.lower();
    for (int child = 0; child < nchild; ++child) {
      p3a::vector3<int> const local = get_local(child);
      p3a::vector3<int> const fine_offset = get_fine_ijk(coarse_offset, local);
      p3a::vector3<int> const from_cell_ijk = from_subgrid.lower() + fine_offset;
      int const from_cell = general_from_grid.index(from_cell_ijk);
      for (int pt = 0; pt < nintr_pts; ++pt) {
        double const wt = b.wt_intr(pt);
        for (int m = 0; m < b.nmodes; ++m) {
          double const phi = b.phi_child_intr(child, pt, m);
          double const mass = b.mass(m);
          for (int eq = 0; eq < neq; ++eq) {
            double const from_eq = interp_scalar_intr(
                from, b, from_cell, pt, eq);
            to(to_cell, eq, m) += factor * from_eq * phi * wt / mass;
          }
        }
      }
    }
  };
  p3a::for_each(p3a::execution::par, generalize(to_subgrid), f);
}

static Tree copy_tree(
    Tree const& tree,
    std::vector<Node*> const& leaves) {
  Tree copy;
  copy.set_dim(tree.dim());
  copy.set_base(tree.base());
  for (Node* leaf : leaves) {
    copy.insert(leaf->pt());
  }
  return copy;
}

static Point get_adj_pt(Point const& pt, int axis, int dir) {
  Point adj_pt = pt;
  std::array<int, 2> offset = {-1,1};
  adj_pt.ijk[axis] += offset[dir];
  return adj_pt;
}

static int get_max_ijk(Point const& pt, Point const& base, int axis) {
  int const nblocks = base.ijk[axis];
  int const diff = pt.depth - base.depth;
  int const max_ijk = int(std::pow(2, diff)*nblocks) - 1;
  return max_ijk;
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

static std::vector<bool> get_2to1_refines(
    Tree& tree,
    p3a::vector3<bool> const& periodic,
    std::vector<Node*> const& leaves) {
  std::vector<bool> refines(leaves.size(), false);
  int const dim = tree.dim();
  Point const base = tree.base();
  for (Node* leaf : leaves) {
    Point const pt = leaf->pt();
    for (int axis = 0; axis < dim; ++axis) {
      for (int dir = 0; dir < ndirs; ++dir) {
        Point adj_pt = get_adj_pt(pt, axis, dir);
        bool const in = contains(base.depth, p3a::subgrid3(base.ijk), adj_pt);
        if (!in) {
          if (periodic[axis]) {
            adj_pt = make_adj_periodic(adj_pt, base, axis, dir);
          } else {
            continue;
          }
        }
        Node* adj = tree.find(adj_pt);
        if (!adj){
          Point const adj_parent_pt = get_parent_point(adj_pt);
          Node* adj_parent = tree.find(adj_parent_pt);
          if (!adj_parent) {
            Point const adj_parent_parent_pt = get_parent_point(adj_parent_pt);
            Node* adj_parent_parent = tree.find(adj_parent_parent_pt);
            verify_2to1_parent(adj_parent_parent);
            int const id = adj_parent_parent->block.id();
            refines[id] = true;
          }
        }
      }
    }
  }
  return refines;
}

void ensure_tree_is_2to1(
    mpicpp::comm* comm,
    Tree& tree,
    p3a::vector3<bool> const& periodic) {
  bool modified = true;
  while (modified) {
    modified = false;
    std::vector<Node*> const leaves = collect_leaves(tree);
    partition_leaves(comm, leaves); // this is to assign ids to leaves;
    std::vector<bool> const refines = get_2to1_refines(tree, periodic, leaves);
    for (size_t id = 0; id < refines.size(); ++id) {
      bool const should_refine = refines[id];
      if (should_refine) {
        Node* leaf = leaves[id];
        refine(tree.dim(), leaf);
        modified = true;
      }
    }
  }
}

std::vector<Point> collect_refine_pts(
    std::vector<int8_t> const& marks,
    std::vector<Node*> const& leaves) {
  verify_marks(marks, leaves);
  std::vector<Point> refine_pts;
  for (size_t i = 0; i < marks.size(); ++i) {
    if (marks[i] == REFINE) {
      Node* node = leaves[i];
      Point const pt = node->pt();
      refine_pts.push_back(pt);
    }
  }
  return refine_pts;
}

std::vector<Point> collect_coarsen_pts(
    std::vector<int8_t> const& marks,
    std::vector<Node*> const& leaves) {
  verify_marks(marks, leaves);
  std::vector<Point> coarsen_pts;
  for (size_t i = 0; i < marks.size(); ++i) {
    if (marks[i] == COARSEN) {
      Node* node = leaves[i];
      Point const pt = node->pt();
      coarsen_pts.push_back(pt);
    }
  }
  return coarsen_pts;
}

static void refine_tree(
    Tree& tree,
    std::vector<Point> const& refine_pts) {
  int const dim = tree.dim();
  for (Point const pt : refine_pts) {
    Node* node = tree.find(pt);
    refine(dim, node);
  }
}

static bool can_coarsen(int dim, Node* parent) {
  bool can = true;
  auto f = [&] (p3a::vector3<int> const& local) {
    Node* child = parent->child(local);
    if (!child->is_leaf()) can = false;
  };
  p3a::for_each(p3a::execution::seq, generalize(get_child_grid(dim)), f);
  return can;
}

static void coarsen_tree(
    Tree& tree,
    std::vector<Point> const& coarsen_pts) {
  int const dim = tree.dim();
  int const nchild = num_child(dim);
  std::list<Node*> unique_parents;
  std::list<Node*> sorted_parents;
  for (Point const pt : coarsen_pts) {
    Node* node = tree.find(pt);
    Node* parent = node->parent();
    unique_parents.push_back(parent);
  }
  unique_parents.sort();
  sorted_parents = unique_parents;
  unique_parents.unique();
  for (Node* parent :  unique_parents) {
    // we can only coarsen if every single child was marked for
    // coarsening AND if ensuring 2-1 interfaces at the refinement
    // step didn't make the children non-leaves. (Additionally, ensuring 2-1
    // interfaces after this step may just undo this coarsening as well)
    int const count = std::count(sorted_parents.begin(), sorted_parents.end(), parent);
    if ((count == nchild) && (can_coarsen(dim, parent))) {
      coarsen(dim, parent);
    }
  }
}

enum {ORIGINAL=0, MODIFIED=1};

struct Transfer {
  public:
    int op;
    Node* leaf[2];
    p3a::vector3<int> local = {0, 0, 0};
    Message<double***> msg;
  public:
    Transfer(int op_in, Node* original_leaf, Node* modified_leaf) {
      op = op_in;
      leaf[ORIGINAL] = original_leaf;
      leaf[MODIFIED] = modified_leaf;
      get_local();
    }
  private:
    void get_local() {
      if (op == REMAIN) return;
      int type = (op == REFINE) ? MODIFIED : ORIGINAL;
      Node* child = leaf[type];
      p3a::vector3<int> const child_ijk = child->pt().ijk;
      local = get_local_from_fine_ijk(child_ijk);
    }

};

struct Transfers {
  std::vector<Transfer> on_rank;
  std::vector<Transfer> send;
  std::vector<Transfer> recv;
};

template <class OP>
void collect(
    OP const& op,
    int mark,
    Node* original_leaf,
    Node* modified_leaf,
    std::vector<Transfer>& xfers) {
  verify_leaf(original_leaf);
  verify_leaf(modified_leaf);
  int const old_owner = original_leaf->block.owner();
  int const new_owner = modified_leaf->block.owner();
  if (op(old_owner, new_owner)) {
    xfers.push_back(Transfer(mark, original_leaf, modified_leaf));
  }
}

template <class OP>
void collect_children(
    OP const& op,
    int mark,
    int type,
    Node* leaf,
    Tree& tree,
    std::vector<Transfer>& xfers) {
  int const itype = invert_dir(type);
  int const dim = leaf->block.dim();
  Node* leaves[2] = {nullptr};
  leaves[type] = leaf;
  auto f = [&] (p3a::vector3<int> const& local) {
    Point const pt = get_child_point(leaf->pt(), local);
    leaves[itype] = tree.find(pt);
    collect(op, mark, leaves[ORIGINAL], leaves[MODIFIED], xfers);
  };
  p3a::for_each(p3a::execution::seq, generalize(get_child_grid(dim)), f);
}

template <class OP>
void collect_parent(
    OP const& op,
    int mark,
    int type,
    Node* leaf,
    Tree& tree,
    std::vector<Transfer>& xfers) {
  int const itype = invert_dir(type);
  Node* leaves[2] = {nullptr};
  leaves[type] = leaf;
  Point const pt = get_parent_point(leaf->pt());
  leaves[itype] = tree.find(pt);
  collect(op, mark, leaves[ORIGINAL], leaves[MODIFIED], xfers);
}

template <class OP>
void collect_same(
    OP const& op,
    int type,
    Node* leaf,
    Tree& tree,
    std::vector<Transfer>& xfers) {
  int const itype = invert_dir(type);
  Node* leaves[2] = {nullptr};
  leaves[type] = leaf;
  Point const pt = leaf->pt();
  leaves[itype] = tree.find(pt);
  collect(op, REMAIN, leaves[ORIGINAL], leaves[MODIFIED], xfers);
}

static void collect_on_rank(
    Mesh const& mesh,
    Tree& modified_tree,
    std::vector<Transfer>& xfers) {
  auto is_same_rank = [](int a, int b) { return a == b; };
  for (Node* leaf : mesh.owned_leaves()) {
    Point const pt = leaf->pt();
    Node* modified_node = modified_tree.find(pt);
    if (!modified_node) {
      collect_parent(is_same_rank, COARSEN, ORIGINAL, leaf, modified_tree, xfers);
    } else if (!(modified_node->is_leaf())) {
      collect_children(is_same_rank, REFINE, ORIGINAL, leaf, modified_tree, xfers);
    } else {
      collect_same(is_same_rank, ORIGINAL, leaf, modified_tree, xfers);
      continue;
    }
  }
}

static void collect_sends(
    Mesh const& mesh,
    Tree& modified_tree,
    std::vector<Transfer>& xfers) {
  auto is_diff_rank = [](int a, int b) { return a != b; };
  for (Node* leaf : mesh.owned_leaves()) {
    Point const pt = leaf->pt();
    Node* modified_node = modified_tree.find(pt);
    if (!modified_node) {
      collect_parent(is_diff_rank, COARSEN, ORIGINAL, leaf, modified_tree, xfers);
    } else if (!(modified_node->is_leaf())) {
      collect_children(is_diff_rank, REFINE, ORIGINAL, leaf, modified_tree, xfers);
    } else {
      collect_same(is_diff_rank, ORIGINAL, leaf, modified_tree, xfers);
    }
  }
}

static void collect_recvs(
    Mesh& mesh,
    std::vector<Node*> const& new_owned_leaves,
    std::vector<Transfer>& xfers) {
  auto is_diff_rank = [](int a, int b) { return a != b; };
  Tree& tree = mesh.tree();
  for (Node* leaf : new_owned_leaves) {
    Point const pt = leaf->pt();
    Node* original_node = tree.find(pt);
    if (!original_node) {
      collect_parent(is_diff_rank, REFINE, MODIFIED, leaf, tree, xfers);
    } else if (!original_node->is_leaf()) {
      collect_children(is_diff_rank, COARSEN, MODIFIED, leaf, tree, xfers);
    } else {
      collect_same(is_diff_rank, MODIFIED, leaf, tree, xfers);
    }
  }
}

static void allocate_block(
    Block& block,
    int nsoln,
    int nmodal_eq,
    int nflux_eq) {
  if (block.nsoln() == 0) {
    block.allocate(nsoln, nmodal_eq, nflux_eq);
  }
}

static void allocate_msg(
    Transfer& xfer,
    p3a::grid3 const& g,
    int nmodes,
    int nmodal_eq) {
  int ncells = 0;
  int const dim = xfer.leaf[ORIGINAL]->block.dim();
  int const nall_cells = generalize(g).size();
  int const nhalf_cells = nall_cells / ipow(2, dim);
  if (xfer.op == REMAIN) ncells = nall_cells;
  else ncells = nhalf_cells;
  Kokkos::resize(xfer.msg.val, ncells, nmodal_eq, nmodes);
}

static void allocate(Mesh const& mesh, Transfers& xfers) {
  Block const& block0 = mesh.owned_leaves()[0]->block;
  int const nsoln = block0.nsoln();
  int const nmodal_eq = mesh.nmodal_eq();
  int const nflux_eq = mesh.nflux_eq();
  int const nmodes = mesh.basis().nmodes;
  for (Transfer& xfer : xfers.on_rank) {
    allocate_block(xfer.leaf[MODIFIED]->block, nsoln, nmodal_eq, nflux_eq);
  }
  for (Transfer& xfer : xfers.send) {
    allocate_msg(xfer, mesh.cell_grid(), nmodes, nmodal_eq);
  }
  for (Transfer& xfer : xfers.recv) {
    allocate_block(xfer.leaf[MODIFIED]->block, nsoln, nmodal_eq, nflux_eq);
    allocate_msg(xfer, mesh.cell_grid(), nmodes, nmodal_eq);
  }
}

static void set_for_remain(Transfer& xfer, int type) {
  xfer.msg.val = xfer.leaf[type]->block.soln(0);
}

static void extract_for_refine(Transfer& xfer) {
  Block& parent = xfer.leaf[ORIGINAL]->block;
  View<double***> U_from = parent.soln(0);
  View<double***> U_to = xfer.msg.val;
  p3a::grid3 const from_grid = parent.cell_grid();
  p3a::subgrid3 const from_subgrid = get_local_subgrid(from_grid, xfer.local);
  p3a::grid3 const to_grid(from_subgrid.extents());
  p3a::subgrid3 const to_subgrid(to_grid);
  do_insertion(
      U_from, U_to,
      from_grid, to_grid,
      from_subgrid, to_subgrid);
}

static void extract_for_coarsen(Transfer& xfer) {
  Block& child = xfer.leaf[ORIGINAL]->block;
  View<double***> U_from = child.soln(0);
  View<double***> U_to = xfer.msg.val;
  p3a::grid3 const from_grid = child.cell_grid();
  p3a::subgrid3 const from_subgrid(from_grid);
  p3a::subgrid3 const local_subgrid = get_local_subgrid(from_grid, xfer.local);
  p3a::grid3 const to_grid(local_subgrid.extents());
  p3a::subgrid3 const to_subgrid(to_grid);
  do_restriction(
      child.basis(),
      U_from, U_to,
      from_grid, to_grid,
      from_subgrid, to_subgrid);
}

static void extract_msg_vals(Transfers& xfers) {
  for (Transfer& xfer : xfers.send) {
    if (xfer.op == REMAIN) set_for_remain(xfer, send);
    if (xfer.op == REFINE) extract_for_refine(xfer);
    if (xfer.op == COARSEN) extract_for_coarsen(xfer);
  }
  for (Transfer& xfer : xfers.recv) {
    if (xfer.op == REMAIN) set_for_remain(xfer, recv);
  }
}

static void begin_msg(
    mpicpp::comm* comm,
    Transfer& xfer,
    int type) {
  int const op = xfer.op;
  int const id = xfer.leaf[ORIGINAL]->block.id();
  int const itype = invert_dir(type);
  int const rank = xfer.leaf[itype]->block.owner();
  int const dim = xfer.leaf[ORIGINAL]->block.dim();
  int const nchild = num_child(dim);
  int const which_child = get_which_child(xfer.local);
  int const tag = (id*3 + op)*nchild + which_child;
  if (type == send) xfer.msg.send(comm, rank, tag);
  else xfer.msg.recv(comm, rank, tag);
}

static void begin_msgs(mpicpp::comm* comm, Transfers& xfers) {
  for (Transfer& xfer : xfers.send) {
    begin_msg(comm, xfer, send);
  }
  for (Transfer& xfer : xfers.recv) {
    begin_msg(comm, xfer, recv);
  }
}

static void copy_on_rank(Transfer& xfer) {
  // this is trivial copy for POD and assigns views (shared pointers)
  // to the new blocks as well
  xfer.leaf[MODIFIED]->block = xfer.leaf[ORIGINAL]->block;
}

static void prolong_on_rank(Transfer& xfer) {
  Block const& parent = xfer.leaf[ORIGINAL]->block;
  Block& child = xfer.leaf[MODIFIED]->block;
  p3a::vector3<int> const local = xfer.local;
  p3a::grid3 const from_grid = parent.cell_grid();
  p3a::grid3 const to_grid = child.cell_grid();
  p3a::subgrid3 const from_subgrid = get_local_subgrid(from_grid, local);
  p3a::subgrid3 const to_subgrid(to_grid);
  do_prolongation(
      parent.basis(),
      parent.soln(0), child.soln(0),
      from_grid, to_grid,
      from_subgrid, to_subgrid);
}

static void restrict_on_rank(Transfer& xfer) {
  Block& parent = xfer.leaf[MODIFIED]->block;
  Block const& child = xfer.leaf[ORIGINAL]->block;
  p3a::vector3<int> const local = xfer.local;
  p3a::grid3 const from_grid = child.cell_grid();
  p3a::grid3 const to_grid = parent.cell_grid();
  p3a::subgrid3 const from_subgrid(from_grid);
  p3a::subgrid3 const to_subgrid = get_local_subgrid(from_grid, local);
  do_restriction(
      child.basis(),
      child.soln(0), parent.soln(0),
      from_grid, to_grid,
      from_subgrid, to_subgrid);
}

static void transfer_on_rank(std::vector<Transfer>& xfers) {
  for (Transfer& xfer : xfers) {
    if (xfer.op == REMAIN) copy_on_rank(xfer);
    if (xfer.op == REFINE) prolong_on_rank(xfer);
    if (xfer.op == COARSEN) restrict_on_rank(xfer);
  }
}

static void end_msgs(Transfers& xfers) {
  for (Transfer& xfer : xfers.send) {
    xfer.msg.wait();
  }
  for (Transfer& xfer : xfers.recv) {
    xfer.msg.wait();
  }
}

static void inject_for_refine(Transfer& xfer) {
  Block& child = xfer.leaf[MODIFIED]->block;
  View<double***> U_from = xfer.msg.val;
  View<double***> U_to = child.soln(0);
  p3a::grid3 const to_grid = child.cell_grid();
  p3a::subgrid3 const to_subgrid(to_grid);
  p3a::subgrid3 const local_subgrid = get_local_subgrid(to_grid, xfer.local);
  p3a::grid3 const from_grid(local_subgrid.extents());
  p3a::subgrid3 const from_subgrid(from_grid);
  do_prolongation(
      child.basis(),
      U_from, U_to,
      from_grid, to_grid,
      from_subgrid, to_subgrid);
}

static void inject_for_coarsen(Transfer& xfer) {
  Block& parent = xfer.leaf[MODIFIED]->block;
  View<double***> U_from = xfer.msg.val;
  View<double***> U_to = parent.soln(0);
  p3a::grid3 const to_grid = parent.cell_grid();
  p3a::subgrid3 const to_subgrid = get_local_subgrid(to_grid, xfer.local);
  p3a::grid3 const from_grid(to_subgrid.extents());
  p3a::subgrid3 const from_subgrid(from_grid);
  do_insertion(
      U_from, U_to,
      from_grid, to_grid,
      from_subgrid, to_subgrid);
}

static void inject_msg_vals(std::vector<Transfer>& xfers) {
  for (Transfer& xfer : xfers) {
    if (xfer.op == REFINE) inject_for_refine(xfer);
    if (xfer.op == COARSEN) inject_for_coarsen(xfer);
  }
}

static void transfer_data(
    Mesh& mesh,
    std::vector<Node*> const& new_owned_leaves,
    Tree& modified_tree) {
  Transfers xfers;
  collect_on_rank(mesh, modified_tree, xfers.on_rank);
  collect_sends(mesh, modified_tree, xfers.send);
  collect_recvs(mesh, new_owned_leaves, xfers.recv);
  allocate(mesh, xfers);
  extract_msg_vals(xfers);
  p3a::execution::par.synchronize();
  begin_msgs(mesh.comm(), xfers);
  transfer_on_rank(xfers.on_rank);
  end_msgs(xfers);
  inject_msg_vals(xfers.recv);
}

static void cleanup(Mesh& mesh) {
  for (Node* leaf : mesh.owned_leaves()) {
    Block& block = leaf->block;
    Kokkos::deep_copy(block.resid(), 0.);
    int const nmodal_eq = mesh.nmodal_eq();
    int const nflux_eq = mesh.nflux_eq();
    for (int axis = 0; axis < mesh.dim(); ++axis) {
      Kokkos::deep_copy(block.flux(axis), 0.);
      for (int dir = 0; dir < ndirs; ++dir) {
        Border& border = block.border(axis, dir);
        border.deallocate();
        border.allocate(nmodal_eq, nflux_eq);
      }
    }
  }
}

void modify(Mesh& mesh, std::vector<int8_t> const& in_marks) {
  CALI_CXX_MARK_FUNCTION;
  verify_marks(mesh, in_marks);
  verify_mesh(mesh);
  mpicpp::comm* comm = mesh.comm();
  std::vector<int8_t> marks = reduce_marks(mesh, in_marks);
  Tree copy = copy_tree(mesh.tree(), mesh.leaves());
  std::vector<Node*> leaves = collect_leaves(copy);
  std::vector<Point> refine_pts = collect_refine_pts(marks, leaves);
  std::vector<Point> coarsen_pts = collect_coarsen_pts(marks, leaves);
  refine_tree(copy, refine_pts);
  ensure_tree_is_2to1(comm, copy, mesh.periodic());
  coarsen_tree(copy, coarsen_pts);
  ensure_tree_is_2to1(comm, copy, mesh.periodic());
  leaves = collect_leaves(copy);
  partition_leaves(comm, leaves);
  init_leaves(&mesh, copy, mesh.periodic(), leaves);
  std::vector<Node*> owned_leaves = collect_owned_leaves(comm, leaves);
  transfer_data(mesh, owned_leaves, copy);
  mesh.set_tree(copy);
  mesh.rebuild();
  mesh.clean();
  cleanup(mesh);
}

}
