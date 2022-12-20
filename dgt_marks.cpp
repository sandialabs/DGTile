#include "dgt_amr.hpp"
#include "dgt_marks.hpp"
#include "dgt_mesh.hpp"
#include "dgt_spatial.hpp"

namespace dgt {

static void verify_marks(Mesh const& mesh, std::vector<int8_t> const& marks) {
  if (marks.size() != mesh.leaves().size()) {
    throw std::runtime_error("marks- invalid marks");
  }
}

static void verify_leaf(Node const* leaf) {
  if (!leaf) {
    throw std::runtime_error("marks - leaf doesn't exist");
  }
  if (!leaf->is_leaf()) {
    throw std::runtime_error("marks - not a leaf");
  }
}

std::vector<int8_t> reduce_marks(
    Mesh const& mesh,
    std::vector<int8_t> const& in_marks) {
  verify_marks(mesh, in_marks);
  mpicpp::comm* comm = mesh.comm();
  std::vector<int8_t> marks = in_marks;
  MPI_Allreduce(
      in_marks.data(), marks.data(), marks.size(),
      MPI_BYTE, MPI_MAX, comm->get());
  return marks;
}

static bool ensure_2to1_refine(
    Border const& border,
    std::vector<int8_t>& marks) {
  Block const& block = border.node()->block;
  int const id = block.id();
  int const dim = block.dim();
  int const axis = border.axis();
  int const dir = border.dir();
  int const nchild = num_child(dim-1);
  bool changed = false;
  for (int which_child = 0; which_child < nchild; ++which_child) {
    vector3<int> local = get_local(axis, which_child);
    local[axis] = invert_dir(dir);
    Node const* adj_node = border.adj()->child(local);
    verify_leaf(adj_node);
    int const adj_id = adj_node->block.id();
    if (marks[adj_id] == REFINE) {
      if (marks[id] != REFINE) {
        marks[id] = REFINE;
        changed = true;
      }
    }
  }
  return changed;
}

static void ensure_2to1_refine(
    std::vector<Node*> const& leaves,
    std::vector<int8_t>& marks) {
  bool changed = true;
  while (changed) {
    changed = false;
    for (Node* leaf : leaves) {
      int const dim = leaf->block.dim();
      for (int axis = 0; axis < dim; ++axis) {
        for (int dir = 0; dir < ndirs; ++dir) {
          Border const& border = leaf->block.border(axis, dir);
          if (border.type() == COARSE_TO_FINE) {
            bool const flip = ensure_2to1_refine(border, marks);
            if (flip) changed = true;
          }
        }
      }
    }
  }
}

static void ensure_2to1_refine_coarsen(
    std::vector<Node*> const& leaves,
    std::vector<int8_t>& marks) {
  for (Node* leaf : leaves) {
    Block const& block = leaf->block;
    int const dim = block.dim();
    int const id = block.id();
    if (marks[id] == COARSEN) {
      for (int axis = 0; axis < dim; ++axis) {
        for (int dir = 0; dir < ndirs; ++dir) {
          Border const& border = block.border(axis, dir);
          if (border.type() == STANDARD) {
            int const adj_id = border.adj()->block.id();
            if (marks[adj_id] == REFINE) {
              marks[id] = REMAIN;
            }
          }
        }
      }
    }
  }
}

static void ensure_2to1_coarsen(
    Tree const& tree,
    std::vector<Node*> const& leaves,
    std::vector<int8_t>& marks) {
  int const base_depth = tree.base().depth;
  for (Node* leaf : leaves) {
    Block const& block = leaf->block;
    int const dim = block.dim();
    int const id = block.id();
    if (marks[id] == COARSEN) {
      if (leaf->pt().depth == base_depth) {
        marks[id] = REMAIN;
        continue;
      }
      for (int axis = 0; axis < dim; ++axis) {
        for (int dir = 0; dir < ndirs; ++dir) {
          Border const& border = block.border(axis, dir);
          if (border.type() == COARSE_TO_FINE) {
            marks[id] = REMAIN;
          }
        }
      }
      Node* parent = leaf->parent();
      for (int which_child = 0; which_child < num_child(dim); ++which_child) {
        vector3<int> const local = get_local(which_child);
        Node* child = parent->child(local);
        if (!child->is_leaf()) {
          marks[id] = REMAIN;
          continue;
        }
        int const child_id = child->block.id();
        if (marks[child_id] != COARSEN) {
          marks[id] = REMAIN;
        }
      }
    }
  }
}

void ensure_2to1(
    Mesh const& mesh,
    std::vector<int8_t>& marks) {
  Tree const& tree = mesh.tree();
  std::vector<Node*> const& leaves = mesh.leaves();
  ensure_2to1_refine(leaves, marks);
  ensure_2to1_refine_coarsen(leaves, marks);
  ensure_2to1_coarsen(tree, leaves, marks);
}

}
