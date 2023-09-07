#include "dgt_mesh.hpp"

namespace dgt {

void Mesh::set_comm(mpicpp::comm* comm)
{
  m_comm = comm;
}

void Mesh::set_cell_grid(Grid3 const& cell_grid)
{
  m_cell_grid = cell_grid;
  int const dim = infer_dimension(cell_grid);
  for (int axis = 0; axis < dim; ++axis) {
    m_cell_grid.extents()[axis] += 2;
  }
}

void Mesh::set_periodic(Vec3<bool> const& periodic)
{
  m_periodic = periodic;
}

void Mesh::set_domain(Box3<real> const& domain)
{
  m_domain = domain;
}

void Mesh::set_basis(int const dim, int const p, int const q, bool const tensor)
{
  m_basis = build_basis<View>(dim, p, q, tensor);
}

void Mesh::set_tree(tree::Leaves const& leaves)
{
  m_leaves = leaves;
}

}
