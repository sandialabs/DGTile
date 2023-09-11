#include "dgt_mesh.hpp"

namespace dgt {

static void verify_comm(mpicpp::comm* comm)
{
  if (!comm) {
    throw std::runtime_error("dgt::Mesh - invalid comm");
  }
}

static void verify_cell_grid(Grid3 const& cell_grid)
{
  int const dim = infer_dimension(cell_grid);
  if (dim < 0) {
    throw std::runtime_error("dgt::Mesh - invalid cell grid");
  }
}

void Mesh::set_comm(mpicpp::comm* comm)
{
  m_comm = comm;
}

void Mesh::set_domain(Vec3<real> const& domain)
{
  m_domain = domain;
}

void Mesh::set_cell_grid(Grid3 const& cell_grid)
{
  m_cell_grid = cell_grid;
}

void Mesh::set_periodic(Vec3<bool> const& periodic)
{
  m_periodic = periodic;
}

void Mesh::set_basis(int const p, int const q, bool const tensor)
{
  verify_cell_grid(m_cell_grid);
  int const dim = infer_dimension(m_cell_grid);
  m_basis = build_basis<View>(dim, p, q, tensor);
}

void Mesh::init(Grid3 const& block_grid)
{
  verify_comm(m_comm);
  verify_cell_grid(m_cell_grid);
  int const dim = infer_dimension(m_cell_grid);
  m_leaves = tree::create(dim, block_grid);
  m_zleaves = tree::order(dim, m_leaves);
}

}
