#include "dgt_mesh.hpp"
#include "dgt_partitioning.hpp"

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


std::vector<tree::ID> get_owned_leaves(
    mpicpp::comm* comm,
    tree::ZLeaves const& z_leaves)
{
  using namespace linear_partitioning;
  int const rank = comm->rank();
  int const nranks = comm->size();
  int const nleaves = z_leaves.size();
  int const nlocal = get_num_local(nleaves, nranks, rank);
  int const offset = get_local_offset(nleaves, nranks, rank);
  auto begin = z_leaves.begin() + offset;
  auto end = z_leaves.end() + offset + nlocal;
  return std::vector<tree::ID>(begin, end);
}

void Mesh::init(Grid3 const& block_grid)
{
  verify_comm(m_comm);
  verify_cell_grid(m_cell_grid);
  int const dim = infer_dimension(m_cell_grid);
  m_leaves = tree::create(dim, block_grid);
  m_zleaves = tree::order(dim, m_leaves);
  m_owned_leaves = get_owned_leaves(m_comm, m_zleaves);
  tree::Point const base_pt = tree::get_base_point(dim, m_zleaves);
  m_block_info = create_block_info(dim, m_domain, m_owned_leaves, base_pt);
}

}
