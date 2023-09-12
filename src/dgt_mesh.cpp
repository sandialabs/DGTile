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

static void verify_basis(Basis<View> const& basis)
{
  if ((basis.dim < 1) || (basis.dim > 3)) {
    throw std::runtime_error("dgt::Mesh - invalid basis.dim");
  }
  if ((basis.p < 0) || (basis.p > max_polynomial_order)) {
    throw std::runtime_error("dgt::Mesh - invalid basis.p");
  }
  if ((basis.q < 1) || (basis.q > max_1D_quadrature_points)) {
    throw std::runtime_error("dgt::Mesh - invalid basis.q");
  }
}

static void verify_dimensions(Grid3 const& cell_grid, Basis<View> const& basis)
{
  if (basis.dim != infer_dimension(cell_grid)) {
    throw std::runtime_error("dgt::Mesh - cell_grid/basis dimension mismatch");
  }
}

void Mesh::set_comm(mpicpp::comm* comm)
{
  m_comm = comm;
}

void Mesh::set_domain(Box3<real> const& domain)
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

void Mesh::set_basis(Basis<View> basis)
{
  m_basis = basis;
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
  auto end = z_leaves.begin() + offset + nlocal;
  return std::vector<tree::ID>(begin, end);
}

void Mesh::verify()
{
  verify_comm(m_comm);
  verify_cell_grid(m_cell_grid);
  verify_basis(m_basis);
  verify_dimensions(m_cell_grid, m_basis);
}

void Mesh::init(Grid3 const& block_grid)
{
  verify();
  int const dim = infer_dimension(m_cell_grid);
  m_leaves = tree::create(dim, block_grid);
  m_zleaves = tree::order(dim, m_leaves);
  m_owned_leaves = get_owned_leaves(m_comm, m_zleaves);
  tree::Point const base_pt = tree::get_base_point(dim, m_zleaves);
  m_block_info = create_block_info(dim, m_domain, m_owned_leaves, base_pt);
}

static bool is_modal(std::string const& name, std::vector<Modal> const& fields)
{
  for (auto f : fields) {
    if (f.name == name) return true;
  }
  return false;
}

static void verify_doesnt_exist(
    std::string const& name,
    std::vector<Modal> const& fields)
{
  if (is_modal(name, fields)) {
    std::string const msg = fmt::format(
        "dgt::Mesh::add_modal - field {} exists", name);
    throw std::runtime_error(msg);
  }
}

void Mesh::add_modal(std::string const& name, bool const with_flux)
{
  verify_doesnt_exist(name, m_fields);
  Modal f;
  f.name = name;
  f.has_flux = with_flux;
  m_fields.push_back(f);
}

}
