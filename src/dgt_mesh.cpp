#include "dgt_cartesian.hpp"
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

tree::OwnedLeaves get_owned_leaves(
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

int Mesh::modal_index(std::string const& name)
{
  for (std::size_t i = 0; i < m_modal.size(); ++i) {
    if (m_modal[i].name == name) return int(i);
  }
  return -1;
}

void Mesh::add_modal(
    std::string const& name,
    int const nstored,
    int const ncomps,
    bool const with_flux)
{
  verify();
  if (modal_index(name) >= 0) {
    std::string const msg = fmt::format(
        "dgt::Mesh::add_modal -> field {} already exists", name);
    throw std::runtime_error(msg);
  }
  int const dim = m_basis.dim;
  int const nmodes = m_basis.num_modes;
  int const nface_pts = m_basis.num_face_pts;
  int const ncells = get_num_cells(m_cell_grid);
  int const nblocks = int(m_owned_leaves.size());
  SolutionField f;
  f.name = name;
  f.solution.resize(nstored);
  std::string const rname = fmt::format("{}_residual", name);
  f.residual.create(rname, nblocks, ncells, ncomps, nmodes);
  for (int soln_idx = 0; soln_idx < nstored; ++soln_idx) {
    std::string const sname = fmt::format("{}_{}", name, soln_idx);
    f.solution[soln_idx].create(sname, nblocks, ncells, ncomps, nmodes);
  }
  if (with_flux) {
    for (int axis = 0; axis < dim; ++axis) {
      int const nfaces = get_num_faces(m_cell_grid, axis);
      std::string const axis_name = get_axis_name(axis);
      std::string const fname = fmt::format("{}_flux_{}", name, axis_name);
      f.fluxes[axis].create(fname, nblocks, nfaces, nface_pts, ncomps);
    }
  }
  m_modal.push_back(f);
  //TODO also initialize the communication pattern here
}

Field<real***> Mesh::get_solution(std::string const& name, int const soln_idx)
{
  int const modal_idx = modal_index(name);
  if (modal_idx < 0) {
    std::string const msg = fmt::format(
        "dgt::Mesh::get_solution -> field {} doesn't exist", name);
    throw std::runtime_error(msg);
  }
  return m_modal[modal_idx].solution[soln_idx];
}

Vec3<Field<real***>> Mesh::get_flux(std::string const& name)
{
  int const modal_idx = modal_index(name);
  if (modal_idx < 0) {
    std::string const msg = fmt::format(
        "dgt::Mesh::get_flux -> field {} doesn't exist", name);
  }
  return m_modal[modal_idx].fluxes;
}

Field<real***> Mesh::get_residual(std::string const& name)
{
  int const modal_idx = modal_index(name);
  if (modal_idx < 0) {
    std::string const msg = fmt::format(
        "dgt::Mesh::get_residual -> field {} doesn't exist", name);
  }
  return m_modal[modal_idx].residual;
}

}
