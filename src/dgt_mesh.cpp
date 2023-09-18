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
  m_cell_grid.extents() = Vec3<int>::zero();
  int const dim = infer_dimension(cell_grid);
  for (int axis = 0; axis < dim; ++axis) {
    m_cell_grid.extents()[axis] = cell_grid.extents()[axis] + 2;
  }
}

void Mesh::set_periodic(Vec3<bool> const& periodic)
{
  m_periodic = periodic;
}

void Mesh::set_basis(int const p, int const q, bool const tensor)
{
  m_basis = build_basis<View>(dim(), p, q, tensor);
  m_basis_h = build_basis<HostView>(dim(), p, q, tensor);
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
  m_block_info = build_block_info<View>(
      dim, m_domain, m_owned_leaves, base_pt);
  m_block_info_h = build_block_info<HostView>(
      dim, m_domain, m_owned_leaves, base_pt);
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
      std::string const fname = fmt::format("{}_fluxes_{}", name, axis_name);
      f.fluxes[axis].create(fname, nblocks, nfaces, nface_pts, ncomps);
    }
  }
  m_modal.push_back(f);
  //TODO also initialize the communication pattern here
}

int Mesh::dim() const
{
  return infer_dimension(m_cell_grid);
}

int Mesh::num_total_blocks() const
{
  return int(m_leaves.size());
}

int Mesh::num_owned_blocks() const
{
  return int(m_owned_leaves.size());
}

int Mesh::num_total_cells() const
{
  Subgrid3 const owned_cells = get_owned_cells(m_cell_grid);
  Grid3 const owned_cell_grid(owned_cells.extents());
  int const num_cells_per_block = get_num_cells(owned_cell_grid);
  return num_total_blocks() * num_cells_per_block;
}

int Mesh::num_owned_cells() const
{
  Subgrid3 const owned_cells = get_owned_cells(m_cell_grid);
  Grid3 const owned_cell_grid(owned_cells.extents());
  int const num_cells_per_block = get_num_cells(owned_cell_grid);
  return num_owned_blocks() * num_cells_per_block;
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

Field<real***> Mesh::get_fluxes(std::string const& name, int const axis)
{
  int const modal_idx = modal_index(name);
  if (modal_idx < 0) {
    std::string const msg = fmt::format(
        "dgt::Mesh::get_fluxes -> field {} doesn't exist", name);
  }
  return m_modal[modal_idx].fluxes[axis];
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

static Vec3<real> vec_min(Vec3<real> const& a, Vec3<real> const b)
{
  return (a.x() < b.x()) ? a : b;
}

static Vec3<real> vec_max(Vec3<real> const& a, Vec3<real> const& b)
{
  return (a.x() > b.x()) ? a : b;
}

template <class Op>
static Vec3<real> const reduce_dx(
    mpicpp::comm* comm,
    HostView<Vec3<real>*> dxs,
    Op const& op)
{
  Vec3<real> result = dxs[0];
  for (std::size_t i = 0; i < dxs.size(); ++i) {
    result = op(result, dxs[i]);
  }
  mpicpp::request req = comm->iallreduce(
      &(result[0]), DIMENSIONS, mpicpp::op::sum());
  req.wait();
  return result;
}

void Mesh::print_stats() const
{
  if (m_comm->rank()) return;
  Vec3<real> const min = reduce_dx(m_comm, m_block_info_h.dxs, vec_min);
  Vec3<real> const max = reduce_dx(m_comm, m_block_info_h.dxs, vec_max);
  printf("mesh stats\n");
  printf("----------\n");
  printf("> blocks: %d\n", num_total_blocks());
  printf("> cells: %d\n", num_total_cells());
  printf("> minimum dx: [%e, %e, %e]\n", min.x(), min.y(), min.z());
  printf("> maximum dx: [%e, %e, %e]\n", max.x(), max.y(), max.z());
}

}
