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

static void verify_modal(std::vector<ModalDescriptor> const& modal)
{
  if (modal.size() == 0) {
    throw std::runtime_error("dgt::Mesh - no modal fields added");
  }
  for (ModalDescriptor const& m : modal) {
    if (m.name == "") {
      throw std::runtime_error("dgt::Mesh - modal field unnamed");
    }
    std::string const base = fmt::format("dgt::Mesh - modal_field[{}]", m.name);
    if (m.num_stored < 0) {
      throw std::runtime_error(base + ".num_stored < 0");
    }
    if (m.num_comps < 0) {
      throw std::runtime_error(base + ".num_comps < 0");
    }
  }
}

static void verify_fields(std::vector<FieldDescriptor> const& fields)
{
  for (FieldDescriptor const& f : fields) {
    if (f.name == "") {
      throw std::runtime_error("dgt::Mesh - field unnamed");
    }
    std::string const base = fmt::format("dgt::Mesh - field[{}]", f.name);
    if (f.entity_dim < 0) {
      throw std::runtime_error(base + ".num_stored < 0");
    }
    if (f.num_comps < 0) {
      throw std::runtime_error(base + ".num_comps < 0");
    }
  }
}

void Mesh::ensure_set()
{
  verify_comm(m_comm);
  verify_cell_grid(m_cell_grid);
  verify_basis(m_basis);
  verify_dimensions(m_cell_grid, m_basis);
}

void Mesh::ensure_added()
{
  verify_modal(m_modal_meta);
  verify_fields(m_field_meta);
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

template <class VecT>
static int index(VecT const& vector, std::string const& name)
{
  for (std::size_t i = 0; i < vector.size(); ++i) {
    if (vector[i].name == name) return int(i);
  }
  return -1;
}

void Mesh::add_modal(ModalDescriptor const modal)
{
  if (index(m_modal_meta, modal.name) >= 0) {
    std::string const msg = fmt::format(
        "dgt::Mesh::add_modal -> field {} already exists", modal.name);
    throw std::runtime_error(msg);
  }
  m_modal_meta.push_back(modal);
}

void Mesh::add_field(FieldDescriptor const field)
{
  if (index(m_field_meta, field.name) >= 0) {
    std::string const msg = fmt::format(
        "dgt::Mesh::add_field -> field {} already exists", field.name);
    throw std::runtime_error(msg);
  }
  m_field_meta.push_back(field);
}

static tree::Partitioning get_partitioning(
    mpicpp::comm* comm,
    tree::ZLeaves const& z_leaves)
{
  using namespace linear_partitioning;
  int const nranks = comm->size();
  int const nleaves = z_leaves.size();
  tree::Partitioning partitioning;
  for (int rank = 0; rank < nranks; ++rank) {
    int const nlocal = get_num_local(nleaves, nranks, rank);
    int const offset = get_local_offset(nleaves, nranks, rank);
    for (int block = 0; block < nlocal; ++block) {
      tree::ID const leaf_id = z_leaves[offset + block];
      partitioning[leaf_id].rank = rank;
      partitioning[leaf_id].block = block;
    }
  }
  return partitioning;
}

static tree::ZLeaves get_owned_leaves(
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

void Mesh::initialize(Grid3 const& block_grid)
{
  ensure_set();
  ensure_added();
  m_leaves = tree::create(infer_dimension(m_cell_grid), block_grid);
  initialize();
}

void Mesh::initialize(tree::Leaves const& leaves)
{
  ensure_set();
  ensure_added();
  m_leaves = leaves;
  initialize();
}

static ModalField make_modal_field(
    ModalDescriptor const& m,
    Grid3 const& cell_grid,
    Basis<View> const& B,
    int const nblocks)
{
  int const dim = B.dim;
  int const nmodes = B.num_modes;
  int const nface_pts = B.num_face_pts;
  int const ncells = get_num_cells(cell_grid);
  ModalField f;
  f.solution.resize(m.num_stored);
  std::string const rname = fmt::format("{}_residual", m.name);
  f.residual.create(rname, nblocks, ncells, m.num_comps, nmodes);
  for (int soln_idx = 0; soln_idx < m.num_stored; ++soln_idx) {
    std::string const sname = fmt::format("{}_{}", m.name, soln_idx);
    f.solution[soln_idx].create(sname, nblocks, ncells, m.num_comps, nmodes);
  }
  if (m.with_flux) {
    for (int axis = 0; axis < dim; ++axis) {
      int const nfaces = get_num_faces(cell_grid, axis);
      std::string const axis_name = get_axis_name(axis);
      std::string const fname = fmt::format("{}_fluxes_{}", m.name, axis_name);
      f.fluxes[axis].create(fname, nblocks, nfaces, nface_pts, m.num_comps);
    }
  }
  return f;
}

void Mesh::initialize()
{
  int const dim = infer_dimension(m_cell_grid);
  m_zleaves = tree::order(dim, m_leaves);
  m_partitioning = get_partitioning(m_comm, m_zleaves);
  m_owned_leaves = get_owned_leaves(m_comm, m_zleaves);
  tree::Point const base_pt = tree::get_base_point(dim, m_zleaves);
  m_owned_adjacencies = tree::get_adjacencies(
      dim, m_owned_leaves, m_leaves, base_pt, m_periodic);
  m_block_info = build_block_info<View>(
      dim, m_cell_grid, m_domain, m_owned_leaves, base_pt);
  m_block_info_h = build_block_info<HostView>(
      dim, m_cell_grid, m_domain, m_owned_leaves, base_pt);
  int const nblocks = int(m_owned_leaves.size());
  for (ModalDescriptor const& meta : m_modal_meta) {
    ModalField f = make_modal_field(meta, m_cell_grid, m_basis, nblocks);
    m_modal.push_back(f);
  }
  for (FieldDescriptor const& f : m_field_meta) {
    (void)f;
    //TODO: Field f = build_field(f, m_cell_grid, nblocks);
    //TODO: m_field.push_back(f);
  }
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

Field<real***>& Mesh::get_solution(std::string const& name, int const soln_idx)
{
  int const modal_idx = index(m_modal_meta, name);
  if (modal_idx < 0) {
    std::string const msg = fmt::format(
        "dgt::Mesh::get_solution -> field {} doesn't exist", name);
    throw std::runtime_error(msg);
  }
  return m_modal[modal_idx].solution[soln_idx];
}

Field<real***>& Mesh::get_fluxes(std::string const& name, int const axis)
{
  int const modal_idx = index(m_modal_meta, name);
  if (modal_idx < 0) {
    std::string const msg = fmt::format(
        "dgt::Mesh::get_fluxes -> field {} doesn't exist", name);
  }
  return m_modal[modal_idx].fluxes[axis];
}

Field<real***>& Mesh::get_residual(std::string const& name)
{
  int const modal_idx = index(m_modal_meta, name);
  if (modal_idx < 0) {
    std::string const msg = fmt::format(
        "dgt::Mesh::get_residual -> field {} doesn't exist", name);
  }
  return m_modal[modal_idx].residual;
}

Field<real***> const& Mesh::get_solution(
    std::string const& name,
    int const soln_idx) const
{
  int const modal_idx = index(m_modal_meta, name);
  if (modal_idx < 0) {
    std::string const msg = fmt::format(
        "dgt::Mesh::get_solution -> field {} doesn't exist", name);
    throw std::runtime_error(msg);
  }
  return m_modal[modal_idx].solution[soln_idx];
}

Field<real***> const& Mesh::get_fluxes(
    std::string const& name,
    int const axis) const
{
  int const modal_idx = index(m_modal_meta, name);
  if (modal_idx < 0) {
    std::string const msg = fmt::format(
        "dgt::Mesh::get_fluxes -> field {} doesn't exist", name);
  }
  return m_modal[modal_idx].fluxes[axis];
}

Field<real***> const& Mesh::get_residual(std::string const& name) const
{
  int const modal_idx = index(m_modal_meta, name);
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
  Vec3<real> const min = reduce_dx(m_comm, m_block_info_h.cell_dxs, vec_min);
  Vec3<real> const max = reduce_dx(m_comm, m_block_info_h.cell_dxs, vec_max);
  printf("mesh stats:\n");
  printf("> blocks: %d\n", num_total_blocks());
  printf("> cells: %d\n", num_total_cells());
  printf("> minimum cell dx: [%e, %e, %e]\n", min.x(), min.y(), min.z());
  printf("> maximum cell dx: [%e, %e, %e]\n", max.x(), max.y(), max.z());
}

}
