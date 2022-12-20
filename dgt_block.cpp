#include <stdexcept>

#include "caliper/cali.h"

#include "dgt_block.hpp"
#include "dgt_grid.hpp"
#include "dgt_mesh.hpp"
#include "dgt_spatial.hpp"

namespace dgt {

static void verify_mesh(Mesh const* m) {
  if (!m) {
    throw std::runtime_error("Block - unset mesh");
  }
}

static void verify_node(Node const* n) {
  if (!n) {
    throw std::runtime_error("Block - unset node");
  }
}

static void verify_axis(int dim, int axis) {
  if ((axis < 0) || (axis >= dim)) {
    throw std::runtime_error("Block - invalid axis");
  }
}

static void verify_dir(int dir) {
  if ((dir != left) && (dir != right)) {
    throw std::runtime_error("Block - invalid dir");
  }
}

static void verify_basis(Basis const& b) {
  if ((b.dim == -1) || (b.p == -1)) {
    throw std::runtime_error("Block - unset basis");
  }
}

static void verify_U_idx(int nsoln, int idx) {
  if ((idx < 0) || (idx >= nsoln)) {
    throw std::runtime_error("Block - invalid U idx");
  }
}

static void verify_no_field(
    std::string const& name,
    std::vector<Field> const& fields) {
  for (Field const& f : fields) {
    if (f.name() == name) {
      throw std::runtime_error("Block - field " + name + " exists");
    }
  }
}

static void verify_field_idx(int idx, std::string const& name) {
  if (idx == -1) {
    throw std::runtime_error("Block - field " + name + " doesn't exist");
  }
}

int Block::id() const {
  return m_id;
}

int Block::owner() const {
  return m_owner;
}

int Block::dim() const {
  verify_mesh(m_mesh);
  return m_mesh->dim();
}

int Block::nsoln() const {
  return m_soln.size();
}

int Block::nfields() const {
  return m_fields.size();
}

p3a::grid3 Block::cell_grid() const {
  verify_mesh(m_mesh);
  return m_mesh->cell_grid();
}

p3a::box3<double> Block::domain() const {
  verify_mesh(m_mesh);
  verify_node(m_node);
  Point const base_pt = m_mesh->tree().base();
  Point const node_pt = m_node->pt();
  p3a::box3<double> const box = m_mesh->domain();
  return get_block_domain(base_pt, node_pt, box);
}

p3a::vector3<double> Block::dx() const {
  return get_dx(domain(), cell_grid());
}

double Block::cell_detJ() const {
  return get_cell_detJ(dim(), dx());
}

double Block::side_detJ(int axis) const {
  return get_side_detJ(dim(), axis, dx());
}

double Block::amr_side_detJ(int axis) const {
  return get_amr_side_detJ(dim(), axis, dx());
}

Mesh const* Block::mesh() const {
  return m_mesh;
}

Node const* Block::node() const {
  return m_node;
}

Basis const& Block::basis() const {
  verify_mesh(m_mesh);
  verify_basis(m_mesh->basis());
  return m_mesh->basis();
}

Border& Block::border(int axis, int dir) {
  verify_axis(dim(), axis);
  verify_dir(dir);
  return m_border[axis][dir];
}

Border const& Block::border(int axis, int dir) const {
  verify_axis(dim(), axis);
  verify_dir(dir);
  return m_border[axis][dir];
}

View<double***> Block::soln(int idx) const {
  verify_U_idx(m_soln.size(), idx);
  return m_soln[idx];
}

View<double***> Block::flux(int axis) const {
  verify_axis(dim(), axis);
  return m_flux[axis];
}

View<double***> Block::resid() const {
  return m_resid;
}

p3a::simd_view<double***> Block::simd_soln(int idx) const {
  verify_U_idx(m_soln.size(), idx);
  return p3a::simd_view<double***>(m_soln[idx]);
}

p3a::simd_view<double***> Block::simd_flux(int axis) const {
  verify_axis(dim(), axis);
  return p3a::simd_view<double***>(m_flux[axis]);
}

p3a::simd_view<double***> Block::simd_resid() const {
  return p3a::simd_view<double***>(m_resid);
}

int Block::field_idx(std::string name) const {
  for (int idx = 0; idx < int(m_fields.size()); ++idx) {
    if (m_fields[idx].name() == name) {
      return idx;
    }
  }
  return -1;
}

Field const& Block::field(std::string name) const {
  int const idx = field_idx(name);
  verify_field_idx(idx, name);
  return m_fields[idx];
}

void Block::set_id(int id) {
  m_id = id;
}

void Block::set_owner(int owner) {
  m_owner = owner;
}

void Block::set_mesh(Mesh* mesh) {
  m_mesh = mesh;
}

void Block::set_node(Node* node) {
  m_node = node;
}

void Block::add_field(FieldInfo const& info) {
  verify_no_field(info.name, m_fields);
  Field field;
  field.set_info(info);
  m_fields.push_back(field);
}

void Block::reset() {
  m_id = -1;
  m_owner = -1;
  m_mesh = nullptr;
  m_node = nullptr;
  for (int axis = 0; axis < DIMS; ++axis) {
    for (int dir = 0; dir < ndirs; ++dir) {
      m_border[axis][dir].reset();
    }
  }
}

static std::string resid_name() {
  return "dgt::Block::m_resid";
}

static std::string soln_name(int i) {
  return "dgt::Block::m_soln[" + std::to_string(i) + "]";
}

static std::string flux_name(int i) {
  return "dgt::Block::m_flux[" + std::to_string(i) + "]";
}

void Block::allocate(int nsoln, int nmodal_eq, int nflux_eq) {
  CALI_CXX_MARK_FUNCTION;
  verify_basis(basis());
  m_soln.resize(nsoln);
  p3a::grid3 const cgrid = generalize(cell_grid());
  int const nmodes = basis().nmodes;
  int const nside_pts = num_pts(dim()-1, basis().p);
  int const ncells = cgrid.size();
  m_resid = View<double***>(resid_name(), ncells, nmodal_eq, nmodes);
  for (int soln = 0; soln < nsoln; ++soln) {
    m_soln[soln] = View<double***>(soln_name(soln), ncells, nmodal_eq, nmodes);
  }
  for (int axis = 0; axis < dim(); ++axis) {
    int const nsides = get_side_grid(cgrid, axis).size();
    m_flux[axis] = View<double***>(flux_name(axis), nsides, nside_pts, nflux_eq);
  }
  for (int axis = 0; axis < dim(); ++axis) {
    for (int dir = 0; dir < ndirs; ++dir) {
      m_border[axis][dir].allocate(nmodal_eq, nflux_eq);
    }
  }
  for (int field = 0; field < nfields(); ++field) {
    m_fields[field].allocate(cell_grid());
  }
}

void Block::deallocate() {
  CALI_CXX_MARK_FUNCTION;
  m_resid = View<double***>();
  m_soln.resize(0);
  m_fields.resize(0);
  for (int axis = 0; axis < DIMS; ++axis) {
    for (int dir = 0; dir < ndirs; ++dir) {
      m_border[axis][dir].deallocate();
    }
  }
}

}
