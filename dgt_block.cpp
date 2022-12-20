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

grid3 Block::cell_grid() const {
  verify_mesh(m_mesh);
  return m_mesh->cell_grid();
}

box3<double> Block::domain() const {
  verify_mesh(m_mesh);
  verify_node(m_node);
  Point const base_pt = m_mesh->tree().base();
  Point const node_pt = m_node->pt();
  box3<double> const box = m_mesh->domain();
  return get_block_domain(base_pt, node_pt, box);
}

vector3<double> Block::dx() const {
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

simd_view<double***> Block::simd_soln(int idx) const {
  verify_U_idx(m_soln.size(), idx);
  return simd_view<double***>(m_soln[idx]);
}

simd_view<double***> Block::simd_flux(int axis) const {
  verify_axis(dim(), axis);
  return simd_view<double***>(m_flux[axis]);
}

simd_view<double***> Block::simd_resid() const {
  return simd_view<double***>(m_resid);
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

void Block::allocate(int nsoln, int neq) {
  CALI_CXX_MARK_FUNCTION;
  using Kokkos::resize;
  verify_basis(basis());
  m_soln.resize(nsoln);
  grid3 const cgrid = generalize(cell_grid());
  int const nmodes = basis().nmodes;
  int const nside_pts = num_pts(dim()-1, basis().p);
  int const ncells = cgrid.size();
  resize(m_resid, ncells, neq, nmodes);
  for (int soln = 0; soln < nsoln; ++soln) {
    resize(m_soln[soln], ncells, neq, nmodes);
  }
  for (int axis = 0; axis < dim(); ++axis) {
    int const nsides = get_side_grid(cgrid, axis).size();
    resize(m_flux[axis], nsides, nside_pts, neq);
  }
  for (int axis = 0; axis < dim(); ++axis) {
    for (int dir = 0; dir < ndirs; ++dir) {
      m_border[axis][dir].allocate(neq);
    }
  }
}

void Block::deallocate() {
  CALI_CXX_MARK_FUNCTION;
  m_resid = View<double***>();
  m_soln.resize(0);
  for (int axis = 0; axis < DIMS; ++axis) {
    for (int dir = 0; dir < ndirs; ++dir) {
      m_border[axis][dir].deallocate();
    }
  }
}

}
