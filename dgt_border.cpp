#include "caliper/cali.h"

#include "p3a_static_array.hpp"

#include "dgt_amr.hpp"
#include "dgt_basis.hpp"
#include "dgt_border.hpp"
#include "dgt_grid.hpp"
#include "dgt_interp.hpp"
#include "dgt_mesh.hpp"
#include "dgt_spatial.hpp"

namespace dgt {

static void verify_node(Node const* n) {
  if (!n) {
    throw std::runtime_error("Border - unset node");
  }
}

static void verify_axis(int axis) {
  if ((axis < 0) || (axis > DIMS)) {
    throw std::runtime_error("Border - invalid axis");
  }
}

static void verify_dir(int dir) {
  if ((dir != left) && (dir != right)) {
    throw std::runtime_error("Border - invalid dir");
  }
}

static void verify_type(int type) {
  if ((type < STANDARD) || (type > BOUNDARY)) {
    throw std::runtime_error("Border - invalid type");
  }
}

static void verify_border(Border const& border) {
  verify_node(border.node());
  verify_axis(border.axis());
  verify_dir(border.dir());
  verify_type(border.type());
}

static void verify_coarse_to_fine_transfer(Node const* adj) {
  if (!adj) {
    throw std::runtime_error("post_transfer - invalid adj");
  }
}

static void verify(
    int dim,
    int p,
    bool tensor,
    int neq,
    p3a::grid3 const& cell_grid,
    p3a::grid3 const& bside_grid,
    View<double***> U,
    p3a::static_array<View<double***>, ndirs> const& U_border,
    p3a::static_array<View<double**>, ndirs> const&  U_avg_border) {
  bool valid = true;
  if (U.extent(0) != size_t(cell_grid.size())) valid = false;
  if (U.extent(1) != size_t(neq)) valid = false;
  if (U.extent(2) != size_t(num_modes(dim,p,tensor))) valid = false;
  for (int msg_dir = 0; msg_dir < ndirs; ++msg_dir) {
    if (U_border[msg_dir].extent(0) != size_t(bside_grid.size())) valid = false;
    if (U_border[msg_dir].extent(1) != size_t(num_pts(dim-1,p))) valid = false;
    if (U_border[msg_dir].extent(2) != size_t(neq)) valid = false;
    if (U_avg_border[msg_dir].extent(0) != size_t(bside_grid.size())) valid = false;
    if (U_avg_border[msg_dir].extent(1) != size_t(neq)) valid = false;
  }
  if (!valid) {
    throw std::runtime_error("fill_border: inconsistent arrays");
  }
}

static void verify(
    int dim,
    int p,
    bool tensor,
    int neq,
    int nchild,
    p3a::grid3 const& cell_grid,
    p3a::grid3 const& bside_grid,
    View<double***> U,
    p3a::static_array<View<double****>, ndirs> const& U_border,
    p3a::static_array<View<double***>, ndirs> const& U_avg_border) {
  bool valid = true;
  if (U.extent(0) != size_t(cell_grid.size())) valid = false;
  if (U.extent(1) != size_t(neq)) valid = false;
  if (U.extent(2) != size_t(num_modes(dim,p,tensor))) valid = false;
  for (int msg_dir = 0; msg_dir < ndirs; ++msg_dir) {
    if (U_border[msg_dir].extent(0) != size_t(bside_grid.size())) valid = false;
    if (U_border[msg_dir].extent(1) != size_t(nchild)) valid = false;
    if (U_border[msg_dir].extent(2) != size_t(num_pts(dim-1,p))) valid = false;
    if (U_border[msg_dir].extent(3) != size_t(neq)) valid = false;
    if (U_avg_border[msg_dir].extent(0) != size_t(bside_grid.size())) valid = false;
    if (U_avg_border[msg_dir].extent(1) != size_t(nchild)) valid = false;
    if (U_avg_border[msg_dir].extent(2) != size_t(neq)) valid = false;
  }
  if (!valid) {
    throw std::runtime_error("fill_amr_border: inconsistent arrays");
  }
}

static void verify(
    int dim,
    int p,
    int neq,
    int nchild,
    p3a::grid3 const& bside_grid,
    View<double****> U,
    View<double***> U_avg,
    View<double***> U_buffer,
    View<double**> U_avg_buffer) {
  bool valid = true;
  if (U.extent(0) != size_t(bside_grid.size())) valid = false;
  if (U.extent(1) != size_t(nchild)) valid = false;
  if (U.extent(2) != size_t(num_pts(dim-1,p))) valid = false;
  if (U.extent(3) != size_t(neq)) valid = false;
  if (U_avg.extent(0) != size_t(bside_grid.size())) valid = false;
  if (U_avg.extent(1) != size_t(nchild)) valid = false;
  if (U_avg.extent(2) != size_t(neq)) valid = false;
  if (U_buffer.extent(0) != size_t(bside_grid.size())) valid = false;
  if (U_buffer.extent(1) != size_t(num_pts(dim-1,p))) valid = false;
  if (U_buffer.extent(2) != size_t(neq)) valid = false;
  if (U_avg_buffer.extent(0) != size_t(bside_grid.size())) valid = false;
  if (U_avg_buffer.extent(1) != size_t(neq)) valid = false;
  if (!valid) {
    throw std::runtime_error("fill_amr_bb: inconsistent arrays");
  }
}

int Border::axis() const {
  return m_axis;
}

int Border::dir() const {
  return m_dir;
}

int Border::type() const {
  return m_type;
}

Node const* Border::node() const {
  return m_node;
}

Node const* Border::adj() const {
  return m_adj;
}

Message<double***>& Border::soln(int msg_dir) {
  return m_soln[msg_dir];
}

Message<double**>& Border::avg_soln(int msg_dir) {
  return m_avg_soln[msg_dir];
}

AMRBorderData& Border::amr(int msg_dir) {
  return m_amr[msg_dir];
}

p3a::static_array<View<double***>, ndirs> Border::soln() const {
  verify_border(*this);
  p3a::static_array<View<double***>, ndirs> r;
  int const left_send_or_recv  = (m_dir == left) ? recv : send;
  int const right_send_or_recv = (m_dir == left) ? send : recv;
  r[left]  = (m_soln[left_send_or_recv]).val;
  r[right] = (m_soln[right_send_or_recv]).val;
  return r;
}

p3a::static_array<View<double****>, ndirs> Border::amr_soln() const {
  verify_border(*this);
  p3a::static_array<View<double****>, ndirs> r;
  int const left_send_or_recv  = (m_dir == left) ? recv : send;
  int const right_send_or_recv = (m_dir == left) ? send : recv;
  r[left]   = m_amr[left_send_or_recv].soln;
  r[right]  = m_amr[right_send_or_recv].soln;
  return r;
}

View<double****> Border::amr_flux() const {
  return m_amr_flux;
}

void Border::set_axis(int axis) {
  m_axis = axis;
}

void Border::set_type(int type) {
  m_type = type;
}

void Border::set_dir(int dir) {
  m_dir = dir;
}

void Border::set_node(Node* node) {
  m_node = node;
}

void Border::set_adj(Node* adj) {
  m_adj = adj;
}

void Border::reset() {
  m_axis = -1;
  m_dir = -1;
  m_type = -1;
  m_node = nullptr;
  m_adj = nullptr;
}

static std::string name(int i) {
  return "dgt::Border::m_soln[" + std::to_string(i) + "].val";
}

static std::string avg_name(int i) {
  return "dgt::Border::m_avg_soln[" + std::to_string(i) + "].val";
}

static std::string amr_name(int i) {
  return "dgt::Border::m_amr[" + std::to_string(i) + "].soln";
}

static std::string amr_avg_name(int i) {
  return "dgt::Border::m_amr[" + std::to_string(i) + "].avg_soln";
}

static std::string amr_name(int i, int c) {
  return "dgt::Border::m_amr[" +
    std::to_string(i) + "].child_soln[" +
    std::to_string(c) + "]";
}

static std::string amr_avg_name(int i, int c) {
  return "dgt::Border::m_amr[" +
    std::to_string(i) + "].child_avg_soln[" +
    std::to_string(c) + "]";
}

static std::string amr_flux_name() {
  return "dgt::Border::m_amr_flux";
}

void Border::allocate(int nmodal_eq, int nflux_eq) {
  CALI_CXX_MARK_FUNCTION;
  verify_border(*this);
  int const dim = m_node->block.dim();
  int const p = m_node->block.basis().p;
  int const nchild = num_child(dim-1);
  int const npts = num_pts(dim-1, p);
  p3a::grid3 const g = m_node->block.cell_grid();
  p3a::subgrid3 const sides = generalize(get_adj_sides(g, m_axis, m_dir));
  int const nsides = sides.size();
  if (m_type == COARSE_TO_FINE) {
    m_amr_flux = View<double****>(
        amr_flux_name(), nsides, nchild, npts, nflux_eq);
  }
  for (int msg_dir = 0; msg_dir < ndirs; ++msg_dir) {
    if (m_type == COARSE_TO_FINE) {
      m_amr[msg_dir].soln = View<double****>(
          amr_name(msg_dir), nsides, nchild, npts, nmodal_eq);
      m_amr[msg_dir].avg_soln = View<double***>(
          amr_avg_name(msg_dir), nsides, nchild, nmodal_eq);
      for (int which_child = 0; which_child < nchild; ++which_child) {
        m_amr[msg_dir].child_soln[which_child].val = View<double***>(
            amr_name(msg_dir, which_child), nsides, npts, nmodal_eq);
        m_amr[msg_dir].child_avg_soln[which_child].val = View<double**>(
            amr_avg_name(msg_dir, which_child), nsides, nmodal_eq);
      }
    } else {
      m_soln[msg_dir].val = View<double***>(
          name(msg_dir), nsides, npts, nmodal_eq);
      m_avg_soln[msg_dir].val = View<double**>(
          avg_name(msg_dir), nsides, nmodal_eq);
    }
  }
}

void Border::deallocate() {
  CALI_CXX_MARK_FUNCTION;
  m_amr_flux = View<double****>();
  for (int msg_dir = 0; msg_dir < ndirs; ++msg_dir) {
    m_soln[msg_dir].val = View<double***>();
    m_avg_soln[msg_dir].val = View<double**>();
    m_amr[msg_dir].soln = View<double****>();
    m_amr[msg_dir].avg_soln = View<double***>();
    for (int which_child = 0; which_child < NBORDER_CHILD; ++which_child) {
      m_amr[msg_dir].child_soln[which_child].val = View<double***>();
      m_amr[msg_dir].child_avg_soln[which_child].val = View<double**>();
    }
  }
}

static p3a::static_array<View<double***>, ndirs> get_U(Border& border) {
  p3a::static_array<View<double***>, ndirs> r;
  r[send] = border.soln(send).val;
  r[recv] = border.soln(recv).val;
  return r;
}

static p3a::static_array<View<double**>, ndirs> get_U_avg(Border& border) {
  p3a::static_array<View<double**>, ndirs> r;
  r[send] = border.avg_soln(send).val;
  r[recv] = border.avg_soln(recv).val;
  return r;
}

static p3a::static_array<View<double****>, ndirs> get_amr_U(Border& border) {
  p3a::static_array<View<double****>, ndirs> r;
  r[send] = border.amr(send).soln;
  r[recv] = border.amr(recv).soln;
  return r;
}

static p3a::static_array<View<double***>, ndirs> get_amr_U_avg(Border& border) {
  p3a::static_array<View<double***>, ndirs> r;
  r[send] = border.amr(send).avg_soln;
  r[recv] = border.amr(recv).avg_soln;
  return r;
}

static void fill_border(Border& border, int soln_idx) {
  CALI_CXX_MARK_FUNCTION;
  if (border.type() == COARSE_TO_FINE) return;
  verify_border(border);
  Block const& block = border.node()->block;
  int const dim = block.dim();
  int const axis = border.axis();
  int const dir = border.dir();
  int const idir = invert_dir(dir);
  int const p = block.basis().p;
  bool const tensor = block.basis().tensor;
  int const npts = num_pts(dim-1, p);
  int const neq = block.soln(0).extent(1);
  Basis const b = block.basis();
  p3a::grid3 const g = block.cell_grid();
  p3a::grid3 const cell_grid = generalize(g);
  p3a::subgrid3 const sides = generalize(get_adj_sides(g, axis, dir));
  p3a::grid3 const border_side_grid(sides.extents());
  View<double***> U = block.soln(soln_idx);
  p3a::static_array<View<double***>, ndirs> U_border = get_U(border);
  p3a::static_array<View<double**>, ndirs> U_avg_border = get_U_avg(border);
  verify(dim, p, tensor, neq, cell_grid, border_side_grid, U, U_border, U_avg_border);

  auto f = [=] P3A_DEVICE (p3a::vector3<int> const& side_ijk, int const& eq, int const& pt) {
    p3a::vector3<int> const border_side_ijk = get_border_ijk(side_ijk, axis);
    p3a::vector3<int> const cell_ijk = get_sides_adj_cell(side_ijk, axis, idir);
    int const cell = cell_grid.index(cell_ijk);
    int const border_side = border_side_grid.index(border_side_ijk);
    if ( pt == 0 ) {
      auto const avg = U(cell, eq, 0);
      for (int msg_dir = 0; msg_dir < ndirs; ++msg_dir)
        U_avg_border[msg_dir](border_side, eq) = avg;
    }
    double const val = interp_scalar_side(U, b, cell, axis, dir, pt, eq);
    for (int msg_dir = 0; msg_dir < ndirs; ++msg_dir)
      U_border[msg_dir](border_side, pt, eq) = val;
  };
  p3a::for_each(p3a::execution::par, sides, neq, npts, f);
}

static void fill_amr_border(Border& border, int soln_idx) {
  CALI_CXX_MARK_FUNCTION;
  if (border.type() != COARSE_TO_FINE) return;
  verify_border(border);
  Block const& block = border.node()->block;
  int const dim = block.dim();
  int const axis = border.axis();
  int const dir = border.dir();
  int const idir = invert_dir(dir);
  int const p = block.basis().p;
  bool const tensor = block.basis().tensor;
  int const npts = num_pts(dim-1, p);
  int const neq = block.soln(0).extent(1);
  int const nchild = num_child(dim-1);
  Basis const b = block.basis();
  p3a::grid3 const g = block.cell_grid();
  p3a::grid3 const cell_grid = generalize(g);
  p3a::subgrid3 const sides = generalize(get_adj_sides(g, axis, dir));
  p3a::grid3 const border_side_grid(sides.extents());
  View<double***> U = block.soln(soln_idx);
  p3a::static_array<View<double****>, ndirs> U_border = get_amr_U(border);
  p3a::static_array<View<double***>, ndirs> U_avg_border = get_amr_U_avg(border);
  verify(dim, p, tensor, neq, nchild, cell_grid, border_side_grid, U, U_border, U_avg_border);

  auto f = [=] P3A_DEVICE (p3a::vector3<int> const& side_ijk, int const& eq, int const& pt) {
    p3a::vector3<int> const border_side_ijk = get_border_ijk(side_ijk, axis);
    p3a::vector3<int> const cell_ijk = get_sides_adj_cell(side_ijk, axis, idir);
    int const cell = cell_grid.index(cell_ijk);
    int const border_side = border_side_grid.index(border_side_ijk);
    for (int which_child = 0; which_child < nchild; ++which_child) {
      if ( pt == 0 ) {
        double const avg = U(cell, eq, 0);
        for (int msg_dir = 0; msg_dir < ndirs; ++msg_dir)
          U_avg_border[msg_dir](border_side, which_child, eq) = avg;
      }
      double const val = interp_scalar_child_side(U, b, cell, axis, dir, which_child, pt, eq);
      for (int msg_dir = 0; msg_dir < ndirs; ++msg_dir)
        U_border[msg_dir](border_side, which_child, pt, eq) = val;
    }
  };
  p3a::for_each(p3a::execution::par, sides, neq, npts, f);
}

static void fill_amr_buffer_from_border(
    Border& border,
    int border_which_child) {
  CALI_CXX_MARK_FUNCTION;
  Block const& block = border.node()->block;
  int const dim = block.dim();
  int const axis = border.axis();
  int const dir = border.dir();
  int const p = block.basis().p;
  int const npts = num_pts(dim-1, p);
  int const neq = block.soln(0).extent(1);
  int const nchild = num_child(dim-1);
  p3a::grid3 const g = block.cell_grid();
  p3a::subgrid3 const sides = generalize(get_adj_sides(g, axis, dir));
  p3a::grid3 const border_side_grid(sides.extents());
  p3a::vector3<int> const border_local = get_local(axis, border_which_child);
  View<double****> U = border.amr(send).soln;
  View<double***> U_avg = border.amr(send).avg_soln;
  View<double***> U_buffer = border.amr(send).child_soln[border_which_child].val;
  View<double**> U_avg_buffer = border.amr(send).child_avg_soln[border_which_child].val;
  verify(dim, p, neq, nchild, border_side_grid, U, U_avg, U_buffer, U_avg_buffer);

  auto f = [=] P3A_DEVICE (p3a::vector3<int> const& ijk, int const& pt, int const& eq) {
    p3a::vector3<int> const fine_ijk = map_to_fine(ijk, border_local, border_side_grid);
    p3a::vector3<int> const coarse_ijk = get_coarse_ijk(fine_ijk);
    p3a::vector3<int> const local = get_local_from_fine_ijk(fine_ijk);
    int const which_child = get_which_child(axis, local);
    int const border_side = border_side_grid.index(ijk);
    int const border_side_coarse = border_side_grid.index(coarse_ijk);
    if ( pt == 0 )
      U_avg_buffer(border_side, eq) = U_avg(border_side_coarse, which_child, eq);
    U_buffer(border_side, pt, eq) = U(border_side_coarse, which_child, pt, eq);
  };
  p3a::for_each(p3a::execution::par, border_side_grid, npts, neq, f);
}

static void fill_amr_buffers_from_border(Border& border) {
  CALI_CXX_MARK_FUNCTION;
  if (border.type() != COARSE_TO_FINE) return;
  verify_border(border);
  int const dim = border.node()->block.dim();
  int const nchild = num_child(dim-1);
  for (int which_child = 0; which_child < nchild; ++which_child) {
    fill_amr_buffer_from_border(border, which_child);
  }
}

static void fill_amr_border_from_buffer(
    Border& border,
    int border_which_child) {
  CALI_CXX_MARK_FUNCTION;
  Block const& block = border.node()->block;
  int const dim = block.dim();
  int const axis = border.axis();
  int const dir = border.dir();
  int const p = block.basis().p;
  int const npts = num_pts(dim-1, p);
  int const neq = block.soln(0).extent(1);
  int const nchild = num_child(dim-1);
  p3a::grid3 const g = block.cell_grid();
  p3a::subgrid3 const sides = generalize(get_adj_sides(g, axis, dir));
  p3a::grid3 const border_side_grid(sides.extents());
  p3a::vector3<int> const border_local = get_local(axis, border_which_child);
  View<double****> U = border.amr(recv).soln;
  View<double***> U_avg = border.amr(recv).avg_soln;
  View<double***> U_buffer = border.amr(recv).child_soln[border_which_child].val;
  View<double**> U_avg_buffer = border.amr(recv).child_avg_soln[border_which_child].val;
  verify(dim, p, neq, nchild, border_side_grid, U, U_avg, U_buffer, U_avg_buffer);

  auto f = [=] P3A_DEVICE (p3a::vector3<int> const& ijk, int const& pt, int const& eq) {
    p3a::vector3<int> const fine_ijk = map_to_fine(ijk, border_local, border_side_grid);
    p3a::vector3<int> const coarse_ijk = get_coarse_ijk(fine_ijk);
    p3a::vector3<int> const local = get_local_from_fine_ijk(fine_ijk);
    int const which_child = get_which_child(axis, local);
    int const border_side = border_side_grid.index(ijk);
    int const border_side_coarse = border_side_grid.index(coarse_ijk);
    if ( pt == 0 )
      U_avg(border_side_coarse, which_child, eq) = U_avg_buffer(border_side, eq);
    U(border_side_coarse, which_child, pt, eq) = U_buffer(border_side, pt, eq);
  };
  p3a::for_each(p3a::execution::par, border_side_grid, npts, neq, f);
}

static void fill_amr_border_from_buffers(Border& border) {
  CALI_CXX_MARK_FUNCTION;
  if (border.type() != COARSE_TO_FINE) return;
  verify_border(border);
  int const dim = border.node()->block.dim();
  int const nchild = num_child(dim-1);
  for (int which_child = 0; which_child < nchild; ++which_child) {
    fill_amr_border_from_buffer(border, which_child);
  }
}

static int get_tag(int block_id, int axis, int dir, int which_child, int data) {
  static constexpr int ndata = 2;
  static constexpr int nborder = DIMS * ndirs;
  static constexpr int nchild = NBORDER_CHILD;
  int const border = axis * ndirs + dir;
  return ((block_id * nborder + border) * nchild + which_child) * ndata + data;
}

template <class T>
void post_transfer(
    Border& border,
    Block const& adj,
    int data,
    int which_child,
    Message<T>& send_msg,
    Message<T>& recv_msg) {
  CALI_CXX_MARK_FUNCTION;
  Block const& block = border.node()->block;
  mpicpp::comm* comm = block.mesh()->comm();
  int const axis = border.axis();
  int const dir = border.dir();
  int const idir = invert_dir(dir);
  int const owner = adj.owner();
  int const send_tag = get_tag(block.id(), axis, dir, which_child, data);
  int const recv_tag = get_tag(adj.id(), axis, idir, which_child, data);
  send_msg.send(comm, owner, send_tag);
  recv_msg.recv(comm, owner, recv_tag);
}

static void post_transfer(Border& border) {
  CALI_CXX_MARK_FUNCTION;
  int const type = border.type();
  if (type == BOUNDARY) {
    return;
  } else if (type == STANDARD) {
    Block const& adj = border.adj()->block;
    post_transfer(border, adj, 0, 0,
        border.soln(send),
        border.soln(recv));
    post_transfer(border, adj, 1, 0,
        border.avg_soln(send),
        border.avg_soln(recv));
  } else if (type == FINE_TO_COARSE) {
    Block const& adj = border.adj()->block;
    int const axis = border.axis();
    p3a::vector3<int> const ijk = border.node()->pt().ijk;
    p3a::vector3<int> const local = get_local_from_fine_ijk(ijk);
    int const which_child = get_which_child(axis, local);
    post_transfer(border, adj, 0, which_child,
        border.soln(send),
        border.soln(recv));
    post_transfer(border, adj, 1, which_child,
        border.avg_soln(send),
        border.avg_soln(recv));
  } else if (type == COARSE_TO_FINE) {
    Block const& block = border.node()->block;
    int const dim = block.dim();
    int const axis = border.axis();
    int const dir = border.dir();
    for (int which_child = 0; which_child < num_child(dim-1); ++which_child) {
      p3a::vector3<int> local = get_local(axis, which_child);
      local[axis] = invert_dir(dir);
      Node const* adj_node = border.adj()->child(local);
      verify_coarse_to_fine_transfer(adj_node);
      Block const& adj = adj_node->block;
      post_transfer(border, adj, 0, which_child,
          border.amr(send).child_soln[which_child],
          border.amr(recv).child_soln[which_child]);
      post_transfer(border, adj, 1, which_child,
          border.amr(send).child_avg_soln[which_child],
          border.amr(recv).child_avg_soln[which_child]);
    }
  }
}

void begin_border_transfer(Mesh& mesh, int soln_idx) {
  CALI_CXX_MARK_FUNCTION;
  for (Node* leaf : mesh.owned_leaves()) {
    for (int axis = 0; axis < mesh.dim(); ++axis) {
      for (int dir = 0; dir < ndirs; ++dir) {
        Border& border = leaf->block.border(axis, dir);
        fill_border(border, soln_idx);
        fill_amr_border(border, soln_idx);
        fill_amr_buffers_from_border(border);
        p3a::execution::par.synchronize();
        post_transfer(border);
      }
    }
  }
}

static void end_transfer(Border& border) {
  CALI_CXX_MARK_FUNCTION;
  border.soln(send).wait();
  border.soln(recv).wait();
  border.avg_soln(send).wait();
  border.avg_soln(recv).wait();
  for (int which_child = 0; which_child < NBORDER_CHILD; ++which_child) {
    border.amr(send).child_soln[which_child].wait();
    border.amr(recv).child_soln[which_child].wait();
    border.amr(send).child_avg_soln[which_child].wait();
    border.amr(recv).child_avg_soln[which_child].wait();
  }
}

void end_border_transfer(Mesh& mesh) {
  CALI_CXX_MARK_FUNCTION;
  for (Node* leaf : mesh.owned_leaves()) {
    for (int axis = 0; axis < mesh.dim(); ++axis) {
      for (int dir = 0; dir < ndirs; ++dir) {
        Border& border = leaf->block.border(axis, dir);
        end_transfer(border);
        fill_amr_border_from_buffers(border);
      }
    }
  }
}

}
