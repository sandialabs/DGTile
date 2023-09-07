#include <fmt/format.h>

#include "dgt_cartesian.hpp"
#include "dgt_dg.hpp"
#include "dgt_for_each.hpp"
#include "dgt_view.hpp"

namespace dgt {

static std::string basis_name(BasisInfo const& in)
{
  return fmt::format("Basis(dim={},p={},q={},tensor={})",
      in.dim, in.p, in.q, in.tensor);
}

static std::string loc_name(int const loc)
{
  std::string names[basis_locations::NUM] = {
    "CELL",
    "VERTICES",
    "XMIN_FACE",
    "XMAX_FACE",
    "YMIN_FACE",
    "YMAX_FACE",
    "ZMIN_FACE",
    "ZMAX_FACE",
    "EVALUATION"
  };
  return names[loc];
}

static HostView<real*> get_cell_weights(int const dim, int const q)
{
  HostView<real*> wts("", num_gauss_points(dim, q));
  Grid3 const bounds = tensor_bounds(dim, q-1);
  auto functor = [&] (Vec3<int> const& ijk) {
    int const pt = bounds.index(ijk);
    wts(pt)               = get_gauss_weight(q, ijk.x());
    if (dim > 1) wts(pt) *= get_gauss_weight(q, ijk.y());
    if (dim > 2) wts(pt) *= get_gauss_weight(q, ijk.z());
  };
  seq_for_each(bounds, functor);
  return wts;
}

static HostView<real*> get_face_weights(int const dim, int const q)
{
  HostView<real*> wts("", num_gauss_points(dim-1, q));
  Grid3 const bounds = tensor_bounds(dim-1, q-1);
  auto functor = [&] (Vec3<int> const& ijk) {
    int const pt = bounds.index(ijk);
    wts(pt) = 1.;
    if (dim > 1) wts(pt) *= get_gauss_weight(q, ijk.x());
    if (dim > 2) wts(pt) *= get_gauss_weight(q, ijk.y());
  };
  seq_for_each(bounds, functor);
  return wts;
}

static HostView<real**> get_cell_points(int const dim, int const q)
{
  HostView<real**> pts("", num_gauss_points(dim, q), dim);
  Grid3 const bounds = tensor_bounds(dim, q-1);
  auto functor = [&] (Vec3<int> const& ijk) {
    int const pt = bounds.index(ijk);
    for (int axis = 0; axis < dim; ++axis) {
      pts(pt, axis) = get_gauss_point(q, ijk[axis]);
    }
  };
  seq_for_each(bounds, functor);
  return pts;
}

static HostView<real**> get_vert_points(int const dim)
{
  HostView<real**> pts("", num_vertices(dim), dim);
  static double vertex_pt[2] = {-1., 1.};
  Grid3 const bounds = tensor_bounds(dim, 1);
  auto functor = [&] (Vec3<int> const& ijk) {
    int const pt = bounds.index(ijk);
    for (int axis = 0; axis < dim; ++axis) {
      pts(pt, axis) = vertex_pt[ijk[axis]];
    }
  };
  seq_for_each(bounds, functor);
  return pts;
}

static HostView<real**> get_face_points(
    int const dim,
    int const q,
    int const axis,
    int const dir)
{
  HostView<real**> pts("", num_gauss_points(dim-1, q), dim);
  Grid3 const bounds = tensor_bounds(dim-1, q-1);
  auto functor = [&] (Vec3<int> const& ijk) {
    int const pt = bounds.index(ijk);
    int const ii = (axis == X) ? Y : X;
    int const jj = (axis == Z) ? Y : Z;
    pts(pt, axis) = get_dir_sign(dir);
    if (dim > 1) pts(pt, ii) = get_gauss_point(q, ijk.x());
    if (dim > 2) pts(pt, jj) = get_gauss_point(q, ijk.y());
  };
  seq_for_each(bounds, functor);
  return pts;
}

static HostView<real**> get_eval_points(int const dim, int const q)
{
  int const npts = num_evaluation_points(dim, q);
  HostView<real**> pts("", npts, dim);
  int eval_pt = 0;
  HostView<real**> cell = get_cell_points(dim, q);
  for (size_t pt = 0; pt < cell.extent(0); ++pt) {
    for (int d = 0; d < dim; ++d) {
      pts(eval_pt, d) = cell(pt, d);
    }
    eval_pt++;
  }
  for (int axis = 0; axis < dim; ++axis) {
    for (int dir = 0; dir < DIRECTIONS; ++dir) {
      HostView<real**> face = get_face_points(dim, q, axis, dir);
      for (size_t pt = 0; pt < face.extent(0); ++pt) {
        for (int d = 0; d < dim; ++d) {
          pts(eval_pt, d) = face(pt, d);
        }
        eval_pt++;
      }
    }
  }
  return pts;
}

static HostView<real**> get_phis(
    int const dim,
    int const p,
    bool const tensor,
    HostView<real**> const pts)
{
  int const npts = pts.extent(0);
  int const nmodes = num_modes(dim, p, tensor);
  HostView<real**> phis("", npts, nmodes);
  Vec3<real> xi = Vec3<real>::zero();
  for (int pt = 0; pt < npts; ++pt) {
    for (int d = 0;  d < dim; ++d) {
      xi[d] = pts(pt, d);
    }
    ModeVector const phi_xi = get_modes(dim, p, tensor, xi);
    for (int m = 0; m < nmodes; ++m) {
      phis(pt, m) = phi_xi[m];
    }
  }
  return phis;
}

static HostView<real***> get_dphis(
    int const dim,
    int const p,
    bool const tensor,
    HostView<real**> const pts)
{
  int const npts = pts.extent(0);
  int const nmodes = num_modes(dim, p, tensor);
  HostView<real***> dphis("", dim, npts, nmodes);
  Vec3<real> xi = Vec3<real>::zero();
  for (int axis = 0; axis < dim; ++axis) {
    for (int pt = 0; pt < npts; ++pt) {
      for (int d = 0; d < dim; ++d) {
        xi[d] = pts(pt, d);
      }
      Vec3<int> const deriv = Vec3<int>::axis(axis);
      ModeVector const dphi_dxi = get_dmodes(dim, p, tensor, deriv, xi);
      for (int m = 0; m < nmodes; ++m) {
        dphis(axis, pt, m) = dphi_dxi[m];
      }
    }
  }
  return dphis;
}

static HostView<real*> get_mass(
    int const dim,
    int const p,
    bool const tensor)
{
  int const nmodes = num_modes(dim, p, tensor);
  int const ncell_pts = num_gauss_points(dim, p+1);
  HostView<real*> const wts = get_cell_weights(dim, p+1);
  HostView<real**> const pts = get_cell_points(dim, p+1);
  HostView<real**> phis = get_phis(dim, p, tensor, pts);
  HostView<real*> mass("", nmodes);
  for (int pt = 0; pt < ncell_pts; ++pt) {
    for (int m = 0; m < nmodes; ++m) {
      mass(m) += phis(pt, m) * phis(pt, m) * wts(pt);
    }
  }
  return mass;
}

template <template <class> class ViewT>
void build_mass(
    Basis<ViewT>& B,
    std::string const& name,
    int const dim,
    int const p,
    bool const tensor)
{
  int const nmodes = num_modes(dim, p, tensor);
  B.mass = ViewT<real*>(name + ".mass", nmodes);
  HostView<real*> mass = get_mass(dim, p, tensor);
  Kokkos::deep_copy(B.mass, mass);
}

template <template <class> class ViewT>
void build_weights(
    Basis<ViewT>& B,
    std::string const& name,
    int const dim,
    int const q)
{
  int const ncell_pts = num_gauss_points(dim, q);
  int const nface_pts = num_gauss_points(dim-1, q);
  B.cell_weights = ViewT<real*>(name + ".cell_weights", ncell_pts);
  B.face_weights = ViewT<real*>(name + ".face_weights", nface_pts);
  HostView<real*> cell_weights = get_cell_weights(dim, q);
  HostView<real*> face_weights = get_face_weights(dim, q);
  Kokkos::deep_copy(B.cell_weights, cell_weights);
  Kokkos::deep_copy(B.face_weights, face_weights);
}

template <template <class> class ViewT>
void build_mode(
    Basis<ViewT>& B,
    std::string const& name,
    int const location,
    int const dim,
    int const p,
    bool const tensor,
    HostView<real**> const pts)
{
  size_t const npts = pts.extent(0);
  size_t const nmodes = num_modes(dim, p, tensor);
  std::string const lname = "[" + loc_name(location) + "]";
  B.modes[location].points = ViewT<real**>(name + lname + ".points", npts, dim);
  B.modes[location].phis = ViewT<real**>(name + lname + ".phis", npts, nmodes);
  B.modes[location].grad_phis = ViewT<real***>(name + lname + ".grad_phis", dim, npts, nmodes);
  HostView<real**> phis = get_phis(dim, p, tensor, pts);
  HostView<real***> grad_phis = get_dphis(dim, p, tensor, pts);
  Kokkos::deep_copy(B.modes[location].points, pts);
  Kokkos::deep_copy(B.modes[location].phis, phis);
  Kokkos::deep_copy(B.modes[location].grad_phis, grad_phis);
}

template <template <class> class ViewT>
Basis<ViewT> build_basis(BasisInfo const& in)
{
  using namespace dgt::basis_locations;
  Basis<ViewT> B;
  std::string const name = basis_name(in);
  B.dim = in.dim;
  B.p = in.p;
  B.q = in.q;
  B.tensor = in.tensor;
  B.num_modes = num_modes(in.dim, in.p, in.tensor);
  B.num_cell_pts = num_gauss_points(in.dim, in.q);
  B.num_vert_pts = num_vertices(in.dim);
  B.num_face_pts = num_gauss_points(in.dim-1, in.q);
  B.num_eval_pts = num_evaluation_points(in.dim, in.q);
  HostView<real**> cell_pts = get_cell_points(in.dim, in.q);
  HostView<real**> vert_pts = get_vert_points(in.dim);
  HostView<real**> eval_pts = get_eval_points(in.dim, in.q);
  build_weights(B, name, in.dim, in.q);
  build_mass(B, name, in.dim, in.p, in.tensor);
  build_mode(B, name, CELL, in.dim, in.p, in.tensor, cell_pts);
  build_mode(B, name, VERTICES, in.dim, in.p, in.tensor, vert_pts);
  for (int axis = 0; axis < in.dim; ++axis) {
    for (int dir = 0; dir < DIRECTIONS; ++dir) {
      HostView<real**> pts = get_face_points(in.dim, in.q, axis, dir);
      build_mode(B, name, face(axis, dir), in.dim, in.p, in.tensor, pts);
    }
  }
  build_mode(B, name, EVALUATION, in.dim, in.p, in.tensor, eval_pts);
  return B;
}

template Basis<View> build_basis(BasisInfo const&);
template Basis<HostView> build_basis(BasisInfo const&);

}
