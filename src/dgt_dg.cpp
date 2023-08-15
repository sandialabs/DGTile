#include <fmt/format.h>

#include <spdlog/spdlog.h>

#include "dgt_dg.hpp"
#include "dgt_for_each.hpp"
#include "dgt_view.hpp"

namespace dgt {

static std::string basis_name(
    int const dim,
    int const p,
    int const q,
    bool const tensor)
{
  return fmt::format("Basis[dim={}][p={}][q={}][tensor={}]", dim, p, q, tensor);
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

template <template <class> class ViewT>
Basis<ViewT> build_basis(
    int const dim,
    int const p,
    int const q,
    bool const tensor)
{
  using fmt::format;
  using namespace dgt::basis_locations;
  Basis<ViewT> b;
  std::string const name = basis_name(dim, p, q, tensor);
  spdlog::debug("dgt: building -> {}", name);
  b.cell_weights = ViewT<real*>(name + ".cell_weights", num_gauss_points(dim, q));
  b.face_weights = ViewT<real*>(name + ".face_weights", num_gauss_points(dim-1, q));
  b.modes[CELL].points = ViewT<real**>(name + ".modes[CELL].points", num_gauss_points(dim, q), dim);
  b.modes[VERTICES].points = ViewT<real**>(name + ".modes[VERTICES].points", num_vertices(dim), dim);



  HostView<real*> cell_weights = get_cell_weights(dim, q);
  HostView<real*> face_weights = get_face_weights(dim, q);
  HostView<real**> cell_points = get_cell_points(dim, q);
  Kokkos::deep_copy(b.cell_weights, cell_weights);
  Kokkos::deep_copy(b.face_weights, face_weights);
  Kokkos::deep_copy(b.modes[CELL].points, cell_points);



  for (size_t i = 0; i < b.modes[CELL].points.extent(0); ++i) {
    for (size_t j = 0; j < b.modes[CELL].points.extent(1); ++j) {
    std::cout << b.modes[CELL].points(i, j) << "\n";
    }
  }


  return b;


}

template Basis<View> build_basis(int, int, int, bool);
template Basis<HostView> build_basis(int, int, int, bool);

}
