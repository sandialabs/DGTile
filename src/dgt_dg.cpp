#include <fmt/format.h>

#include <spdlog/spdlog.h>

#include "dgt_dg.hpp"
#include "dgt_for_each.hpp"
#include "dgt_view.hpp"

namespace dgt {

namespace names {

static std::string base(
    int const dim,
    int const p,
    int const q,
    bool const tensor)
{
  return fmt::format("Basis[dim={}][p={}][q={}][tensor={}]", dim, p, q, tensor);
}

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

template <template <class> class ViewT>
Basis<ViewT> build_basis(
    int const dim,
    int const p,
    int const q,
    bool const tensor)
{
  using fmt::format;
  Basis<ViewT> b;
  std::string const base = names::base(dim, p, q, tensor);
  spdlog::debug("dgt: building -> {}", base);

  get_cell_weights(dim, q);
  return b;
}

template Basis<View> build_basis(int, int, int, bool);
template Basis<HostView> build_basis(int, int, int, bool);

}
