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
  spdlog::debug("dgitle: building -> {}", base);
  return b;
}

template Basis<View> build_basis(int, int, int, bool);
template Basis<HostView> build_basis(int, int, int, bool);

}
