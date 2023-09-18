#include <fmt/core.h>

#include "dgt_block_info.hpp"
#include "dgt_cartesian.hpp"

namespace dgt {

template <template <class> class ViewT, class T>
ViewT<T*> make_view(std::string const& name, std::vector<T> const& data)
{
  ViewT<T*> view(name, data.size());
  HostView<T*> host_view(name, data.size());
  for (std::size_t i = 0; i < data.size(); ++i) {
    host_view[i] = data[i];
  }
  Kokkos::deep_copy(view, host_view);
  return view;
}

static std::vector<std::int8_t> get_levels(
    int const dim,
    std::vector<tree::ID> const& ids)
{
  std::vector<std::int8_t> levels(ids.size());
  for (std::size_t i = 0; i < ids.size(); ++i) {
    tree::ID const id = ids[i];
    levels[i] = tree::get_level(dim, id);
  }
  return levels;
}

static Vec3<real> get_base_dx(
    int const dim,
    Grid3 const& cell_grid,
    Box3<real> const& domain,
    tree::Point const& base_pt)
{
  Vec3<real> dx = Vec3<real>::zero();
  Vec3<real> const length = domain.extents();
  for (int axis = 0; axis < dim; ++axis) {
    int const nblocks = base_pt.ijk[axis];
    int const ncells = cell_grid.extents()[axis] - 2;
    real const den = real(nblocks*ncells);
    dx[axis] = length[axis] / den;
  }
  return dx;
}

static Vec3<real> get_dx(
    int const dim,
    Grid3 const& cell_grid,
    Box3<real> const& domain,
    tree::Point const& base_pt,
    tree::ID const id)
{
  int const level = int(tree::get_level(dim, id));
  int const diff = level - base_pt.level;
  Vec3<real> const base_dx = get_base_dx(dim, cell_grid, domain, base_pt);
  return base_dx / std::pow(2., diff);
}

static std::vector<Box3<real>> get_domains(
    int const dim,
    Grid3 const& cell_grid,
    Box3<real> const& domain,
    std::vector<tree::ID> const& ids,
    tree::Point const& base_pt)
{
  std::vector<Box3<real>> domains(ids.size());
  for (std::size_t i = 0; i < ids.size(); ++i) {
    tree::ID const id = ids[i];
    Vec3<real> const dx = get_dx(dim, cell_grid, domain, base_pt, id);
    Box3<real> block_domain = tree::get_domain(dim, id, base_pt, domain);
    block_domain.lower() -= dx;
    block_domain.upper() += dx;
    domains[i] = block_domain;
  }
  return domains;
}

static std::vector<Vec3<real>> get_cell_dxs(
    int const dim,
    Grid3 const& cell_grid,
    Box3<real> const& domain,
    std::vector<tree::ID> const& ids,
    tree::Point const& base_pt)
{
  std::vector<Vec3<real>> cell_dxs(ids.size());
  for (std::size_t i = 0; i < ids.size(); ++i) {
    tree::ID const id = ids[i];
    cell_dxs[i] = get_dx(dim, cell_grid, domain, base_pt, id);
  }
  return cell_dxs;
}

static real get_volume(int const dim, Vec3<real> const& v)
{
  real volume = v.x();
  if (dim > 1) volume *= v.y();
  if (dim > 2) volume *= v.z();
  return volume;
}

static real get_cell_detJ(int const dim, Vec3<real> const& dx)
{
  return std::pow(0.5, dim) * get_volume(dim, dx);
}

static real get_face_detJ(int const dim, int const axis, Vec3<real> const& dx) {
  return std::pow(0.5, dim-1) * get_volume(dim, dx) / dx[axis];
}

static std::vector<real> get_cell_detJs(
    int const dim,
    Grid3 const& cell_grid,
    Box3<real> const& domain,
    std::vector<tree::ID> const& ids,
    tree::Point const& base_pt)
{
  std::vector<real> detJs(ids.size());
  for (std::size_t i = 0; i < ids.size(); ++i) {
    tree::ID const id = ids[i];
    Vec3<real> const dx = get_dx(dim, cell_grid, domain, base_pt, id);
    detJs[i] = get_cell_detJ(dim, dx);
  }
  return detJs;
}

static std::vector<real> get_face_detJs(
    int const dim,
    int const axis,
    Grid3 const& cell_grid,
    Box3<real> const& domain,
    std::vector<tree::ID> const& ids,
    tree::Point const& base_pt)
{
  std::vector<real> detJs(ids.size());
  for (std::size_t i = 0; i < ids.size(); ++i) {
    tree::ID const id = ids[i];
    Vec3<real> const dx = get_dx(dim, cell_grid, domain, base_pt, id);
    detJs[i] = get_face_detJ(dim, axis, dx);
  }
  return detJs;
}

template <template <class> class ViewT>
BlockInfo<ViewT> build_block_info(
  int const dim,
  Grid3 const& cell_grid,
  Box3<real> const& domain,
  tree::OwnedLeaves const& ids,
  tree::Point const& base_pt)
{
  BlockInfo<ViewT> B;
  auto const levels = get_levels(dim, ids);
  auto const domains = get_domains(dim, cell_grid, domain, ids, base_pt);
  auto const cell_dxs = get_cell_dxs(dim, cell_grid, domain, ids, base_pt);
  auto const cell_detJs = get_cell_detJs(dim, cell_grid, domain, ids, base_pt);
  B.ids = make_view<ViewT>("block_info.ids", ids);
  B.levels = make_view<ViewT>("block_info.levels", levels);
  B.domains = make_view<ViewT>("block_info.domains", domains);
  B.cell_dxs = make_view<ViewT>("block_info.cell_dxs", cell_dxs);
  B.cell_detJs = make_view<ViewT>("block_info.cell_detJs", cell_detJs);
  for (int axis = 0; axis < dim; ++axis) {
    auto const name = fmt::format("block_info.face_detJs[{}]", get_axis_name(axis));
    auto const face_detJs = get_face_detJs(dim, axis, cell_grid, domain, ids, base_pt);
    B.face_detJs[axis] = make_view<ViewT>(name, face_detJs);
  }
  return B;
}

template BlockInfo<View>
build_block_info(int const, Grid3 const&, Box3<real> const&, tree::OwnedLeaves const&, tree::Point const&);

template BlockInfo<HostView>
build_block_info(int const, Grid3 const&, Box3<real> const&, tree::OwnedLeaves const&, tree::Point const&);

}
