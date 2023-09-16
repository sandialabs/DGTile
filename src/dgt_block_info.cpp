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

static std::vector<tree::ID> get_global_ids(std::vector<tree::ID> const& ids)
{
  std::vector<tree::ID> global_ids(ids.size());
  for (std::size_t i = 0; i < ids.size(); ++i) {
    global_ids[i] = ids[i];
  }
  return global_ids;
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

static std::vector<Box3<real>> get_domains(
    int const dim,
    Box3<real> const& domain,
    std::vector<tree::ID> const& ids,
    tree::Point const& base_pt)
{
  std::vector<Box3<real>> domains(ids.size());
  for (std::size_t i = 0; i < ids.size(); ++i) {
    tree::ID const id = ids[i];
    domains[i] = tree::get_domain(dim, id, base_pt, domain);
  }
  return domains;
}

static std::vector<Vec3<real>> get_dxs(
    int const dim,
    Box3<real> const& domain,
    std::vector<tree::ID> const& ids,
    tree::Point const& base_pt)
{
  std::vector<Vec3<real>> dxs(ids.size());
  for (std::size_t i = 0; i < ids.size(); ++i) {
    tree::ID const id = ids[i];
    Box3<real> const leaf_domain = tree::get_domain(dim, id, base_pt, domain);
    dxs[i] = leaf_domain.extents();
  }
  return dxs;
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
    Box3<real> const& domain,
    std::vector<tree::ID> const& ids,
    tree::Point const& base_pt)
{
  std::vector<real> detJs(ids.size());
  for (std::size_t i = 0; i < ids.size(); ++i) {
    tree::ID const id = ids[i];
    Box3<real> const leaf_domain = tree::get_domain(dim, id, base_pt, domain);
    Vec3<real> const dx = leaf_domain.extents();
    detJs[i] = get_cell_detJ(dim, dx);
  }
  return detJs;
}

static std::vector<real> get_face_detJs(
    int const dim,
    int const axis,
    Box3<real> const& domain,
    std::vector<tree::ID> const& ids,
    tree::Point const& base_pt)
{
  std::vector<real> detJs(ids.size());
  for (std::size_t i = 0; i < ids.size(); ++i) {
    tree::ID const id = ids[i];
    Box3<real> const leaf_domain = tree::get_domain(dim, id, base_pt, domain);
    Vec3<real> const dx = leaf_domain.extents();
    detJs[i] = get_face_detJ(dim, axis, dx);
  }
  return detJs;
}

template <template <class> class ViewT>
BlockInfo<ViewT> build_block_info(
  int const dim,
  Box3<real> const& domain,
  tree::OwnedLeaves const& ids,
  tree::Point const& base_pt)
{
  BlockInfo<ViewT> B;
  B.global_ids = make_view<ViewT>(
      "global_ids", get_global_ids(ids));
  B.levels = make_view<ViewT>(
      "levels", get_levels(dim, ids));
  B.domains = make_view<ViewT>(
      "domains",get_domains(dim, domain, ids, base_pt));
  B.dxs = make_view<ViewT>(
      "dxs", get_dxs(dim, domain, ids, base_pt));
  B.cell_detJs = make_view<ViewT>(
      "cell_detJs", get_cell_detJs(dim, domain, ids, base_pt));
  for (int axis = 0; axis < dim; ++axis) {
    auto const name = "face_detJs" + get_axis_name(axis);
    B.face_detJs[axis] = make_view<ViewT>(
        name, get_face_detJs(dim, axis, domain, ids, base_pt));
  }
  return B;
}

template BlockInfo<View>
build_block_info(int const, Box3<real> const&, tree::OwnedLeaves const&, tree::Point const&);

template BlockInfo<HostView>
build_block_info(int const, Box3<real> const&, tree::OwnedLeaves const&, tree::Point const&);

}
