#include <dgt_block_info.hpp>

#include <gtest/gtest.h>

#include <dgt_print.hpp> // debug

using namespace dgt;

template <class ViewT>
void compare(ViewT view, std::vector<typename ViewT::value_type> const& expected)
{
  HostView<typename ViewT::value_type*> host_view("", view.size());
  Kokkos::deep_copy(host_view, view);
  EXPECT_EQ(view.size(), expected.size());
  for (std::size_t i = 0; i < view.size(); ++i) {
    EXPECT_EQ(host_view[i], expected[i]);
  }
}

TEST(block_info, build_1D_host)
{
  int const dim = 1;
  Box3<real> const domain({0.,0.,0.}, {1.,0.,0.});
  tree::OwnedLeaves const leaves{3,4,2};
  tree::Point const base_pt(1, {2,0,0});
  Grid3 const ghost_cell_grid(3,0,0);
  BlockInfo<HostView> B = build_block_info<HostView>(
    dim, ghost_cell_grid, domain, leaves, base_pt);
  compare(B.ids, {3,4,2});
  compare(B.levels, {2,2,1});
  compare(B.domains, {
      {{-0.25, 0., 0.}, {0.5,  0., 0.}},
      {{ 0.,   0., 0.}, {0.75, 0., 0.}},
      {{ 0.,   0., 0.}, {1.5,  0., 0.}} });
  compare(B.cell_dxs, {
      {0.25, 0., 0.},
      {0.25, 0., 0.},
      {0.5,  0., 0.} });
  compare(B.cell_detJs, {0.125, 0.125, 0.25});
  compare(B.face_detJs[X], {1., 1., 1.});
}

TEST(block_info, build_1D_device)
{
  int const dim = 1;
  Box3<real> const domain({0.,0.,0.}, {1.,0.,0.});
  tree::OwnedLeaves const leaves{3,4,2};
  tree::Point const base_pt(1, {2,0,0});
  Grid3 const ghost_cell_grid(3,0,0);
  BlockInfo<View> B = build_block_info<View>(
    dim, ghost_cell_grid, domain, leaves, base_pt);
  compare(B.ids, {3,4,2});
  compare(B.levels, {2,2,1});
  compare(B.domains, {
      {{-0.25, 0., 0.}, {0.5,  0., 0.}},
      {{ 0.,   0., 0.}, {0.75, 0., 0.}},
      {{ 0.,   0., 0.}, {1.5,  0., 0.}} });
  compare(B.cell_dxs, {
      {0.25, 0., 0.},
      {0.25, 0., 0.},
      {0.5,  0., 0.} });
  compare(B.cell_detJs, {0.125, 0.125, 0.25});
  compare(B.face_detJs[X], {1., 1., 1.});
}

TEST(block_info, build_2D_host)
{
  int const dim = 2;
  Box3<real> const domain({0.,0.,0.}, {1.,1.,0.});
  tree::OwnedLeaves const leaves{4,3,2,9,6,10,5};
  tree::Point const base_pt(1, {2,2,0});
  Grid3 const ghost_cell_grid(3,3,0);
  BlockInfo<HostView> B = build_block_info<HostView>(
      dim, ghost_cell_grid, domain, leaves, base_pt);
  compare(B.ids, {4,3,2,9,6,10,5});
  compare(B.levels, {1,1,1,2,2,2,2});
  compare(B.domains, {
      {{0.,     0.,   0.}, {1.5,  1.5,  0.}},
      {{-0.5,   0.,   0.}, {1.,   1.5,  0.}},
      {{0.,    -0.5,  0.}, {1.5,  1.,   0.}},
      {{-0.25,  0.,   0.}, {0.5,  0.75, 0.}},
      {{0.,    -0.25, 0.}, {0.75, 0.5,  0.}},
      {{0.,     0.,   0.}, {0.75, 0.75, 0.}},
      {{-0.25, -0.25, 0.}, {0.5,  0.5,  0.}} });
  compare(B.cell_dxs, {
      {0.5,  0.5,  0.},
      {0.5,  0.5,  0.},
      {0.5,  0.5,  0.},
      {0.25, 0.25, 0.},
      {0.25, 0.25, 0.},
      {0.25, 0.25, 0.},
      {0.25, 0.25, 0.} });
  compare(B.cell_detJs, {
      0.0625,
      0.0625,
      0.0625,
      0.015625,
      0.015625,
      0.015625,
      0.015625 });
  compare(B.face_detJs[X], {
      0.25,
      0.25,
      0.25,
      0.125,
      0.125,
      0.125,
      0.125 });
  compare(B.face_detJs[Y], {
      0.25,
      0.25,
      0.25,
      0.125,
      0.125,
      0.125,
      0.125 });
}

TEST(block_info, build_2D_device)
{
  int const dim = 2;
  Box3<real> const domain({0.,0.,0.}, {1.,1.,0.});
  tree::OwnedLeaves const leaves{4,3,2,9,6,10,5};
  tree::Point const base_pt(1, {2,2,0});
  Grid3 const ghost_cell_grid(3,3,0);
  BlockInfo<View> B = build_block_info<View>(
      dim, ghost_cell_grid, domain, leaves, base_pt);
  compare(B.ids, {4,3,2,9,6,10,5});
  compare(B.levels, {1,1,1,2,2,2,2});
  compare(B.domains, {
      {{0.,     0.,   0.}, {1.5,  1.5,  0.}},
      {{-0.5,   0.,   0.}, {1.,   1.5,  0.}},
      {{0.,    -0.5,  0.}, {1.5,  1.,   0.}},
      {{-0.25,  0.,   0.}, {0.5,  0.75, 0.}},
      {{0.,    -0.25, 0.}, {0.75, 0.5,  0.}},
      {{0.,     0.,   0.}, {0.75, 0.75, 0.}},
      {{-0.25, -0.25, 0.}, {0.5,  0.5,  0.}} });
  compare(B.cell_dxs, {
      {0.5,  0.5,  0.},
      {0.5,  0.5,  0.},
      {0.5,  0.5,  0.},
      {0.25, 0.25, 0.},
      {0.25, 0.25, 0.},
      {0.25, 0.25, 0.},
      {0.25, 0.25, 0.} });
  compare(B.cell_detJs, {
      0.0625,
      0.0625,
      0.0625,
      0.015625,
      0.015625,
      0.015625,
      0.015625 });
  compare(B.face_detJs[X], {
      0.25,
      0.25,
      0.25,
      0.125,
      0.125,
      0.125,
      0.125 });
  compare(B.face_detJs[Y], {
      0.25,
      0.25,
      0.25,
      0.125,
      0.125,
      0.125,
      0.125 });
}

#if 0
TEST(block_info, build_3D_host)
{
  int const dim = 3;
  Box3<real> const domain({0.,0.,0.}, {1.,1.,1.});
  tree::Point const base_pt(1, {2,2,2});
  tree::OwnedLeaves const leaves{9,10,13,14,25,26,29,30,2,3,4,5,6,7,8};
  BlockInfo<HostView> B = build_block_info<HostView>(dim, domain, leaves, base_pt);
  compare(B.global_ids, {9,10,13,14,25,26,29,30,2,3,4,5,6,7,8});
  compare(B.levels, {2,2,2,2,2,2,2,2,1,1,1,1,1,1,1});
  compare(B.domains, {
      {{0, 0, 0}, {0.25, 0.25, 0.25}},
      {{0.25, 0, 0}, {0.5, 0.25, 0.25}},
      {{0, 0.25, 0}, {0.25, 0.5, 0.25}},
      {{0.25, 0.25, 0}, {0.5, 0.5, 0.25}},
      {{0, 0, 0.25}, {0.25, 0.25, 0.5}},
      {{0.25, 0, 0.25}, {0.5, 0.25, 0.5}},
      {{0, 0.25, 0.25}, {0.25, 0.5, 0.5}},
      {{0.25, 0.25, 0.25}, {0.5, 0.5, 0.5}},
      {{0.5, 0, 0}, {1, 0.5, 0.5}},
      {{0, 0.5, 0}, {0.5, 1, 0.5}},
      {{0.5, 0.5, 0}, {1, 1, 0.5}},
      {{0, 0, 0.5}, {0.5, 0.5, 1}},
      {{0.5, 0, 0.5}, {1, 0.5, 1}},
      {{0, 0.5, 0.5}, {0.5, 1, 1}},
      {{0.5, 0.5, 0.5}, {1, 1, 1}} });
  compare(B.dxs, {
      {0.25, 0.25, 0.25},
      {0.25, 0.25, 0.25},
      {0.25, 0.25, 0.25},
      {0.25, 0.25, 0.25},
      {0.25, 0.25, 0.25},
      {0.25, 0.25, 0.25},
      {0.25, 0.25, 0.25},
      {0.25, 0.25, 0.25},
      {0.5,  0.5,  0.5},
      {0.5,  0.5,  0.5},
      {0.5,  0.5,  0.5},
      {0.5,  0.5,  0.5},
      {0.5,  0.5,  0.5},
      {0.5,  0.5,  0.5},
      {0.5,  0.5,  0.5} });
  compare(B.cell_detJs, {
      0.001953125,
      0.001953125,
      0.001953125,
      0.001953125,
      0.001953125,
      0.001953125,
      0.001953125,
      0.001953125,
      0.015625,
      0.015625,
      0.015625,
      0.015625,
      0.015625,
      0.015625,
      0.015625 });
  std::vector<real> face_detJs{
      0.015625,
      0.015625,
      0.015625,
      0.015625,
      0.015625,
      0.015625,
      0.015625,
      0.015625,
      0.0625,
      0.0625,
      0.0625,
      0.0625,
      0.0625,
      0.0625,
      0.0625};
  compare(B.face_detJs[X], face_detJs);
  compare(B.face_detJs[Y], face_detJs);
  compare(B.face_detJs[Z], face_detJs);
}

TEST(block_info, build_3D_device)
{
  int const dim = 3;
  Box3<real> const domain({0.,0.,0.}, {1.,1.,1.});
  tree::Point const base_pt(1, {2,2,2});
  tree::OwnedLeaves const leaves{9,10,13,14,25,26,29,30,2,3,4,5,6,7,8};
  BlockInfo<View> B = build_block_info<View>(dim, domain, leaves, base_pt);
  compare(B.global_ids, {9,10,13,14,25,26,29,30,2,3,4,5,6,7,8});
  compare(B.levels, {2,2,2,2,2,2,2,2,1,1,1,1,1,1,1});
  compare(B.domains, {
      {{0, 0, 0}, {0.25, 0.25, 0.25}},
      {{0.25, 0, 0}, {0.5, 0.25, 0.25}},
      {{0, 0.25, 0}, {0.25, 0.5, 0.25}},
      {{0.25, 0.25, 0}, {0.5, 0.5, 0.25}},
      {{0, 0, 0.25}, {0.25, 0.25, 0.5}},
      {{0.25, 0, 0.25}, {0.5, 0.25, 0.5}},
      {{0, 0.25, 0.25}, {0.25, 0.5, 0.5}},
      {{0.25, 0.25, 0.25}, {0.5, 0.5, 0.5}},
      {{0.5, 0, 0}, {1, 0.5, 0.5}},
      {{0, 0.5, 0}, {0.5, 1, 0.5}},
      {{0.5, 0.5, 0}, {1, 1, 0.5}},
      {{0, 0, 0.5}, {0.5, 0.5, 1}},
      {{0.5, 0, 0.5}, {1, 0.5, 1}},
      {{0, 0.5, 0.5}, {0.5, 1, 1}},
      {{0.5, 0.5, 0.5}, {1, 1, 1}} });
  compare(B.dxs, {
      {0.25, 0.25, 0.25},
      {0.25, 0.25, 0.25},
      {0.25, 0.25, 0.25},
      {0.25, 0.25, 0.25},
      {0.25, 0.25, 0.25},
      {0.25, 0.25, 0.25},
      {0.25, 0.25, 0.25},
      {0.25, 0.25, 0.25},
      {0.5,  0.5,  0.5},
      {0.5,  0.5,  0.5},
      {0.5,  0.5,  0.5},
      {0.5,  0.5,  0.5},
      {0.5,  0.5,  0.5},
      {0.5,  0.5,  0.5},
      {0.5,  0.5,  0.5} });
  compare(B.cell_detJs, {
      0.001953125,
      0.001953125,
      0.001953125,
      0.001953125,
      0.001953125,
      0.001953125,
      0.001953125,
      0.001953125,
      0.015625,
      0.015625,
      0.015625,
      0.015625,
      0.015625,
      0.015625,
      0.015625 });
  std::vector<real> face_detJs{
      0.015625,
      0.015625,
      0.015625,
      0.015625,
      0.015625,
      0.015625,
      0.015625,
      0.015625,
      0.0625,
      0.0625,
      0.0625,
      0.0625,
      0.0625,
      0.0625,
      0.0625};
  compare(B.face_detJs[X], face_detJs);
  compare(B.face_detJs[Y], face_detJs);
  compare(B.face_detJs[Z], face_detJs);
}
#endif
