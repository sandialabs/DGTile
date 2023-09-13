#include <dgt_dg.hpp>
#include <dgt_view.hpp>

#include <gtest/gtest.h>

using namespace dgt;

template <class BasisT>
void test_integers(
    BasisT const& B,
    int const dim,
    int const p,
    int const q,
    bool const tensor)
{
  EXPECT_EQ(B.dim, dim);
  EXPECT_EQ(B.p, p);
  EXPECT_EQ(B.q, q);
  EXPECT_EQ(B.tensor, tensor);
  EXPECT_EQ(B.num_modes, num_modes(dim, p, tensor));
  EXPECT_EQ(B.num_cell_pts, num_gauss_points(dim, q));
  EXPECT_EQ(B.num_face_pts, num_gauss_points(dim-1, q));
  EXPECT_EQ(B.num_vert_pts, num_vertices(dim));
  EXPECT_EQ(B.num_eval_pts, num_evaluation_points(dim, q));
}

template <class BasisT>
void test_view_sizes(BasisT const& B)
{
  using namespace basis_locations;
  EXPECT_EQ(B.mass.extent(0), B.num_modes);
  EXPECT_EQ(B.cell_weights.extent(0), B.num_cell_pts);
  EXPECT_EQ(B.face_weights.extent(0), B.num_face_pts);
  EXPECT_EQ(B.modes[CELL].points.extent(0), B.num_cell_pts);
  EXPECT_EQ(B.modes[CELL].points.extent(1), B.dim);
  EXPECT_EQ(B.modes[CELL].phis.extent(0), B.num_cell_pts);
  EXPECT_EQ(B.modes[CELL].phis.extent(1), B.num_modes);
  EXPECT_EQ(B.modes[CELL].grad_phis.extent(0), B.dim);
  EXPECT_EQ(B.modes[CELL].grad_phis.extent(1), B.num_cell_pts);
  EXPECT_EQ(B.modes[CELL].grad_phis.extent(2), B.num_modes);
  EXPECT_EQ(B.modes[VERTICES].points.extent(0), num_vertices(B.dim));
  EXPECT_EQ(B.modes[VERTICES].points.extent(1), B.dim);
  EXPECT_EQ(B.modes[VERTICES].phis.extent(0), num_vertices(B.dim));
  EXPECT_EQ(B.modes[VERTICES].phis.extent(1), B.num_modes);
  EXPECT_EQ(B.modes[VERTICES].grad_phis.extent(0), B.dim);
  EXPECT_EQ(B.modes[VERTICES].grad_phis.extent(1), num_vertices(B.dim));
  EXPECT_EQ(B.modes[VERTICES].grad_phis.extent(2), B.num_modes);
  EXPECT_EQ(B.modes[EVALUATION].points.extent(0), num_evaluation_points(B.dim, B.q));
  EXPECT_EQ(B.modes[EVALUATION].points.extent(1), B.dim);
  EXPECT_EQ(B.modes[EVALUATION].phis.extent(0), num_evaluation_points(B.dim, B.q));
  EXPECT_EQ(B.modes[EVALUATION].phis.extent(1), B.num_modes);
  EXPECT_EQ(B.modes[EVALUATION].grad_phis.extent(0), B.dim);
  EXPECT_EQ(B.modes[EVALUATION].grad_phis.extent(1), num_evaluation_points(B.dim, B.q));
  EXPECT_EQ(B.modes[EVALUATION].grad_phis.extent(2), B.num_modes);
}

void test_basis(
    int const dim,
    int const p,
    int const q,
    bool const tensor);
