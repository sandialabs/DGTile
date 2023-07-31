#include "gtest/gtest.h"

#include "p3a_for_each.hpp"

#include "dgt_amr.hpp"
#include "dgt_grid.hpp"
#include "dgt_mesh.hpp"
#include "dgt_interp.hpp"
#include "dgt_spatial.hpp"

TEST(amr, local) {
  ASSERT_EQ(dgt::get_local(0), p3a::vector3<int>(0,0,0));
  ASSERT_EQ(dgt::get_local(1), p3a::vector3<int>(1,0,0));
  ASSERT_EQ(dgt::get_local(2), p3a::vector3<int>(0,1,0));
  ASSERT_EQ(dgt::get_local(3), p3a::vector3<int>(1,1,0));
  ASSERT_EQ(dgt::get_local(4), p3a::vector3<int>(0,0,1));
  ASSERT_EQ(dgt::get_local(5), p3a::vector3<int>(1,0,1));
  ASSERT_EQ(dgt::get_local(6), p3a::vector3<int>(0,1,1));
  ASSERT_EQ(dgt::get_local(7), p3a::vector3<int>(1,1,1));
}

TEST(amr, which_child) {
  ASSERT_EQ(dgt::get_which_child({0,0,0}), 0);
  ASSERT_EQ(dgt::get_which_child({1,0,0}), 1);
  ASSERT_EQ(dgt::get_which_child({0,1,0}), 2);
  ASSERT_EQ(dgt::get_which_child({1,1,0}), 3);
  ASSERT_EQ(dgt::get_which_child({0,0,1}), 4);
  ASSERT_EQ(dgt::get_which_child({1,0,1}), 5);
  ASSERT_EQ(dgt::get_which_child({0,1,1}), 6);
  ASSERT_EQ(dgt::get_which_child({1,1,1}), 7);
}

TEST(amr, axis_local) {
  ASSERT_EQ(dgt::get_local(dgt::X, 0), p3a::vector3<int>(0,0,0));
  ASSERT_EQ(dgt::get_local(dgt::X, 1), p3a::vector3<int>(0,1,0));
  ASSERT_EQ(dgt::get_local(dgt::X, 2), p3a::vector3<int>(0,0,1));
  ASSERT_EQ(dgt::get_local(dgt::X, 3), p3a::vector3<int>(0,1,1));
  ASSERT_EQ(dgt::get_local(dgt::Y, 0), p3a::vector3<int>(0,0,0));
  ASSERT_EQ(dgt::get_local(dgt::Y, 1), p3a::vector3<int>(1,0,0));
  ASSERT_EQ(dgt::get_local(dgt::Y, 2), p3a::vector3<int>(0,0,1));
  ASSERT_EQ(dgt::get_local(dgt::Y, 3), p3a::vector3<int>(1,0,1));
  ASSERT_EQ(dgt::get_local(dgt::Z, 0), p3a::vector3<int>(0,0,0));
  ASSERT_EQ(dgt::get_local(dgt::Z, 1), p3a::vector3<int>(1,0,0));
  ASSERT_EQ(dgt::get_local(dgt::Z, 2), p3a::vector3<int>(0,1,0));
  ASSERT_EQ(dgt::get_local(dgt::Z, 3), p3a::vector3<int>(1,1,0));
}

TEST(amr, axis_which_child) {
  ASSERT_EQ(dgt::get_which_child(dgt::X, {0,0,0}), 0);
  ASSERT_EQ(dgt::get_which_child(dgt::X, {0,1,0}), 1);
  ASSERT_EQ(dgt::get_which_child(dgt::X, {0,0,1}), 2);
  ASSERT_EQ(dgt::get_which_child(dgt::X, {0,1,1}), 3);
  ASSERT_EQ(dgt::get_which_child(dgt::Y, {0,0,0}), 0);
  ASSERT_EQ(dgt::get_which_child(dgt::Y, {1,0,0}), 1);
  ASSERT_EQ(dgt::get_which_child(dgt::Y, {0,0,1}), 2);
  ASSERT_EQ(dgt::get_which_child(dgt::Y, {1,0,1}), 3);
  ASSERT_EQ(dgt::get_which_child(dgt::Z, {0,0,0}), 0);
  ASSERT_EQ(dgt::get_which_child(dgt::Z, {1,0,0}), 1);
  ASSERT_EQ(dgt::get_which_child(dgt::Z, {0,1,0}), 2);
  ASSERT_EQ(dgt::get_which_child(dgt::Z, {1,1,0}), 3);
}

TEST(amr, map_to_fine) {
  ASSERT_EQ(dgt::map_to_fine({0,0,0}, {0,0,0}, {10,10,10}), p3a::vector3<int>(0,0,0));
  ASSERT_EQ(dgt::map_to_fine({0,0,0}, {1,0,0}, {10,10,10}), p3a::vector3<int>(10,0,0));
  ASSERT_EQ(dgt::map_to_fine({0,0,0}, {0,1,0}, {10,10,10}), p3a::vector3<int>(0,10,0));
  ASSERT_EQ(dgt::map_to_fine({0,0,0}, {1,1,0}, {10,10,10}), p3a::vector3<int>(10,10,0));
  ASSERT_EQ(dgt::map_to_fine({0,0,0}, {0,0,1}, {10,10,10}), p3a::vector3<int>(0,0,10));
  ASSERT_EQ(dgt::map_to_fine({0,0,0}, {1,0,1}, {10,10,10}), p3a::vector3<int>(10,0,10));
  ASSERT_EQ(dgt::map_to_fine({0,0,0}, {0,1,1}, {10,10,10}), p3a::vector3<int>(0,10,10));
  ASSERT_EQ(dgt::map_to_fine({0,0,0}, {1,1,1}, {10,10,10}), p3a::vector3<int>(10,10,10));
}

TEST(amr, local_from_fine_ijk) {
  ASSERT_EQ(dgt::get_local_from_fine_ijk({10,10,10}), p3a::vector3<int>(0,0,0));
  ASSERT_EQ(dgt::get_local_from_fine_ijk({11,10,10}), p3a::vector3<int>(1,0,0));
  ASSERT_EQ(dgt::get_local_from_fine_ijk({10,11,10}), p3a::vector3<int>(0,1,0));
  ASSERT_EQ(dgt::get_local_from_fine_ijk({11,11,10}), p3a::vector3<int>(1,1,0));
  ASSERT_EQ(dgt::get_local_from_fine_ijk({10,10,11}), p3a::vector3<int>(0,0,1));
  ASSERT_EQ(dgt::get_local_from_fine_ijk({11,10,11}), p3a::vector3<int>(1,0,1));
  ASSERT_EQ(dgt::get_local_from_fine_ijk({10,11,11}), p3a::vector3<int>(0,1,1));
  ASSERT_EQ(dgt::get_local_from_fine_ijk({11,11,11}), p3a::vector3<int>(1,1,1));
}

TEST(amr, coarse_ijk) {
  ASSERT_EQ(dgt::get_coarse_ijk({5,3,1}), p3a::vector3<int>(2,1,0));
  ASSERT_EQ(dgt::get_coarse_ijk({5,3,0}), p3a::vector3<int>(2,1,0));
  ASSERT_EQ(dgt::get_coarse_ijk({5,0,0}), p3a::vector3<int>(2,0,0));
}

static constexpr int neq = 1;
static constexpr int nsoln = 1;

P3A_ALWAYS_INLINE P3A_HOST_DEVICE inline
double exact_solution(p3a::vector3<double> const& x, int p) {
  return std::pow(x.x(), p) + std::pow(x.y(), p) + std::pow(x.z(), p);
}

static dgt::Mesh* create_mesh(mpicpp::comm* comm, int dim, int p, bool tensor) {
  dgt::Mesh* mesh = new dgt::Mesh;
  p3a::vector3<double> const xmin(0.,0.,0.);
  p3a::vector3<double> const xmax((dim>0)?1.:0., (dim>1)?1.:0., (dim>2)?1.:0.);
  p3a::vector3<int> const nblocks((dim>0)?1:0, (dim>1)?1:0, (dim>2)?1:0);
  p3a::vector3<int> const grid((dim>0)?4:0, (dim>1)?4:0, (dim>2)?4:0);
  mesh->set_comm(comm);
  mesh->set_domain({xmin, xmax});
  mesh->set_cell_grid(grid);
  mesh->set_nsoln(nsoln);
  mesh->set_nmodal_eq(neq);
  mesh->set_nflux_eq(neq);
  mesh->init(nblocks,  p, tensor);
  mesh->rebuild();
  return mesh;
}

static void set_mesh_data(dgt::Mesh* mesh) {
  for (dgt::Node* node : mesh->owned_leaves()) {
    dgt::Block& block = node->block;
    dgt::Basis const b = block.basis();
    int const nintr_pts = dgt::num_pts(b.dim, b.p);
    p3a::grid3 const cell_grid = dgt::generalize(block.cell_grid());
    p3a::vector3<double> const origin = block.domain().lower();
    p3a::vector3<double> const dx = block.dx();
    dgt::View<double***> U = block.soln(0);
    auto f = [=] P3A_DEVICE (p3a::vector3<int> const& cell_ijk) {
      int const cell = cell_grid.index(cell_ijk);
      for (int pt = 0; pt < nintr_pts; ++pt) {
        double const wt = b.wt_intr(pt);
        p3a::vector3<double> const xi = dgt::get_intr_pt(b, pt);
        p3a::vector3<double> const x = dgt::get_x(cell_ijk, origin, dx, xi);
        double const ex = exact_solution(x, b.p);
        for (int m = 0; m < b.nmodes; ++m) {
          double const phi = b.phi_intr(pt, m);
          double const mass = b.mass(m);
          U(cell, 0, m) += ex * phi * wt / mass;
        }
      }
    };
    p3a::for_each(p3a::execution::par, cell_grid, f);
  }
}

static void test_solution(dgt::Block const& block) {
  dgt::Basis const b = block.basis();
  p3a::vector3<double> const dx = block.dx();
  p3a::vector3<double> const origin = block.domain().lower();
  p3a::grid3 const cell_grid = dgt::generalize(block.cell_grid());
  int const nintr_pts = dgt::num_pts(b.dim, b.p);
  dgt::View<double***> U = block.soln(0);
  dgt::View<double**> interped;
  Kokkos::resize(interped, cell_grid.size(), nintr_pts);
  auto f = [=] P3A_DEVICE (p3a::vector3<int> const& cell_ijk) {
    int const cell = cell_grid.index(cell_ijk);
    for (int pt = 0; pt < nintr_pts; ++pt) {
      interped(cell, pt) = interp_scalar_intr(U, b, cell, pt, 0);
    }
  };
  p3a::for_each(p3a::execution::par, cell_grid, f);
  dgt::HView<double**> hinterped;
  dgt::HView<double**> hpt_intr;
  dgt::copy(interped, hinterped);
  dgt::copy(b.pt_intr, hpt_intr);
  p3a::execution::par.synchronize();
  auto g = [&] (p3a::vector3<int> const& cell_ijk) {
    int const cell = cell_grid.index(cell_ijk);
    for (int pt = 0; pt < nintr_pts; ++pt) {
      p3a::vector3<double> xi(0.,0.,0.);
      if (b.dim > 0) xi[dgt::X] = hpt_intr(pt, dgt::X);
      if (b.dim > 1) xi[dgt::Y] = hpt_intr(pt, dgt::Y);
      if (b.dim > 2) xi[dgt::Z] = hpt_intr(pt, dgt::Z);
      p3a::vector3<double> const x = dgt::get_x(cell_ijk, origin, dx, xi);
      double const ex = exact_solution(x, b.p);
      double const val = hinterped(cell, pt);
      ASSERT_NEAR(ex, val, 1.e-12);
    }
  };
  p3a::for_each(p3a::execution::seq, cell_grid, g);
}

static void test_prolong(
    p3a::vector3<int> const& local,
    dgt::Block const& parent,
    dgt::Block& child) {
  p3a::grid3 const from_grid = parent.cell_grid();
  p3a::grid3 const to_grid = child.cell_grid();
  p3a::vector3<int> const half = dgt::get_coarse_ijk(from_grid.extents());
  p3a::vector3<int> const from_start = p3a::hadamard_product(local, half);
  p3a::vector3<int> const from_end = from_start + half;
  p3a::subgrid3 const from_subgrid(from_start, from_end);
  p3a::subgrid3 const to_subgrid(from_grid);
  do_prolongation(
      parent.basis(),
      parent.soln(0), child.soln(0),
      from_grid, to_grid,
      from_subgrid, to_subgrid);
  test_solution(child);
}

static void test_restrict(
    p3a::vector3<int> const& local,
    dgt::Block const& child,
    dgt::Block& parent) {
  p3a::grid3 const from_grid = child.cell_grid();
  p3a::grid3 const to_grid = parent.cell_grid();
  p3a::vector3<int> const half = dgt::get_coarse_ijk(to_grid.extents());
  p3a::vector3<int> const to_start = p3a::hadamard_product(local, half);
  p3a::vector3<int> const to_end = to_start + half;
  p3a::subgrid3 const from_subgrid(from_grid);
  p3a::subgrid3 const to_subgrid(to_start, to_end);
  do_restriction(
      parent.basis(),
      child.soln(0), parent.soln(0),
      from_grid, to_grid,
      from_subgrid, to_subgrid);
}

static void test_prolong(int dim, int p, bool tensor) {
  mpicpp::comm comm = mpicpp::comm::world();
  dgt::Mesh* mesh = create_mesh(&comm, dim, p, tensor);
  mesh->allocate();
  mesh->rebuild();
  set_mesh_data(mesh);
  refine(dim, mesh->owned_leaves()[0]);
  mesh->rebuild();
  dgt::Node* parent = mesh->tree().find({0, {0,0,0}});
  auto f = [&] (p3a::vector3<int> const& local) {
    dgt::Node* child = parent->child(local);
    child->block.allocate(nsoln, neq, neq);
    test_prolong(local, parent->block, child->block);
  };
  p3a::for_each(p3a::execution::seq, dgt::generalize(dgt::get_child_grid(dim)), f);
  delete mesh;
}

static void test_restrict(int dim, int p, bool tensor) {
  mpicpp::comm comm = mpicpp::comm::world();
  dgt::Mesh* mesh = create_mesh(&comm, dim, p, tensor);
  refine(dim, mesh->owned_leaves()[0]);
  mesh->rebuild();
  mesh->allocate();
  set_mesh_data(mesh);
  dgt::Node* parent = mesh->tree().find({0, {0,0,0}});
  parent->block.allocate(nsoln, neq, neq);
  auto f = [&] (p3a::vector3<int> const& local) {
    dgt::Node* child = parent->child(local);
    test_restrict(local, child->block, parent->block);
  };
  p3a::for_each(p3a::execution::seq, dgt::generalize(dgt::get_child_grid(dim)), f);
  test_solution(parent->block);
}

TEST(prolong, d2_p0_tensor) {
  test_prolong(2, 0, true);
}

TEST(prolong, d2_p1_tensor) {
  test_prolong(2, 1, true);
}

TEST(prolong, d2_p2_tensor) {
  test_prolong(2, 2, true);
}

TEST(prolong, d3_p0_tensor) {
  test_prolong(3, 0, true);
}

TEST(prolong, d3_p1_tensor) {
  test_prolong(3, 1, true);
}

TEST(prolong, d3_p2_tensor) {
  test_prolong(3, 2, true);
}

TEST(restrict, d2_p0_tensor) {
  test_restrict(2, 0, true);
}

TEST(restrict, d2_p1_tensor) {
  test_restrict(2, 1, true);
}

TEST(restrict, d2_p2_tensor) {
  test_restrict(2, 2, true);
}

TEST(restrict, d3_p0_tensor) {
  test_restrict(3, 0, true);
}

TEST(restrict, d3_p1_tensor) {
  test_restrict(3, 1, true);
}

TEST(restrict, d3_p2_tensor) {
  test_restrict(3, 2, true);
}
