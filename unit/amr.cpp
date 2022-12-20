#include "gtest/gtest.h"

#include "p3a_for_each.hpp"

#include "dgt_amr.hpp"
#include "dgt_grid.hpp"
#include "dgt_mesh.hpp"
#include "dgt_interp.hpp"
#include "dgt_spatial.hpp"

using namespace dgt;

TEST(amr, local) {
  ASSERT_EQ(get_local(0), vector3<int>(0,0,0));
  ASSERT_EQ(get_local(1), vector3<int>(1,0,0));
  ASSERT_EQ(get_local(2), vector3<int>(0,1,0));
  ASSERT_EQ(get_local(3), vector3<int>(1,1,0));
  ASSERT_EQ(get_local(4), vector3<int>(0,0,1));
  ASSERT_EQ(get_local(5), vector3<int>(1,0,1));
  ASSERT_EQ(get_local(6), vector3<int>(0,1,1));
  ASSERT_EQ(get_local(7), vector3<int>(1,1,1));
}

TEST(amr, which_child) {
  ASSERT_EQ(get_which_child({0,0,0}), 0);
  ASSERT_EQ(get_which_child({1,0,0}), 1);
  ASSERT_EQ(get_which_child({0,1,0}), 2);
  ASSERT_EQ(get_which_child({1,1,0}), 3);
  ASSERT_EQ(get_which_child({0,0,1}), 4);
  ASSERT_EQ(get_which_child({1,0,1}), 5);
  ASSERT_EQ(get_which_child({0,1,1}), 6);
  ASSERT_EQ(get_which_child({1,1,1}), 7);
}

TEST(amr, axis_local) {
  ASSERT_EQ(get_local(X, 0), vector3<int>(0,0,0));
  ASSERT_EQ(get_local(X, 1), vector3<int>(0,1,0));
  ASSERT_EQ(get_local(X, 2), vector3<int>(0,0,1));
  ASSERT_EQ(get_local(X, 3), vector3<int>(0,1,1));
  ASSERT_EQ(get_local(Y, 0), vector3<int>(0,0,0));
  ASSERT_EQ(get_local(Y, 1), vector3<int>(1,0,0));
  ASSERT_EQ(get_local(Y, 2), vector3<int>(0,0,1));
  ASSERT_EQ(get_local(Y, 3), vector3<int>(1,0,1));
  ASSERT_EQ(get_local(Z, 0), vector3<int>(0,0,0));
  ASSERT_EQ(get_local(Z, 1), vector3<int>(1,0,0));
  ASSERT_EQ(get_local(Z, 2), vector3<int>(0,1,0));
  ASSERT_EQ(get_local(Z, 3), vector3<int>(1,1,0));
}

TEST(amr, axis_which_child) {
  ASSERT_EQ(get_which_child(X, {0,0,0}), 0);
  ASSERT_EQ(get_which_child(X, {0,1,0}), 1);
  ASSERT_EQ(get_which_child(X, {0,0,1}), 2);
  ASSERT_EQ(get_which_child(X, {0,1,1}), 3);
  ASSERT_EQ(get_which_child(Y, {0,0,0}), 0);
  ASSERT_EQ(get_which_child(Y, {1,0,0}), 1);
  ASSERT_EQ(get_which_child(Y, {0,0,1}), 2);
  ASSERT_EQ(get_which_child(Y, {1,0,1}), 3);
  ASSERT_EQ(get_which_child(Z, {0,0,0}), 0);
  ASSERT_EQ(get_which_child(Z, {1,0,0}), 1);
  ASSERT_EQ(get_which_child(Z, {0,1,0}), 2);
  ASSERT_EQ(get_which_child(Z, {1,1,0}), 3);
}

TEST(amr, map_to_fine) {
  ASSERT_EQ(map_to_fine({0,0,0}, {0,0,0}, {10,10,10}), vector3<int>(0,0,0));
  ASSERT_EQ(map_to_fine({0,0,0}, {1,0,0}, {10,10,10}), vector3<int>(10,0,0));
  ASSERT_EQ(map_to_fine({0,0,0}, {0,1,0}, {10,10,10}), vector3<int>(0,10,0));
  ASSERT_EQ(map_to_fine({0,0,0}, {1,1,0}, {10,10,10}), vector3<int>(10,10,0));
  ASSERT_EQ(map_to_fine({0,0,0}, {0,0,1}, {10,10,10}), vector3<int>(0,0,10));
  ASSERT_EQ(map_to_fine({0,0,0}, {1,0,1}, {10,10,10}), vector3<int>(10,0,10));
  ASSERT_EQ(map_to_fine({0,0,0}, {0,1,1}, {10,10,10}), vector3<int>(0,10,10));
  ASSERT_EQ(map_to_fine({0,0,0}, {1,1,1}, {10,10,10}), vector3<int>(10,10,10));
}

TEST(amr, local_from_fine_ijk) {
  ASSERT_EQ(get_local_from_fine_ijk({10,10,10}), vector3<int>(0,0,0));
  ASSERT_EQ(get_local_from_fine_ijk({11,10,10}), vector3<int>(1,0,0));
  ASSERT_EQ(get_local_from_fine_ijk({10,11,10}), vector3<int>(0,1,0));
  ASSERT_EQ(get_local_from_fine_ijk({11,11,10}), vector3<int>(1,1,0));
  ASSERT_EQ(get_local_from_fine_ijk({10,10,11}), vector3<int>(0,0,1));
  ASSERT_EQ(get_local_from_fine_ijk({11,10,11}), vector3<int>(1,0,1));
  ASSERT_EQ(get_local_from_fine_ijk({10,11,11}), vector3<int>(0,1,1));
  ASSERT_EQ(get_local_from_fine_ijk({11,11,11}), vector3<int>(1,1,1));
}

TEST(amr, coarse_ijk) {
  ASSERT_EQ(get_coarse_ijk({5,3,1}), vector3<int>(2,1,0));
  ASSERT_EQ(get_coarse_ijk({5,3,0}), vector3<int>(2,1,0));
  ASSERT_EQ(get_coarse_ijk({5,0,0}), vector3<int>(2,0,0));
}

static constexpr int neq = 1;
static constexpr int nsoln = 1;

P3A_ALWAYS_INLINE P3A_HOST_DEVICE inline
double exact_solution(vector3<double> const& x, int p) {
  return std::pow(x.x(), p) + std::pow(x.y(), p) + std::pow(x.z(), p);
}

static Mesh* create_mesh(mpicpp::comm* comm, int dim, int p, bool tensor) {
  Mesh* mesh = new Mesh;
  vector3<double> const xmin(0.,0.,0.);
  vector3<double> const xmax((dim>0)?1.:0., (dim>1)?1.:0., (dim>2)?1.:0.);
  vector3<int> const nblocks((dim>0)?1:0, (dim>1)?1:0, (dim>2)?1:0);
  vector3<int> const grid((dim>0)?4:0, (dim>1)?4:0, (dim>2)?4:0);
  mesh->set_comm(comm);
  mesh->set_domain({xmin, xmax});
  mesh->set_cell_grid(grid);
  mesh->set_nsoln(neq);
  mesh->set_neq(nsoln);
  mesh->init(nblocks,  p, tensor);
  mesh->rebuild();
  return mesh;
}

static void set_mesh_data(Mesh* mesh) {
  for (Node* node : mesh->owned_leaves()) {
    Block& block = node->block;
    Basis const b = block.basis();
    int const nintr_pts = num_pts(b.dim, b.p);
    grid3 const cell_grid = generalize(block.cell_grid());
    vector3<double> const origin = block.domain().lower();
    vector3<double> const dx = block.dx();
    View<double***> U = block.soln(0);
    auto f = [=] P3A_DEVICE (vector3<int> const& cell_ijk) {
      int const cell = cell_grid.index(cell_ijk);
      for (int pt = 0; pt < nintr_pts; ++pt) {
        double const wt = b.wt_intr(pt);
        vector3<double> const xi = get_intr_pt(b, pt);
        vector3<double> const x = get_x(cell_ijk, origin, dx, xi);
        double const ex = exact_solution(x, b.p);
        for (int m = 0; m < b.nmodes; ++m) {
          double const phi = b.phi_intr(pt, m);
          double const mass = b.mass(m);
          U(cell, 0, m) += ex * phi * wt / mass;
        }
      }
    };
    for_each(execution::par, cell_grid, f);
  }
}

static void test_solution(Block const& block) {
  using DView2D = Kokkos::DualView<double**, Kokkos::LayoutLeft>;
  Basis const b = block.basis();
  vector3<double> const dx = block.dx();
  vector3<double> const origin = block.domain().lower();
  grid3 const cell_grid = generalize(block.cell_grid());
  int const nintr_pts = num_pts(b.dim, b.p);
  View<double***> U = block.soln(0);
  DView2D interped;
  Kokkos::resize(interped, cell_grid.size(), nintr_pts);
  auto interped_d = interped.d_view;
  auto interped_h = interped.h_view;
  auto f = [=] P3A_DEVICE (vector3<int> const& cell_ijk) {
    int const cell = cell_grid.index(cell_ijk);
    for (int pt = 0; pt < nintr_pts; ++pt) {
      interped_d(cell, pt) = interp_scalar_intr(U, b, cell, pt, 0);
    }
  };
  for_each(execution::par, cell_grid, f);
  interped.modify<DView2D::execution_space>();
  interped.sync<DView2D::host_mirror_space>();
  p3a::execution::par.synchronize();
  auto g = [&] (vector3<int> const& cell_ijk) {
    int const cell = cell_grid.index(cell_ijk);
    for (int pt = 0; pt < nintr_pts; ++pt) {
      vector3<double> const xi = get_intr_pt(b, pt);
      vector3<double> const x = get_x(cell_ijk, origin, dx, xi);
      double const ex = exact_solution(x, b.p);
      double const val = interped_h(cell, pt);
      ASSERT_NEAR(ex, val, 1.e-12);
    }
  };
  for_each(execution::seq, cell_grid, g);
}

static void test_prolong(
    vector3<int> const& local,
    Block const& parent,
    Block& child) {
  grid3 const from_grid = parent.cell_grid();
  grid3 const to_grid = child.cell_grid();
  vector3<int> const half = get_coarse_ijk(from_grid.extents());
  vector3<int> const from_start = hadamard_product(local, half);
  vector3<int> const from_end = from_start + half;
  subgrid3 const from_subgrid(from_start, from_end);
  subgrid3 const to_subgrid(from_grid);
  do_prolongation(
      parent.basis(),
      parent.soln(0), child.soln(0),
      from_grid, to_grid,
      from_subgrid, to_subgrid);
  test_solution(child);
}

static void test_restrict(
    vector3<int> const& local,
    Block const& child,
    Block& parent) {
  grid3 const from_grid = child.cell_grid();
  grid3 const to_grid = parent.cell_grid();
  vector3<int> const half = get_coarse_ijk(to_grid.extents());
  vector3<int> const to_start = hadamard_product(local, half);
  vector3<int> const to_end = to_start + half;
  subgrid3 const from_subgrid(from_grid);
  subgrid3 const to_subgrid(to_start, to_end);
  do_restriction(
      parent.basis(),
      child.soln(0), parent.soln(0),
      from_grid, to_grid,
      from_subgrid, to_subgrid);
}

static void test_prolong(int dim, int p, bool tensor) {
  mpicpp::comm comm = mpicpp::comm::world();
  Mesh* mesh = create_mesh(&comm, dim, p, tensor);
  mesh->allocate();
  mesh->rebuild();
  set_mesh_data(mesh);
  refine(dim, mesh->owned_leaves()[0]);
  mesh->rebuild();
  Node* parent = mesh->tree().find({0, {0,0,0}});
  auto f = [&] (vector3<int> const& local) {
    Node* child = parent->child(local);
    child->block.allocate(nsoln, neq);
    test_prolong(local, parent->block, child->block);
  };
  for_each(execution::seq, generalize(get_child_grid(dim)), f);
  delete mesh;
}

static void test_restrict(int dim, int p, bool tensor) {
  mpicpp::comm comm = mpicpp::comm::world();
  Mesh* mesh = create_mesh(&comm, dim, p, tensor);
  refine(dim, mesh->owned_leaves()[0]);
  mesh->rebuild();
  mesh->allocate();
  set_mesh_data(mesh);
  Node* parent = mesh->tree().find({0, {0,0,0}});
  parent->block.allocate(nsoln, neq);
  auto f = [&] (vector3<int> const& local) {
    Node* child = parent->child(local);
    test_restrict(local, child->block, parent->block);
  };
  for_each(execution::seq, generalize(get_child_grid(dim)), f);
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
