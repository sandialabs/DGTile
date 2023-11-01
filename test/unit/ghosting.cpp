#include <dgt_for_each.hpp>
#include <dgt_ghosting.hpp>
#include <dgt_mesh.hpp>

#include <gtest/gtest.h>

#include <dgt_print.hpp> // debug

using namespace dgt;

static Mesh get_single_block_mesh(
    mpicpp::comm* comm,
    int const dim,
    int const p,
    int const q,
    bool const tensor)
{
  Grid3 cell_grid(0,0,0);
  Grid3 block_grid(0,0,0);
  Vec3<bool> periodic(false, false, false);
  Box3<real> domain({0,0,0}, {0,0,0});
  for (int axis = 0; axis < dim; ++axis) {
    block_grid.extents()[axis] = 1;
    cell_grid.extents()[axis] = 4;
    periodic[axis] = true;
    domain.upper()[axis] = 1.;
  }
  Mesh mesh;
  mesh.set_comm(comm);
  mesh.set_domain(domain);
  mesh.set_cell_grid(cell_grid);
  mesh.set_periodic(periodic);
  mesh.set_basis(p, q, tensor);
  mesh.add_modal({"hydro", 2, 5, true});
  mesh.initialize(block_grid);
  return mesh;
}

static void fill_U0(Mesh& mesh)
{
  Grid3 const cell_grid = mesh.cell_grid();
  auto U = (mesh.get_solution("hydro", 0)).get();
  auto functor = [=] DGT_DEVICE (
      int const block,
      Vec3<int> const& cell_ijk) DGT_ALWAYS_INLINE
  {
    int const cell = cell_grid.index(cell_ijk);

    // debug
    std::cout << cell_ijk << "\n";
    std::cout << cell << "\n";
    std::cout << "\n";

    for (int eq = 0; eq < 5; ++eq) {
      U[block](cell, eq, 0) = cell;
    }
  };
  for_each("fill_U0", mesh.owned_leaves().size(), cell_grid, functor);
}

TEST(ghosting, build_single_block_2D)
{
  mpicpp::comm comm = mpicpp::comm::world();
  Mesh mesh = get_single_block_mesh(&comm, 2, 1, 2, true);
  Ghosting ghosting;
  ghosting.build(mesh);
  fill_U0(mesh);
  auto const U = mesh.get_solution("hydro", 0);
  auto const B = mesh.basis();
  ghosting.begin_transfer(U, B, 0, 5);
}
