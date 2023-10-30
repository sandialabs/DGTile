#include <dgt_ghosting.hpp>
#include <dgt_mesh.hpp>

#include <gtest/gtest.h>

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
    cell_grid.extents()[axis] = 6;
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

TEST(ghosting, single_block)
{
  mpicpp::comm comm = mpicpp::comm::world();
  Mesh mesh = get_single_block_mesh(&comm, 2, 1, 2, true);
  Ghosting ghosting;
  ghosting.build(mesh);
}
